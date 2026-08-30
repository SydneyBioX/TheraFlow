import os
import sys
import argparse
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import selfies as sf
from rdkit import Chem

# Paths setup
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
WORKSPACE_ROOT = os.path.dirname(os.path.dirname(BASE_DIR))
DATA_DIR = os.path.join(BASE_DIR, "data")
REST_METHODS_DIR = os.path.join(DATA_DIR, "test_datasets/rest_methods")
OUTPUTS_DIR = os.path.join(BASE_DIR, "outputs/fame")
REIMPL_ROOT = os.path.join(WORKSPACE_ROOT, "reimplementation")

DATASETS = [
    "lincs1000_1",
    "lincs1000_2",
    "cardiomyopathy",
    "parkinsons",
    "hcc",
    "glaucoma",
    "glaucoma_bulk"
]

def main():
    parser = argparse.ArgumentParser(description="Orchestrator to run FAME on test datasets")
    parser.add_argument("--dataset", choices=DATASETS, help="Run on a specific dataset")
    parser.add_argument("--run-all", action="store_true", help="Run all datasets sequentially")
    
    args = parser.parse_args()
    if not (args.dataset or args.run_all):
        parser.print_help()
        sys.exit(1)
        
    datasets = DATASETS if args.run_all else [args.dataset]
    
    # Change working directory to reimplementation directory context for relative loading
    original_cwd = os.getcwd()
    os.chdir(REIMPL_ROOT)
    
    print("Loading cyclemol cache and checkpoints...")
    cache = torch.load("cyclemol_cache.pt", map_location="cpu")
    weights = torch.load("checkpoints/fame.pt", map_location="cpu")["model"]
    vocab = cache["vocab"]
    norm = cache["normalization"]
    
    BOS_ID = cache["config"]["bos_id"]
    EOS_ID = cache["config"]["eos_id"]
    PAD_ID = cache["config"]["pad_id"]
    
    default_cell = int(torch.bincount(cache["cell"]).argmax())
    default_time = int(torch.bincount(cache["time"]).argmax())
    default_dose = float(cache["dose"].median())
    
    # Restore CWD
    os.chdir(original_cwd)
    
    # Create outputs directory
    os.makedirs(OUTPUTS_DIR, exist_ok=True)
    
    # Model components for manual forward pass
    def gru_layer(x, h, layer):
        wi = weights[f"decoder.gru.weight_ih_l{layer}"]
        wh = weights[f"decoder.gru.weight_hh_l{layer}"]
        bi = weights[f"decoder.gru.bias_ih_l{layer}"]
        bh = weights[f"decoder.gru.bias_hh_l{layer}"]
        
        gates_i = F.linear(x, wi, bi)
        gates_h = F.linear(h, wh, bh)
        
        r_i, z_i, n_i = gates_i.chunk(3, dim=-1)
        r_h, z_h, n_h = gates_h.chunk(3, dim=-1)
        
        r = torch.sigmoid(r_i + r_h)
        z = torch.sigmoid(z_i + z_h)
        n = torch.tanh(n_i + r * n_h)
        return (1 - z) * n + z * h

    def encode_profile(y, dose_norm, cell_id, time_id):
        # 1. Gene expression stream
        h = F.layer_norm(y, normalized_shape=[978], weight=weights["profile.inp.0.weight"], bias=weights["profile.inp.0.bias"])
        h = F.linear(h, weights["profile.inp.1.weight"], weights["profile.inp.1.bias"])
        h = F.gelu(h)
        
        for i in range(4):
            pre = F.layer_norm(h, normalized_shape=[384], weight=weights[f"profile.blocks.{i}.net.0.weight"], bias=weights[f"profile.blocks.{i}.net.0.bias"])
            mid = F.linear(pre, weights[f"profile.blocks.{i}.net.1.weight"], weights[f"profile.blocks.{i}.net.1.bias"])
            mid = F.gelu(mid)
            out = F.linear(mid, weights[f"profile.blocks.{i}.net.4.weight"], weights[f"profile.blocks.{i}.net.4.bias"])
            h = h + out
            
        h = F.layer_norm(h, normalized_shape=[384], weight=weights["profile.out.weight"], bias=weights["profile.out.bias"])
        
        # 2. Dose embedding
        d = torch.tensor([[dose_norm]], dtype=torch.float32)
        d = F.linear(d, weights["dose.0.weight"], weights["dose.0.bias"])
        d = F.gelu(d)
        d = F.linear(d, weights["dose.2.weight"], weights["dose.2.bias"])
        d = F.layer_norm(d, normalized_shape=[384], weight=weights["dose.3.weight"], bias=weights["dose.3.bias"])
        
        # 3. Cell & Time embeddings
        c = weights["cell.weight"][cell_id].unsqueeze(0)
        t = weights["time.weight"][time_id].unsqueeze(0)
        
        # 4. Fuse
        fuse_in = torch.cat([h, d, c, t, h, d], dim=-1)
        fuse_out = F.linear(fuse_in, weights["fuse.0.weight"], weights["fuse.0.bias"])
        fuse_out = F.gelu(fuse_out)
        fuse_out = F.linear(fuse_out, weights["fuse.3.weight"], weights["fuse.3.bias"])
        cond = F.layer_norm(fuse_out, normalized_shape=[384], weight=weights["fuse.4.weight"], bias=weights["fuse.4.bias"])
        return cond

    def sample_z(cond, fixed=False):
        p = F.linear(cond, weights["prior.0.weight"], weights["prior.0.bias"])
        p = F.gelu(p)
        p = F.linear(p, weights["prior.2.weight"], weights["prior.2.bias"])
        mu, logvar = p.chunk(2, dim=-1)
        if fixed:
            return mu
        std = torch.exp(0.5 * logvar.clamp(-4, 4))
        return mu + std * torch.randn_like(std)

    def get_h0(cond, z):
        h_in = torch.cat([cond, z], dim=-1)
        h0 = F.linear(h_in, weights["h0.0.weight"], weights["h0.0.bias"])
        h0 = torch.tanh(h0)
        return h0.view(4, 1, 384)

    @torch.no_grad()
    def generate_one(cond, z, max_len=128, temperature=1.0, top_k=80, top_p=0.97):
        h0 = get_h0(cond, z)
        hidden = [h0[i] for i in range(4)]
        token_ids = []
        curr_input = torch.tensor([[BOS_ID]])
        
        for _ in range(max_len):
            x = F.embedding(curr_input.squeeze(), weights["decoder.token_emb.weight"]).unsqueeze(0)
            for layer in range(4):
                hidden[layer] = gru_layer(x, hidden[layer], layer)
                x = hidden[layer]
            x_norm = F.layer_norm(x, normalized_shape=[384], weight=weights["decoder.norm.weight"], bias=weights["decoder.norm.bias"])
            logits = F.linear(x_norm, weights["decoder.out.weight"])
            
            logits = logits / max(temperature, 1e-8)
            if top_k > 0:
                thresh = torch.topk(logits, top_k).values[0, -1]
                logits[logits < thresh] = -1e9
            if top_p < 1.0:
                sorted_logits, sorted_idx = torch.sort(logits, descending=True)
                cum_probs = torch.cumsum(F.softmax(sorted_logits, dim=-1), dim=-1)
                remove = cum_probs - F.softmax(sorted_logits, dim=-1) > top_p
                sorted_logits[remove] = -1e9
                logits = torch.zeros_like(logits).scatter_(1, sorted_idx, sorted_logits)
            
            next_id = torch.multinomial(F.softmax(logits, dim=-1), 1).item()
            if next_id == EOS_ID:
                break
            token_ids.append(next_id)
            curr_input = torch.tensor([[next_id]])
        return token_ids

    def ids_to_smiles(token_ids):
        selfies_str = "".join(vocab["itos"][t] for t in token_ids if t not in (BOS_ID, EOS_ID, PAD_ID))
        try:
            smiles = sf.decoder(selfies_str)
            if Chem.MolFromSmiles(smiles):
                return smiles
        except Exception:
            pass
        return None

    def prepare_input(csv_path):
        df = pd.read_csv(csv_path, index_col=0)
        mean = torch.tensor(norm["ctp_mean"], dtype=torch.float32).numpy()
        std = torch.tensor(norm["ctp_std"], dtype=torch.float32).numpy()
        std = np.where(std < 1e-8, 1.0, std)
        X = (df.values.astype(np.float32) - mean) / std
        return torch.tensor(X, dtype=torch.float32)

    for d in datasets:
        src_path = os.path.join(REST_METHODS_DIR, f"{d}.csv")
        if not os.path.exists(src_path):
            print(f"Error: Preprocessed file {src_path} does not exist. Run preprocessing first.")
            continue
            
        print(f"\n[{d}] Running FAME inference...")
        y_all = prepare_input(src_path)
        
        all_results = []
        num_samples = 100
        temperature = 1.0
        top_k = 80
        top_p = 0.97
        
        for row_idx in range(len(y_all)):
            y = y_all[row_idx].unsqueeze(0)
            cond = encode_profile(y, default_dose, default_cell, default_time)
            
            valid, seen = [], set()
            for _ in range(num_samples):
                z = sample_z(cond)
                tok_ids = generate_one(cond, z, temperature=temperature, top_k=top_k, top_p=top_p)
                smi = ids_to_smiles(tok_ids)
                if smi and smi not in seen:
                    seen.add(smi)
                    valid.append(smi)
            all_results.extend(valid)
            
        df_out = pd.DataFrame(all_results, columns=['SMILES'])
        dest_path = os.path.join(OUTPUTS_DIR, f"{d}.csv")
        df_out.to_csv(dest_path, index=False)
        print(f"[{d}] Saved {len(df_out)} unique valid SMILES to {dest_path}")

if __name__ == "__main__":
    main()
