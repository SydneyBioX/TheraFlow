import os
import sys
import argparse
import math
import warnings
from pathlib import Path
import pandas as pd
import numpy as np
import torch
import torch.nn.functional as F
import selfies as sf
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore")

# Paths setup
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
WORKSPACE_ROOT = os.path.dirname(os.path.dirname(BASE_DIR))
DATA_DIR = os.path.join(BASE_DIR, "data")
REST_METHODS_DIR = os.path.join(DATA_DIR, "test_datasets/rest_methods")
OUTPUTS_DIR = os.path.join(BASE_DIR, "outputs/gex2sgen")
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

def load_and_align_expression(csv_path: str, cache: dict, gene_info_path: str = None) -> torch.Tensor:
    df = pd.read_csv(csv_path)
    df = df.select_dtypes(include=[np.number])

    if gene_info_path:
        gene_info = pd.read_csv(gene_info_path, sep="\t")
        lm = gene_info[gene_info["pr_is_lm"].astype(int) == 1]
        landmark_entrez = lm["pr_gene_id"].astype(str).tolist()

        df.columns = df.columns.astype(str)
        overlap = len(set(df.columns) & set(landmark_entrez))
        print(f"  Gene overlap: {overlap}/978")

        df = df.reindex(columns=landmark_entrez, fill_value=0.0)
        n_zeros = int((df == 0.0).all(axis=0).sum())
        if n_zeros > 0:
            print(f"  Warning: {n_zeros} genes missing — filled with 0")
    else:
        if df.shape[1] != 978:
            raise ValueError(
                f"CSV has {df.shape[1]} numeric columns but expected 978. "
                "Pass gene_info_path for automatic Entrez ID alignment."
            )
        print("  No gene_info_path — assuming columns already in landmark order")

    ctp_mean = torch.tensor(cache["normalization"]["ctp_mean"], dtype=torch.float32).numpy()
    ctp_std = torch.tensor(cache["normalization"]["ctp_std"], dtype=torch.float32).numpy()
    ctp_std = np.where(ctp_std < 1e-8, 1.0, ctp_std)

    X = (df.values.astype(np.float32) - ctp_mean) / ctp_std
    print(f"  Expression shape: {X.shape}  mean={X.mean():.3f}  std={X.std():.3f}")
    return torch.tensor(X, dtype=torch.float32)

def resolve_conditions(cache: dict, cell_name: str = None, time_name: str = None, dose_value: float = None):
    cell_vocab = cache["cell_vocab"]
    if cell_name and cell_name in cell_vocab["stoi"]:
        cell_id = cell_vocab["stoi"][cell_name]
        print(f"  Cell line : {cell_name} (id={cell_id})")
    else:
        if cell_name:
            print(f"  Warning: '{cell_name}' not in vocab. Available: {cell_vocab['itos']}")
        cell_id = int(torch.bincount(cache["cell"]).argmax())
        print(f"  Cell line : {cell_vocab['itos'][cell_id]} (most common, id={cell_id})")

    time_vocab = cache["time_vocab"]
    if time_name and time_name in time_vocab["stoi"]:
        time_id = time_vocab["stoi"][time_name]
        print(f"  Time point: {time_name}h (id={time_id})")
    else:
        if time_name:
            print(f"  Warning: '{time_name}' not in time vocab. Available: {time_vocab['itos']}")
        time_id = time_vocab["stoi"].get("24", int(torch.bincount(cache["time"]).argmax()))
        print(f"  Time point: {time_vocab['itos'][time_id]}h (default, id={time_id})")

    if dose_value is not None:
        dose_log_mean = float(cache["normalization"]["dose_log_mean"])
        dose_log_std  = float(cache["normalization"]["dose_log_std"])
        dose_norm = (math.log1p(dose_value) - dose_log_mean) / max(dose_log_std, 1e-8)
        print(f"  Dose      : {dose_value} µM (normalized={dose_norm:.3f})")
    else:
        dose_norm = float(cache["dose"].median())
        print(f"  Dose      : training median ({dose_norm:.3f} normalized)")

    return cell_id, time_id, dose_norm

def encode_profile(y, dose_norm, cell_id, time_id, weights):
    w = weights
    h = F.layer_norm(y, [978], w["profile_enc.y_in.0.weight"], w["profile_enc.y_in.0.bias"])
    h = F.linear(h, w["profile_enc.y_in.1.weight"], w["profile_enc.y_in.1.bias"])
    h = F.gelu(h)
    h = F.linear(h, w["profile_enc.y_in.4.weight"], w["profile_enc.y_in.4.bias"])

    for i in range(3):
        residual = h
        h_n = F.layer_norm(h, [384], w[f"profile_enc.y_blocks.{i}.norm.weight"], w[f"profile_enc.y_blocks.{i}.norm.bias"])
        h_n = F.linear(h_n, w[f"profile_enc.y_blocks.{i}.net.0.weight"], w[f"profile_enc.y_blocks.{i}.net.0.bias"])
        h_n = F.gelu(h_n)
        h_n = F.linear(h_n, w[f"profile_enc.y_blocks.{i}.net.3.weight"], w[f"profile_enc.y_blocks.{i}.net.3.bias"])
        h = residual + h_n

    d = torch.tensor([[dose_norm]], dtype=torch.float32)
    d = F.linear(d, w["profile_enc.dose_net.0.weight"], w["profile_enc.dose_net.0.bias"])
    d = F.gelu(d)
    d = F.linear(d, w["profile_enc.dose_net.2.weight"], w["profile_enc.dose_net.2.bias"])
    d = F.layer_norm(d, [384], w["profile_enc.dose_net.3.weight"], w["profile_enc.dose_net.3.bias"])

    c = w["profile_enc.cell_emb.weight"][cell_id].unsqueeze(0)
    t = w["profile_enc.time_emb.weight"][time_id].unsqueeze(0)

    fuse_in  = torch.cat([h, d, c, t, h, d], dim=-1)
    fuse_out = F.linear(fuse_in, w["profile_enc.fuse.0.weight"], w["profile_enc.fuse.0.bias"])
    fuse_out = F.gelu(fuse_out)
    fuse_out = F.linear(fuse_out, w["profile_enc.fuse.3.weight"], w["profile_enc.fuse.3.bias"])
    cond = F.layer_norm(fuse_out, [384], w["profile_enc.fuse.4.weight"], w["profile_enc.fuse.4.bias"])
    return cond

def sample_z(cond, weights, fixed=False):
    p = F.linear(cond, weights["prior.0.weight"], weights["prior.0.bias"])
    p = F.gelu(p)
    p = F.linear(p, weights["prior.2.weight"], weights["prior.2.bias"])
    mu, logvar = p.chunk(2, dim=-1)
    if fixed:
        return mu
    std = torch.exp(0.5 * logvar.clamp(-4, 4))
    return mu + std * torch.randn_like(std)

def get_h0(cond, z, weights):
    h_in = torch.cat([cond, z], dim=-1)
    h0   = F.linear(h_in, weights["h0.0.weight"], weights["h0.0.bias"])
    h0   = torch.tanh(h0)
    return h0.view(4, 1, 384)

def gru_step(x, h, layer, weights):
    wi = weights[f"decoder.gru.weight_ih_l{layer}"]
    wh = weights[f"decoder.gru.weight_hh_l{layer}"]
    bi = weights[f"decoder.gru.bias_ih_l{layer}"]
    bh = weights[f"decoder.gru.bias_hh_l{layer}"]
    gi = F.linear(x, wi, bi)
    gh = F.linear(h, wh, bh)
    r_i, z_i, n_i = gi.chunk(3, dim=-1)
    r_h, z_h, n_h = gh.chunk(3, dim=-1)
    r = torch.sigmoid(r_i + r_h)
    z = torch.sigmoid(z_i + z_h)
    n = torch.tanh(n_i + r * n_h)
    return (1 - z) * n + z * h

def generate_one(cond, z, weights, vocab, bos_id=1, eos_id=2, pad_id=0, max_len=128, temperature=1.0, top_k=80, top_p=0.97):
    h0     = get_h0(cond, z, weights)
    hidden = [h0[l] for l in range(4)]
    token_ids  = []
    curr_input = torch.tensor([bos_id])

    for _ in range(max_len):
        x = F.embedding(curr_input, weights["decoder.token_emb.weight"])
        for layer in range(4):
            hidden[layer] = gru_step(x, hidden[layer], layer, weights)
            x = hidden[layer]

        x_n    = F.layer_norm(x, [384], weights["decoder.norm.weight"], weights["decoder.norm.bias"])
        logits = F.linear(x_n, weights["decoder.out.weight"])
        logits = logits / max(temperature, 1e-8)

        if top_k > 0:
            thresh = torch.topk(logits, top_k).values[0, -1]
            logits[logits < thresh] = -1e9

        if top_p < 1.0:
            sl, si = torch.sort(logits, descending=True)
            cp     = torch.cumsum(F.softmax(sl, dim=-1), dim=-1)
            sl[cp - F.softmax(sl, dim=-1) > top_p] = -1e9
            logits = torch.zeros_like(logits).scatter_(1, si, sl)

        next_id    = torch.multinomial(F.softmax(logits, dim=-1), 1).item()
        curr_input = torch.tensor([next_id])

        if next_id == eos_id:
            break
        token_ids.append(next_id)

    selfies_str = "".join(vocab["itos"][t] for t in token_ids if t not in (bos_id, eos_id, pad_id))
    try:
        smiles = sf.decoder(selfies_str)
        if smiles and Chem.MolFromSmiles(smiles):
            return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
    except Exception:
        pass
    return None

@torch.no_grad()
def run_gex2sgen(csv_path: str, ckpt_path: str, cache_path: str, gene_info_path: str = None, num_samples: int = 100, cell_name: str = None, time_name: str = None, dose_value: float = None, temperature: float = 1.0, top_k: int = 80, top_p: float = 0.97) -> list[str]:
    print("=" * 55)
    print("Gex2SGen Inference")
    print("=" * 55)

    print("\n[1/4] Loading cache and checkpoint...")
    cache   = torch.load(cache_path,  map_location="cpu", weights_only=False)
    ckpt    = torch.load(ckpt_path,   map_location="cpu", weights_only=False)
    weights = ckpt["model"]
    vocab   = cache["vocab"]
    cfg     = cache["config"]

    print(f"  Method : {ckpt['method']} | epoch={ckpt['epoch']} | val_ce={ckpt['val']['ce']:.4f} | params={ckpt['n_params']:,}")
    print(f"  Vocab  : {cfg['vocab_size']} tokens | bos={cfg['bos_id']} eos={cfg['eos_id']}")
    
    print("\n[2/4] Loading and aligning expression data...")
    y_norm = load_and_align_expression(csv_path, cache, gene_info_path)

    print("\n[3/4] Resolving conditions...")
    cell_id, time_id, dose_norm = resolve_conditions(cache, cell_name, time_name, dose_value)

    bos_id = cfg["bos_id"]
    eos_id = cfg["eos_id"]
    pad_id = cfg["pad_id"]

    print(f"\n[4/4] Generating {num_samples} molecules × {len(y_norm)} sample(s)...")
    all_results = []

    for row_idx in range(len(y_norm)):
        y    = y_norm[row_idx:row_idx+1]
        cond = encode_profile(y, dose_norm, cell_id, time_id, weights)
        valid, seen = [], set()

        for s in range(num_samples):
            z     = sample_z(cond, weights)
            smiles = generate_one(
                cond, z, weights, vocab,
                bos_id=bos_id, eos_id=eos_id, pad_id=pad_id,
                max_len=cfg["max_len"],
                temperature=temperature,
                top_k=top_k,
                top_p=top_p,
            )
            if smiles and smiles not in seen:
                seen.add(smiles)
                valid.append(smiles)

        print(f"  Row {row_idx}: {len(valid)}/{num_samples} unique valid SMILES")
        all_results.extend(valid)

    total      = num_samples * len(y_norm)
    n_valid    = len(all_results)
    n_unique   = len(set(all_results))
    print(f"\n{'='*55}")
    print(f"Total generated : {total}")
    print(f"Valid           : {n_valid}  ({n_valid/total*100:.1f}%)")
    print(f"Unique valid    : {n_unique}")
    return list(set(all_results))

def main():
    parser = argparse.ArgumentParser(description="Orchestrator to run GEx2SGen on test datasets")
    parser.add_argument("--dataset", choices=DATASETS, help="Run on a specific dataset")
    parser.add_argument("--run-all", action="store_true", help="Run all datasets sequentially")
    
    args = parser.parse_args()
    if not (args.dataset or args.run_all):
        parser.print_help()
        sys.exit(1)
        
    datasets = DATASETS if args.run_all else [args.dataset]
    
    os.makedirs(OUTPUTS_DIR, exist_ok=True)
    
    ckpt_path = os.path.join(REIMPL_ROOT, "checkpoints/gex2sgen.pt")
    cache_path = os.path.join(REIMPL_ROOT, "cyclemol_cache.pt")
    gene_info_path = os.path.join(DATA_DIR, "LINCS/GSE70138/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz")
    
    for d in datasets:
        src_path = os.path.join(REST_METHODS_DIR, f"{d}.csv")
        if not os.path.exists(src_path):
            print(f"Error: Preprocessed file {src_path} does not exist. Run preprocessing first.")
            continue
            
        print(f"\n[{d}] Running GEx2SGen inference...")
        smiles_list = run_gex2sgen(
            csv_path=src_path,
            ckpt_path=ckpt_path,
            cache_path=cache_path,
            gene_info_path=gene_info_path,
            num_samples=100,
            temperature=1.0,
            top_k=80,
            top_p=0.97
        )
        
        df_out = pd.DataFrame(smiles_list, columns=['SMILES'])
        dest_path = os.path.join(OUTPUTS_DIR, f"{d}.csv")
        df_out.to_csv(dest_path, index=False)
        print(f"[{d}] Saved {len(df_out)} unique valid SMILES to {dest_path}")

if __name__ == "__main__":
    main()
