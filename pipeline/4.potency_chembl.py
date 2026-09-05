#!/usr/bin/env python3
"""
Step 4: ChEMBL Target and Potency Retrieval
Performs similarity searches against ChEMBL to identify predicted gene targets (>=60% similarity)
and detailed activity/IC50/EC50 values and Mechanism of Action (>=80% similarity).

Corresponds to Chunks 18-20 of metrics_calculation/7.1. results_5.ipynb.
"""

import os
import time
import argparse
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor


def get_chembl_client():
    """Lazily initializes and returns the ChEMBL client."""
    try:
        from chembl_webresource_client.new_client import new_client
        return new_client
    except ImportError:
        raise ImportError(
            "chembl_webresource_client is not installed. "
            "Please install it using: pip install chembl_webresource_client"
        )


def get_targets_and_neighbors_from_chembl(smiles, similarity_threshold=60):
    """Retrieves target genes and neighboring compound IDs from ChEMBL at similarity threshold (default 60%)."""
    new_client = get_chembl_client()
    try:
        similarity = new_client.similarity
        activity = new_client.activity
        
        res = list(similarity.filter(smiles=smiles, similarity=similarity_threshold).only(
            ['molecule_chembl_id', 'similarity']
        )[:20])
        
        if not res:
            return {'targets': "None", 'neighbors': "None"}
            
        m_ids = [r['molecule_chembl_id'] for r in res if r.get('molecule_chembl_id')]
        sim_map = {r['molecule_chembl_id']: r.get('similarity') for r in res if r.get('molecule_chembl_id')}
        
        target_hits = list(activity.filter(
            molecule_chembl_id__in=m_ids[:15],
            standard_type='IC50',
            standard_value__lt=10000
        ).only(['molecule_chembl_id', 'target_pref_name'])[:30])
        
        targets = set()
        neighbors = set()
        for t in target_hits:
            t_name = t.get('target_pref_name')
            m_id = t.get('molecule_chembl_id')
            if t_name and t_name not in ['Unchecked', 'Cell', 'Organism', 'Homo sapiens', 'Mus musculus']:
                targets.add(t_name)
            if m_id:
                sim_val = sim_map.get(m_id, "")
                neighbors.add(f"{m_id} ({sim_val}%)" if sim_val else m_id)
                
        if not targets:
            # Fallback to any activity
            target_hits_fallback = list(activity.filter(
                molecule_chembl_id__in=m_ids[:10]
            ).only(['molecule_chembl_id', 'target_pref_name'])[:20])
            for t in target_hits_fallback:
                t_name = t.get('target_pref_name')
                m_id = t.get('molecule_chembl_id')
                if t_name and t_name not in ['Unchecked', 'Cell', 'Organism', 'Homo sapiens', 'Mus musculus']:
                    targets.add(t_name)
                if m_id:
                    sim_val = sim_map.get(m_id, "")
                    neighbors.add(f"{m_id} ({sim_val}%)" if sim_val else m_id)

        return {
            'targets': ", ".join(sorted(targets)) if targets else "None",
            'neighbors': ", ".join(sorted(neighbors)) if neighbors else "None"
        }
    except Exception as e:
        print(f"Error fetching targets from ChEMBL for {smiles}: {e}")
        return {'targets': "None", 'neighbors': "None"}


def get_detailed_activity(smiles, similarity_threshold=80):
    """Retrieves detailed potency values (IC50/EC50) and mechanism actions from ChEMBL at similarity threshold (default 80%)."""
    new_client = get_chembl_client()
    results = []
    try:
        similarity = new_client.similarity
        activity = new_client.activity
        mechanism = new_client.mechanism
        
        search_results = similarity.filter(smiles=smiles, similarity=similarity_threshold).only(['molecule_chembl_id', 'similarity'])
        
        for res in search_results:
            chembl_id = res['molecule_chembl_id']
            sim_score = res['similarity']
            
            moas = mechanism.filter(molecule_chembl_id=chembl_id).only(['action_type'])
            actions = list(set([m['action_type'] for m in moas if m.get('action_type')]))
            action_label = ", ".join(actions) if actions else "Unknown"
            
            acts = activity.filter(molecule_chembl_id=chembl_id, standard_type__in=['IC50', 'EC50']).only(
                ['target_chembl_id', 'standard_type', 'standard_value', 'standard_units', 'description', 'assay_type']
            )
            
            for a in acts:
                results.append({
                    'query_smiles': smiles,
                    'target_compound': chembl_id,
                    'similarity': sim_score,
                    'type': a.get('standard_type'),
                    'value': a.get('standard_value'),
                    'unit': a.get('standard_units'),
                    'action': action_label
                })
        return results
    except Exception as e:
        print(f"Error fetching detailed activity from ChEMBL for {smiles}: {e}")
        return []


def get_best_action(series):
    """Helper to determine consensus mechanism of action (ignores 'Unknown' if valid action exists)."""
    valid_actions = [a for a in series if a != "Unknown"]
    if not valid_actions:
        return "Unknown"
    return max(set(valid_actions), key=valid_actions.count)


def run_potency_chembl_pipeline(properties_csv, output_targets_csv, output_potency_csv, max_workers=5, checkpoint_frequency=10):
    """Full execution: queries ChEMBL targets, neighbors, and detailed potency activities."""
    df = pd.read_csv(properties_csv)
    distinct_smiles = df["SMILES"].dropna().unique().tolist()
    print(f"Total unique SMILES for ChEMBL query: {len(distinct_smiles)}")
    
    # -------------------------------------------------------------------------
    # Phase A: Get Targets and Neighbors (Similarity >= 60%)
    # -------------------------------------------------------------------------
    print("Phase A: Retrieving ChEMBL targets and similar neighbors (>=60% similarity)...")
    targets_results = []
    done_smiles = set()
    
    if os.path.exists(output_targets_csv):
        try:
            done_df = pd.read_csv(output_targets_csv)
            done_smiles = set(done_df["SMILES"].tolist())
            targets_results = done_df.to_dict("records")
            print(f"Resuming targets search: {len(done_smiles)} already completed, {len(distinct_smiles) - len(done_smiles)} remaining.")
        except Exception as e:
            print(f"Warning: Could not parse target checkpoint CSV: {e}. Starting fresh.")
            
    remaining_smiles = [s for s in distinct_smiles if s not in done_smiles]
    
    if remaining_smiles:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            target_iter = executor.map(
                lambda s: get_targets_and_neighbors_from_chembl(s, similarity_threshold=60),
                remaining_smiles
            )
            for i, (smiles, data) in enumerate(zip(remaining_smiles, tqdm(target_iter, total=len(remaining_smiles), desc="Querying Targets"))):
                targets_results.append({
                    "SMILES": smiles,
                    "Predicted_Targets": data["targets"],
                    "ChEMBL_Neighbors": data["neighbors"]
                })
                if (i + 1) % checkpoint_frequency == 0:
                    os.makedirs(os.path.dirname(output_targets_csv) if os.path.dirname(output_targets_csv) else ".", exist_ok=True)
                    pd.DataFrame(targets_results).to_csv(output_targets_csv, index=False)
                    
        os.makedirs(os.path.dirname(output_targets_csv) if os.path.dirname(output_targets_csv) else ".", exist_ok=True)
        pd.DataFrame(targets_results).to_csv(output_targets_csv, index=False)
        
    print(f"Phase A Complete. Targets saved to {output_targets_csv}")
    
    # -------------------------------------------------------------------------
    # Phase B: Get Detailed Potency Activity (Similarity >= 80%)
    # -------------------------------------------------------------------------
    print("\nPhase B: Retrieving detailed potency values (>=80% similarity)...")
    activity_results = []
    done_potency_smiles = set()
    
    temp_activity_file = str(output_potency_csv) + ".raw_activities.csv"
    if os.path.exists(temp_activity_file):
        try:
            done_df = pd.read_csv(temp_activity_file)
            done_potency_smiles = set(done_df["query_smiles"].unique().tolist())
            activity_results = done_df.to_dict("records")
            print(f"Resuming potency search: {len(done_potency_smiles)} already completed, {len(distinct_smiles) - len(done_potency_smiles)} remaining.")
        except Exception as e:
            print(f"Warning: Could not parse potency raw activities CSV: {e}. Starting fresh.")
            
    remaining_potency_smiles = [s for s in distinct_smiles if s not in done_potency_smiles]
    
    if remaining_potency_smiles:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            pot_iter = executor.map(
                lambda s: get_detailed_activity(s, similarity_threshold=80),
                remaining_potency_smiles
            )
            for i, res_list in enumerate(tqdm(pot_iter, total=len(remaining_potency_smiles), desc="Querying Potency")):
                activity_results.extend(res_list)
                if (i + 1) % checkpoint_frequency == 0:
                    os.makedirs(os.path.dirname(temp_activity_file) if os.path.dirname(temp_activity_file) else ".", exist_ok=True)
                    pd.DataFrame(activity_results).to_csv(temp_activity_file, index=False)
                    
        os.makedirs(os.path.dirname(temp_activity_file) if os.path.dirname(temp_activity_file) else ".", exist_ok=True)
        pd.DataFrame(activity_results).to_csv(temp_activity_file, index=False)
        
    print("Consolidating potency data and generating consensus scores...")
    os.makedirs(os.path.dirname(output_potency_csv) if os.path.dirname(output_potency_csv) else ".", exist_ok=True)
    
    if not activity_results:
        print("No ChEMBL activity matches found.")
        pd.DataFrame(columns=["query_smiles", "IC50", "EC50", "similarity", "target_compound", "action"]).to_csv(output_potency_csv, index=False)
        return pd.DataFrame(columns=["query_smiles", "IC50", "EC50", "similarity", "target_compound", "action"])
        
    df_act = pd.DataFrame(activity_results)
    df_act['value'] = pd.to_numeric(df_act['value'], errors='coerce')
    df_act = df_act[df_act['unit'] == "nM"].dropna(subset=['value'])
    
    if len(df_act) == 0:
        print("No activities found in 'nM' units. Consensus potency scores will be empty.")
        pd.DataFrame(columns=["query_smiles", "IC50", "EC50", "similarity", "target_compound", "action"]).to_csv(output_potency_csv, index=False)
        return pd.DataFrame(columns=["query_smiles", "IC50", "EC50", "similarity", "target_compound", "action"])
        
    final_consensus = df_act.pivot_table(
        index='query_smiles', 
        columns='type', 
        values='value', 
        aggfunc='mean'
    ).reset_index()
    
    best_metadata = df_act.groupby('query_smiles').agg({
        'similarity': 'max',
        'target_compound': 'first',
        'action': get_best_action,
    }).reset_index()
    
    final_consensus = final_consensus.merge(best_metadata, on='query_smiles', how='left')
    final_consensus.columns.name = None
    
    final_consensus.to_csv(output_potency_csv, index=False)
    print(f"Phase B Complete. Consensus potency saved to {output_potency_csv}")
    
    if os.path.exists(temp_activity_file):
        try:
            os.remove(temp_activity_file)
        except:
            pass
            
    return final_consensus


def main():
    parser = argparse.ArgumentParser(description="Query ChEMBL for similar compounds, targets, and potency activity.")
    parser.add_argument(
        "--input", 
        type=str, 
        default="TheraFlow/outputs/molecular_properties_filtered.csv",
        help="Input CSV properties file (filtered or full)."
    )
    parser.add_argument(
        "--output_targets", 
        type=str, 
        default="TheraFlow/outputs/chembl_targets.csv",
        help="Output CSV for predicted targets and neighbors."
    )
    parser.add_argument(
        "--output_potency", 
        type=str, 
        default="TheraFlow/outputs/potency_scores.csv",
        help="Output CSV for consensus potency and mechanism of action."
    )
    parser.add_argument(
        "--max_workers", 
        type=int, 
        default=5,
        help="Number of concurrent query workers."
    )
    args = parser.parse_args()
    
    input_path = Path(args.input)
    if not input_path.exists():
        fallback = Path("TheraFlow/outputs/molecular_properties_metrics.csv")
        if fallback.exists():
            input_path = fallback
        else:
            raise FileNotFoundError(f"Input file not found: {args.input}")
            
    run_potency_chembl_pipeline(
        properties_csv=str(input_path),
        output_targets_csv=args.output_targets,
        output_potency_csv=args.output_potency,
        max_workers=args.max_workers
    )


if __name__ == "__main__":
    main()
