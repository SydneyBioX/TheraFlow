import os
import time
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
from chembl_webresource_client.new_client import new_client


def get_targets_and_neighbors_from_chembl(smiles, similarity_threshold=60):
    """Queries ChEMBL similarity search to identify targets and similar molecule IDs."""
    try:
        similarity = new_client.similarity
        res = similarity.filter(smiles=smiles, similarity=similarity_threshold).only(['molecule_chembl_id', 'similarity'])
        
        if not res:
            return {"targets": [], "neighbors": []}

        m_ids = [r['molecule_chembl_id'] for r in res]
        
        # Get targets where the similar molecules are active (IC50 < 10 uM)
        activities = new_client.activity
        target_hits = activities.filter(
            molecule_chembl_id__in=m_ids, 
            standard_type='IC50', 
            standard_value__lt=10000
        ).only(['molecule_chembl_id', 'target_pref_name'])
        
        # Extract unique targets and neighbors
        targets = list(set([t['target_pref_name'] for t in target_hits if t.get('target_pref_name')]))
        neighbors = list(set([t['molecule_chembl_id'] for t in target_hits]))
        
        return {"targets": targets, "neighbors": neighbors}
    except Exception as e:
        print(f"Error fetching targets/neighbors from ChEMBL for {smiles}: {e}")
        return {"targets": [], "neighbors": []}


def get_detailed_activity(smiles, similarity_threshold=80):
    """Retrieves detailed potency values (IC50/EC50) and mechanism actions from ChEMBL at similarity threshold."""
    results = []
    try:
        similarity = new_client.similarity
        activity = new_client.activity
        mechanism = new_client.mechanism
        
        search_results = similarity.filter(smiles=smiles, similarity=similarity_threshold).only(['molecule_chembl_id', 'similarity'])
        
        for res in search_results:
            chembl_id = res['molecule_chembl_id']
            sim_score = res['similarity']
            
            # Query mechanism of action (action type)
            moas = mechanism.filter(molecule_chembl_id=chembl_id).only(['action_type'])
            actions = list(set([m['action_type'] for m in moas if m.get('action_type')]))
            action_label = ", ".join(actions) if actions else "Unknown"
            
            # Query IC50/EC50 activity values
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
    """Helper to determine the consensus mechanism action (ignores 'Unknown' if possible)."""
    valid_actions = [a for a in series if a != "Unknown"]
    if not valid_actions:
        return "Unknown"
    return max(set(valid_actions), key=valid_actions.count)


def run_potency_chembl_pipeline(properties_csv, output_targets_csv, output_potency_csv, max_workers=5, checkpoint_frequency=10):
    """Main orchestrator for step 5: fetches ChEMBL targets, neighbors, and potency detailed activity."""
    df = pd.read_csv(properties_csv)
    distinct_smiles = df["SMILES"].unique().tolist()
    print(f"Total unique SMILES for ChEMBL search: {len(distinct_smiles)}")
    
    # -------------------------------------------------------------------------
    # Phase A: Get Targets and Neighbors (Similarity >= 60)
    # -------------------------------------------------------------------------
    print("Phase A: Retrieving ChEMBL targets and similar neighbors...")
    targets_results = []
    done_smiles = set()
    
    # Resume support
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
            futures = list(tqdm(executor.map(
                lambda s: get_targets_and_neighbors_from_chembl(s, similarity_threshold=60),
                remaining_smiles
            ), total=len(remaining_smiles)))
            
            for i, data in enumerate(futures):
                smiles = remaining_smiles[i]
                targets_results.append({
                    "SMILES": smiles,
                    "Predicted_Targets": data["targets"],
                    "ChEMBL_Neighbors": data["neighbors"]
                })
                
                # Checkpoint save
                if (i + 1) % checkpoint_frequency == 0:
                    pd.DataFrame(targets_results).to_csv(output_targets_csv, index=False)
                    time.sleep(0.5)
                    
        pd.DataFrame(targets_results).to_csv(output_targets_csv, index=False)
        
    print(f"Phase A Complete. Targets results saved to {output_targets_csv}")
    
    # -------------------------------------------------------------------------
    # Phase B: Get Detailed Potency (Similarity >= 80)
    # -------------------------------------------------------------------------
    print("Phase B: Retrieving detailed potency values...")
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
            futures = list(tqdm(executor.map(
                lambda s: get_detailed_activity(s, similarity_threshold=80),
                remaining_potency_smiles
            ), total=len(remaining_potency_smiles)))
            
            for i, res_list in enumerate(futures):
                activity_results.extend(res_list)
                
                # Checkpoint save
                if (i + 1) % checkpoint_frequency == 0:
                    pd.DataFrame(activity_results).to_csv(temp_activity_file, index=False)
                    time.sleep(0.5)
                    
        pd.DataFrame(activity_results).to_csv(temp_activity_file, index=False)
        
    print("Consolidating potency data and generating consensus scores...")
    if not activity_results:
        # Create empty potency file with appropriate headers if no activities found
        print("No ChEMBL activity matches found.")
        pd.DataFrame(columns=["query_smiles", "IC50", "EC50", "similarity", "target_compound", "action"]).to_csv(output_potency_csv, index=False)
        return
        
    df_act = pd.DataFrame(activity_results)
    df_act['value'] = pd.to_numeric(df_act['value'], errors='coerce')
    
    # Filter for standard nM units and drop empty values
    df_act = df_act[df_act['unit'] == "nM"].dropna(subset=['value'])
    
    if len(df_act) == 0:
        print("No activities found in 'nM' units. Consensus potency scores will be empty.")
        pd.DataFrame(columns=["query_smiles", "IC50", "EC50", "similarity", "target_compound", "action"]).to_csv(output_potency_csv, index=False)
        return
        
    # Pivot to get consensus IC50 and EC50 columns
    final_consensus = df_act.pivot_table(
        index='query_smiles', 
        columns='type', 
        values='value', 
        aggfunc='mean'
    ).reset_index()
    
    # Extract metadata metrics
    best_metadata = df_act.groupby('query_smiles').agg({
        'similarity': 'max',
        'target_compound': 'first',
        'action': get_best_action,
    }).reset_index()
    
    final_consensus = final_consensus.merge(best_metadata, on='query_smiles', how='left')
    final_consensus.columns.name = None
    
    # Save the final consensus potency file
    final_consensus.to_csv(output_potency_csv, index=False)
    print(f"Phase B Complete. Consensus potency scores saved to {output_potency_csv}")
    
    # Remove raw temporary file
    if os.path.exists(temp_activity_file):
        try:
            os.remove(temp_activity_file)
        except:
            pass
            
    return final_consensus
