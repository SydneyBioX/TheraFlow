#!/usr/bin/env python3
import sys
import math
import pathlib

if len(sys.argv) < 5:
    print("Usage: python take_interface.py input.pqr ligand_resname cutoff_A output_prefix [ligand_charge]")
    sys.exit(1)

infile = sys.argv[1]
lig_resname = sys.argv[2].upper()
cutoff = float(sys.argv[3])
prefix = sys.argv[4]
ligand_charge = int(sys.argv[5]) if len(sys.argv) > 5 else 0


def guess_element(atom_name, resname):
    """Guess element from atom name and residue name."""
    resname = resname.upper()
    atom_name = atom_name.upper().strip()
    
    standard_residues = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "HIP", "HID", "HIE", "HSP", "HSD", "HSE", "ASH", "GLH", "LYN", "CYM", "TYM"
    }
    
    s = ''.join([c for c in atom_name if c.isalpha()])
    if not s:
        return "C"
        
    if resname in standard_residues:
        return s[0] if s[0] in ["C", "H", "N", "O", "S", "P"] else "C"
        
    two_letter_elements = {
        "CL", "BR", "NA", "MG", "CA", "FE", "ZN", "CU",
        "MN", "CO", "NI", "SE", "SI", "AL", "LI", "HG",
        "CD", "PT", "AU", "AG", "AS", "SN", "PB", "I", "F"
    }
    
    if len(s) >= 2 and s[:2] in two_letter_elements:
        return s[:2].capitalize()
    if s[0] in ["C", "H", "N", "O", "S", "P", "F", "I", "B"]:
        return s[0]
    return "C"


def try_float(val):
    try:
        return float(val)
    except Exception:
        return None


def try_int(val, default=1):
    try:
        return int(val)
    except Exception:
        return default


def parse_atom_line(line):
    """Parse ATOM/HETATM line flexibly from PDB, PQR, or CIF-derived outputs."""
    line_clean = line.rstrip("\n")
    if not line_clean.startswith(("ATOM", "HETATM")):
        return None
        
    f = line_clean.split()
    if len(f) < 7:
        return None
        
    rec = f[0]
    atom_name = f[2] if len(f) > 2 else f[1]
    resname = f[3].upper() if len(f) > 3 else "UNK"
    
    # Check if ligand
    is_lig = lig_resname in resname or "LIG" in resname
    
    # Try fixed column first if length >= 54
    if len(line_clean) >= 54:
        try:
            aname = line_clean[12:16].strip()
            rname = line_clean[17:20].strip().upper()
            ch = line_clean[21].strip() or ("B" if is_lig else "A")
            rseq_str = line_clean[22:26].strip()
            rseq = try_int(rseq_str, 1)
            
            x = float(line_clean[30:38].strip())
            y = float(line_clean[38:46].strip())
            z = float(line_clean[46:54].strip())
            
            elem = line_clean[76:78].strip() if len(line_clean) >= 78 else ""
            if not elem:
                elem = guess_element(aname, rname)
                
            return {
                "record": rec,
                "atom_name": aname,
                "resname": rname,
                "chain": ch,
                "resseq": rseq,
                "icode": "",
                "x": x,
                "y": y,
                "z": z,
                "element": elem,
                "charge": 0.0,
                "radius": 1.5,
                "line": line_clean
            }
        except Exception:
            pass
            
    # Flexible token search for coordinates (find 3 consecutive floats)
    floats = []
    float_indices = []
    for i, token in enumerate(f):
        val = try_float(token)
        if val is not None and "." in token:
            floats.append(val)
            float_indices.append(i)
            
    if len(floats) >= 3:
        # Find consecutive 3 floats
        x, y, z = floats[0], floats[1], floats[2]
        c_idx = float_indices[0]
        
        # Everything before c_idx is metadata
        meta = f[:c_idx]
        atom_name = meta[2] if len(meta) > 2 else meta[1]
        resname = meta[3].upper() if len(meta) > 3 else "UNK"
        chain = "B" if is_lig else "A"
        resseq = 1
        
        if len(meta) > 4:
            if len(meta[4]) == 1 and meta[4].isalpha():
                chain = meta[4]
            else:
                resseq = try_int(meta[4], 1)
        if len(meta) > 5:
            resseq = try_int(meta[5], resseq)
            
        elem = guess_element(atom_name, resname)
        return {
            "record": rec,
            "atom_name": atom_name,
            "resname": resname,
            "chain": chain,
            "resseq": resseq,
            "icode": "",
            "x": x,
            "y": y,
            "z": z,
            "element": elem,
            "charge": 0.0,
            "radius": 1.5,
            "line": line_clean
        }
    return None


def dist(a, b):
    return math.sqrt((a["x"] - b["x"])**2 + (a["y"] - b["y"])**2 + (a["z"] - b["z"])**2)


def main():
    atoms = []
    with open(infile, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                parsed = parse_atom_line(line)
                if parsed:
                    atoms.append(parsed)
                    
    if not atoms:
        print(f"Error: No atoms parsed from {infile}")
        sys.exit(1)
        
    ligand_atoms = [a for a in atoms if (lig_resname in a["resname"] or "LIG" in a["resname"])]
    if not ligand_atoms:
        print(f"Warning: No ligand atoms found with resname {lig_resname}, trying HETATM records...")
        ligand_atoms = [a for a in atoms if a["record"] == "HETATM"]
        
    if not ligand_atoms:
        print(f"Error: No ligand atoms found in {infile}")
        sys.exit(1)
        
    # Find all protein residues within cutoff of any ligand atom
    selected_reskeys = set()
    for lig in ligand_atoms:
        selected_reskeys.add((lig["chain"], lig["resseq"], lig["icode"], lig["resname"]))
        
    protein_atoms = [a for a in atoms if a not in ligand_atoms]
    for p in protein_atoms:
        for lig in ligand_atoms:
            if dist(p, lig) <= cutoff:
                selected_reskeys.add((p["chain"], p["resseq"], p["icode"], p["resname"]))
                break
                
    cluster_atoms = [a for a in atoms if (a["chain"], a["resseq"], a["icode"], a["resname"]) in selected_reskeys]
    
    # Write cluster.xyz
    xyz_out = f"{prefix}.xyz"
    with open(xyz_out, "w") as f:
        f.write(f"{len(cluster_atoms)}\n")
        f.write(f"Active site cluster ({cutoff} A around {lig_resname})\n")
        for a in cluster_atoms:
            f.write(f"{a['element']:<2} {a['x']:12.6f} {a['y']:12.6f} {a['z']:12.6f}\n")
            
    # Write cluster.pdb
    pdb_out = f"{prefix}.pdb"
    with open(pdb_out, "w") as f:
        for i, a in enumerate(cluster_atoms, 1):
            f.write(f"{a['record']:<6}{i:>5} {a['atom_name']:<4} {a['resname']:<3} {a['chain']}{a['resseq']:>4}    {a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}  1.00 50.00          {a['element']:>2}\n")
        f.write("END\n")
        
    # Write cluster_charge.txt
    chrg_out = f"{prefix}_charge.txt"
    with open(chrg_out, "w") as f:
        f.write(f"{ligand_charge}\n")
        
    print(f"Extracted {len(cluster_atoms)} cluster atoms ({len(ligand_atoms)} ligand atoms) around {lig_resname}")


if __name__ == "__main__":
    main()
