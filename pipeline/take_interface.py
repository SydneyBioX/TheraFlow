import sys
import math
import pathlib

if len(sys.argv) != 6:
    print("Usage: python make_cluster_ligand.py input.pqr ligand_resname cutoff_A output_prefix ligand_charge")
    print("Example: python make_cluster_ligand.py complex.pqr LIG 4.0 drug_cluster 0")
    print("Example: python make_cluster_ligand.py complex.pqr STI 5.0 imatinib_cluster 1")
    sys.exit(1)

infile = sys.argv[1]
lig_resname = sys.argv[2].upper()
cutoff = float(sys.argv[3])
prefix = sys.argv[4]
ligand_charge = int(sys.argv[5])


def guess_element(atom_name, resname):
    """
    Guess element from atom name and residue name.
    """
    resname = resname.upper()
    atom_name = atom_name.upper().strip()
    
    # Standard protein/nucleic residues only contain H, C, N, O, S, P
    standard_residues = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "HIP", "HID", "HIE", "HSP", "HSD", "HSE", "ASH", "GLH", "LYN", "CYM", "TYM",
        "DA", "DC", "DG", "DT", "A", "C", "G", "U"
    }
    
    if resname in standard_residues:
        # For standard residues, the first letter of the atom name is the element
        # (with a few exceptions like hydrogen names starting with numbers, e.g., 1H)
        # Remove numbers first
        s = ''.join([c for c in atom_name if c.isalpha()])
        if not s:
            return "H" # Default for numbers-only
        first_char = s[0]
        if first_char in ["C", "H", "N", "O", "S", "P"]:
            return first_char
        return "C" # Default fallback
        
    # For non-standard residues (ligands, ions, water):
    # Remove numbers
    s = ''.join([c for c in atom_name if c.isalpha()])
    if not s:
        return "C"
        
    two_letter_elements = {
        "CL", "BR", "NA", "MG", "CA", "FE", "ZN", "CU",
        "MN", "CO", "NI", "SE", "SI", "AL", "LI", "HG",
        "CD", "PT", "AU", "AG", "AS", "SN", "PB", "I", "F"
    }
    
    if len(s) >= 2 and s[:2] in two_letter_elements:
        return s[:2].capitalize()
    
    if s[0] in ["C", "H", "N", "O", "S", "P", "F", "I", "B"]:
        return s[0]
        
    return s[0]


def try_float(x):
    try:
        return float(x)
    except Exception:
        return None


def parse_atom_line(line):
    """
    Parse ATOM/HETATM line.

    Supports:
    1. Fixed-column PDB/PQR style
    2. Whitespace-split PQR style

    Expected PQR-like format with chain:
    ATOM serial atom resname chain resseq x y z charge radius

    Example:
    HETATM 1234 C1 LIG A 501 12.345 6.789 1.234 -0.123 1.700
    """

    # Try fixed-column parsing first
    try:
        rec = line[0:6].strip()
        atom_name = line[12:16].strip()
        resname = line[17:20].strip().upper()
        chain = line[21].strip() or "_"
        resseq = int(line[22:26].strip())
        icode = line[26].strip() or ""

        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())

        element = line[76:78].strip() if len(line) >= 78 else ""
        if not element:
            element = guess_element(atom_name, resname)

        charge = None
        radius = None

        # Try to read PQR charge/radius from whitespace fields
        f = line.split()
        if len(f) >= 11:
            q = try_float(f[-2])
            r = try_float(f[-1])
            if q is not None and r is not None:
                charge = q
                radius = r

        return {
            "record": rec,
            "atom_name": atom_name,
            "resname": resname,
            "chain": chain,
            "resseq": resseq,
            "icode": icode,
            "x": x,
            "y": y,
            "z": z,
            "element": element,
            "charge": charge,
            "radius": radius,
            "line": line.rstrip("\n")
        }

    except Exception:
        pass

    # Fallback: whitespace-split PQR-like format
    f = line.split()

    if len(f) < 9:
        raise ValueError(f"Cannot parse atom line:\n{line}")

    rec = f[0]
    atom_name = f[2]
    resname = f[3].upper()

    # Format with chain:
    # ATOM serial atom resname chain resseq x y z charge radius
    try:
        chain = f[4] or "_"
        resseq = int(f[5])
        x, y, z = map(float, f[6:9])
        charge = try_float(f[9]) if len(f) > 9 else None
        radius = try_float(f[10]) if len(f) > 10 else None

    except Exception:
        # Format without chain:
        # ATOM serial atom resname resseq x y z charge radius
        chain = "_"
        resseq = int(f[4])
        x, y, z = map(float, f[5:8])
        charge = try_float(f[8]) if len(f) > 8 else None
        radius = try_float(f[9]) if len(f) > 9 else None

    element = guess_element(atom_name, resname)

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
        "element": element,
        "charge": charge,
        "radius": radius,
        "line": line.rstrip("\n")
    }


# ----------------------------
# Read atoms
# ----------------------------

atoms = []

with open(infile) as fh:
    for line in fh:
        if line.startswith(("ATOM", "HETATM")):
            atoms.append(parse_atom_line(line))

if not atoms:
    print(f"ERROR: No ATOM/HETATM records found in {infile}")
    sys.exit(1)


# ----------------------------
# Group atoms by residue
# ----------------------------

res_atoms = {}

for a in atoms:
    key = (a["chain"], a["resseq"], a["icode"], a["resname"])
    res_atoms.setdefault(key, []).append(a)


# ----------------------------
# Find ligand by residue name
# ----------------------------

lig_res_keys = [
    k for k in res_atoms
    if k[3].upper() == lig_resname
]

if not lig_res_keys:
    print(f"ERROR: No ligand residue with resname '{lig_resname}' found.")
    print()
    print("Available HETATM residue names:")

    het_resnames = sorted({
        a["resname"]
        for a in atoms
        if a["record"] == "HETATM"
    })

    if het_resnames:
        for r in het_resnames:
            n_atoms = sum(
                1 for a in atoms
                if a["record"] == "HETATM" and a["resname"] == r
            )
            print(f"  {r:>4s}  atoms={n_atoms}")
    else:
        print("  No HETATM residues found.")

    sys.exit(1)

lig_atoms = [
    a
    for k in lig_res_keys
    for a in res_atoms[k]
]

print("Found ligand residue(s):")
for k in lig_res_keys:
    print(f"  chain={k[0]} resseq={k[1]} icode={k[2] or '-'} resname={k[3]}")

print(f"Ligand atoms: {len(lig_atoms)}")


# ----------------------------
# Define protein residue names
# ----------------------------

protein_resnames = {
    # Standard amino acids
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",

    # Common protonation variants
    "HIP", "HID", "HIE",
    "HSP", "HSD", "HSE",
    "ASH", "GLH",
    "LYN",
    "CYM", "TYM"
}


# ----------------------------
# Find protein residues close to ligand
# ----------------------------

iface_keys = set()
cut2 = cutoff * cutoff

for k, ratoms in res_atoms.items():
    chain, resseq, icode, resname = k

    # Skip ligand itself
    if k in lig_res_keys:
        continue

    # Only select protein residues
    if resname.upper() not in protein_resnames:
        continue

    close = False

    for a in ratoms:
        ax, ay, az = a["x"], a["y"], a["z"]

        for l in lig_atoms:
            dx = ax - l["x"]
            dy = ay - l["y"]
            dz = az - l["z"]

            if dx * dx + dy * dy + dz * dz <= cut2:
                close = True
                break

        if close:
            iface_keys.add(k)
            break


selected_keys = set(lig_res_keys) | iface_keys

selected_atoms = [
    a for a in atoms
    if (a["chain"], a["resseq"], a["icode"], a["resname"]) in selected_keys
]


# ----------------------------
# Estimate cluster charge
# ----------------------------

charge_map = {
    "ASP": -1,
    "GLU": -1,
    "ARG": +1,
    "LYS": +1,

    # Histidine variants
    "HIP": +1,
    "HSP": +1,
    "HIS": 0,
    "HID": 0,
    "HIE": 0,
    "HSD": 0,
    "HSE": 0,

    # Neutral variants
    "ASH": 0,
    "GLH": 0,
    "LYN": 0,

    # Deprotonated variants
    "CYM": -1,
    "TYM": -1
}

protein_charge = sum(
    charge_map.get(k[3].upper(), 0)
    for k in iface_keys
)

approx_charge = protein_charge + ligand_charge


# ----------------------------
# Output files
# ----------------------------

out_pdb = prefix + ".pdb"
out_xyz = prefix + ".xyz"
out_charge = prefix + "_charge.txt"
out_list = prefix + "_residues.txt"


with open(out_pdb, "w") as fh:
    serial = 1

    for a in selected_atoms:
        fh.write(
            f"{a['record']:<6}"
            f"{serial:>5} "
            f"{a['atom_name']:<4} "
            f"{a['resname']:>3} "
            f"{a['chain']:1}"
            f"{a['resseq']:>4}"
            f"{a['icode']:1}   "
            f"{a['x']:>8.3f}"
            f"{a['y']:>8.3f}"
            f"{a['z']:>8.3f}"
            f"  1.00  0.00          "
            f"{a['element']:>2}\n"
        )
        serial += 1

    fh.write("END\n")


with open(out_xyz, "w") as fh:
    fh.write(f"{len(selected_atoms)}\n")
    fh.write(
        f"{prefix} "
        f"ligand={lig_resname} "
        f"cutoff={cutoff:.2f}A "
        f"protein_charge={protein_charge} "
        f"ligand_charge={ligand_charge} "
        f"approx_total_charge={approx_charge}\n"
    )

    for a in selected_atoms:
        fh.write(
            f"{a['element']:2s} "
            f"{a['x']: .6f} "
            f"{a['y']: .6f} "
            f"{a['z']: .6f}\n"
        )


with open(out_charge, "w") as fh:
    fh.write(str(approx_charge) + "\n")


with open(out_list, "w") as fh:
    fh.write("# Ligand residue(s)\n")
    for k in sorted(lig_res_keys):
        fh.write(f"LIGAND  {k[0]} {k[1]} {k[2] or '-'} {k[3]}\n")

    fh.write("\n# Protein residues within cutoff\n")
    for k in sorted(iface_keys):
        fh.write(f"PROTEIN {k[0]} {k[1]} {k[2] or '-'} {k[3]}\n")


# ----------------------------
# Summary
# ----------------------------

print()
print(f"Wrote {out_pdb}")
print(f"Wrote {out_xyz}")
print(f"Wrote {out_charge}")
print(f"Wrote {out_list}")
print()
print(f"Cutoff = {cutoff:.2f} Å")
print(f"Selected ligand residues = {len(lig_res_keys)}")
print(f"Selected protein interface residues = {len(iface_keys)}")
print(f"Selected atoms = {len(selected_atoms)}")
print(f"Protein charge estimate = {protein_charge}")
print(f"Ligand formal charge = {ligand_charge}")
print(f"Approx cluster charge = {approx_charge}")
