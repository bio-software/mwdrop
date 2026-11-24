#!/usr/bin/env python3
import argparse
import os
import sys
import tempfile

try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
except ImportError:
    print("[WARNING] PDBFixer not installed; trying without it!", file=sys.stderr)
    sys.exit(0)


def count_altloc_collapse(path):
    """Inspect the raw input for alternate-location conformers.

    PDBFixer/OpenMM keep a single conformer per atom on load, discarding the
    other alternate locations (column 17). Returns (positions, discarded):
    how many atom positions carried alternates, and how many alternate-location
    atoms will be dropped by collapsing each to one conformer. The discarded
    count is independent of which conformer is kept.
    """
    groups = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM  ", "HETATM")) and line[16] != " ":
                # Group by atom identity, ignoring resName so point-mutation
                # alternates at the same position still group together.
                key = (line[21], line[22:26], line[26], line[12:16])
                groups.setdefault(key, set()).add(line[16])
    positions = sum(1 for locs in groups.values() if len(locs) > 1)
    discarded = sum(len(locs) - 1 for locs in groups.values())
    return positions, discarded


def restore_b_factors(source_pdb, fixed_pdb):
    """Copy B-factors from the input PDB into PDBFixer's output, in place.

    OpenMM's PDBFile.writeFile writes every B-factor as 0.00. MobyWat counts a
    reference water only if its B-factor is below its limit (-bmax, default 30),
    so without this every crystal water would count. Atoms are matched by
    (chain, resSeq, iCode, atom name), which keepIds=True preserves; for
    alternate locations the first one listed in the input wins, as in PDBFixer.
    Atoms PDBFixer added keep 0.00. Returns (restored, total) atom counts.
    """
    def key(line):
        return (line[21], line[22:26], line[26], line[12:16].strip())

    b_factors = {}
    with open(source_pdb) as fh:
        for line in fh:
            if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 66:
                b_factors.setdefault(key(line), line[60:66])
    restored = total = 0
    with open(fixed_pdb) as fh:
        lines = fh.readlines()
    for i, line in enumerate(lines):
        if line.startswith(("ATOM  ", "HETATM")):
            total += 1
            b_factor = b_factors.get(key(line))
            if b_factor is not None:
                lines[i] = line[:60] + b_factor + line[66:]
                restored += 1
    with open(fixed_pdb, "w") as fh:
        fh.writelines(lines)
    return restored, total


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fix a PDB file (missing residues/atoms, nonstandard residues) in place."
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--target",
        action="store_true",
        help="Target PDB mode: remove existing waters and other heterogens.",
    )
    mode.add_argument(
        "--reference",
        action="store_true",
        help="Reference PDB mode: keep waters, remove other heterogens.",
    )
    parser.add_argument("input_pdb", help="Path to the PDB file to fix in place.")
    return parser.parse_args()


def main():
    args = parse_args()
    keep_water = args.reference

    input_pdb = os.path.abspath(args.input_pdb)
    if not os.path.isfile(input_pdb):
        print(f"Error: file not found: {input_pdb}", file=sys.stderr)
        sys.exit(1)

    input_dir = os.path.dirname(input_pdb)
    input_name = os.path.basename(input_pdb)

    fd, tmp_path = tempfile.mkstemp(
        prefix=input_name + ".",
        suffix=".tmp",
        dir=input_dir,
        text=True,
    )
    os.close(fd)

    alt_positions, alt_discarded = count_altloc_collapse(input_pdb)

    print(f"[INFO] Mode      : {'reference' if keep_water else 'target'}")
    print(f"[INFO] Input PDB : {input_pdb}")
    print(f"[INFO] TMP file  : {tmp_path}")

    try:
        fixer = PDBFixer(filename=input_pdb)

        # Find missing residues, but do not build missing residues at chain termini
        fixer.findMissingResidues()
        chains = list(fixer.topology.chains())
        for key in list(fixer.missingResidues.keys()):
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                del fixer.missingResidues[key]

        # Replace nonstandard residues if present
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()

        # Keep waters only in reference mode; always drop other heterogens
        fixer.removeHeterogens(keepWater=keep_water)

        # Add missing heavy atoms and residues
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        with open(tmp_path, "w") as handle:
            PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
            handle.flush()
            os.fsync(handle.fileno())

        if keep_water:
            b_restored, b_total = restore_b_factors(input_pdb, tmp_path)

        os.replace(tmp_path, input_pdb)

        if keep_water:
            print("[DONE] Waters kept; other heterogens removed.")
            print(f"[DONE] B-factors restored from the input for {b_restored} of {b_total} atoms "
                  f"(the rest were added by PDBFixer and keep 0.00).")
        else:
            print("[DONE] Waters and other heterogens removed.")
        print("[DONE] Missing heavy atoms and internal residues repaired where possible.")
        if alt_discarded:
            print(f"[DONE] Alternate conformers collapsed to a single conformer: "
                  f"kept 1 at each of {alt_positions} atom position(s), "
                  f"discarded {alt_discarded} alternate-location atom(s).")
        print(f"[DONE] Replaced original PDB with fixed version: {input_pdb}")

    except Exception:
        print(f"[ERROR] Failed while fixing PDB: {input_pdb}", file=sys.stderr)
        print(f"[ERROR] TMP file kept at: {tmp_path}", file=sys.stderr)
        raise


if __name__ == "__main__":
    main()
