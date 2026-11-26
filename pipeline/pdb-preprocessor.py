#!/usr/bin/env python3
import argparse
import os
import sys
import tempfile

try:
    from pdbfixer import PDBFixer
    from pdbfixer.pdbfixer import dnaResidues, proteinResidues, rnaResidues
    from openmm.app import Modeller, PDBFile, Topology
    from openmm.unit import nanometer
except ImportError:
    print("[WARNING] PDBFixer not installed; trying without it!", file=sys.stderr)
    sys.exit(0)

# Standard nucleotide residue names (after PDBFixer's replaceNonstandardResidues).
NUCLEOTIDE_RESNAMES = {"A", "C", "G", "U", "I", "DA", "DC", "DG", "DT", "DU", "DI"}

# Phosphate atoms a 5'-terminal nucleotide may carry (OpenMM maps legacy names
# such as O1P/O2P to OP1/OP2 on load).
FIVE_PRIME_PHOSPHATE_ATOMS = {"P", "OP1", "OP2", "OP3"}

# Monatomic ions kept in --target mode: the ones the force fields know by residue name
# (GROMACS amber99sb-ildn aminoacids.rtp / ions.itp, Tinker amber99.prm biotypes).
FORCE_FIELD_IONS = {"MG", "NA", "K", "CA", "ZN", "CL"}

# Candidate chain IDs for the kept ions, tried in order; the first unused one wins.
ION_CHAIN_ID_CANDIDATES = "ZYXWVUTSRQPONMLKJIHGFEDCBA0123456789"


def strip_five_prime_phosphates(fixer):
    """Remove the phosphate group from the first nucleotide of every chain.

    The AMBER force fields have no 5'-phosphate terminal residue: every 5'
    building block (GROMACS rtp DA5/RU5/..., single-residue DAN/RUN/...) starts
    at a 5'-OH, so pdb2gmx aborts with "Atom P in residue ... was not found in
    rtp entry ..5" when a chain starts with its phosphate. Tinker's pdbxyz drops
    those atoms silently; doing it here makes both engines see the same system.
    OpenMM starts a new chain at every chain-ID change and at every TER record,
    the same rule as pdb2gmx (-chainsep id_or_ter), so each strand is handled.
    PDBFixer itself never adds P/OP1/OP2 to a 5' terminus, so nothing re-adds
    them. Returns a list of "chain:resname resid" labels that were stripped.
    """
    to_delete = []
    stripped = []
    for chain in fixer.topology.chains():
        residues = list(chain.residues())
        if not residues or residues[0].name not in NUCLEOTIDE_RESNAMES:
            continue
        first = residues[0]
        atoms = [a for a in first.atoms() if a.name in FIVE_PRIME_PHOSPHATE_ATOMS]
        if atoms:
            to_delete.extend(atoms)
            stripped.append(f"{chain.id}:{first.name}{first.id}{first.insertionCode.strip()}")
    if to_delete:
        modeller = Modeller(fixer.topology, fixer.positions)
        modeller.delete(to_delete)
        fixer.topology = modeller.topology
        fixer.positions = modeller.positions
    return stripped


DNA_TO_RNA = {"DA": "A", "DC": "C", "DG": "G", "DT": "U", "DU": "U"}
RNA_TO_DNA = {"A": "DA", "C": "DC", "G": "DG", "U": "DT"}
# Atoms the converted residue's template lacks; PDBFixer's addMissingAtoms then
# builds the ones it needs (O2' for DNA->RNA, the thymine C7 for U->DT).
RNA_ONLY_ATOMS = {"O2'", "HO2'"}
DNA_TO_RNA_DROP = {"C7", "H71", "H72", "H73"}


def unify_nucleic_chain_types(fixer):
    """Make every nucleic-acid chain all-DNA or all-RNA.

    pdb2gmx refuses a chain whose residues are not one type ("The residues in
    the chain ... do not have a consistent type"). The usual cause is a modified
    nucleotide whose MODRES parent is of the other type, e.g. 1J9H: RNA with
    5-bromo-2'-deoxycytidine (CBR -> DC), which replaceNonstandardResidues turns
    into a DC inside an RNA chain. Tinker's pdbxyz accepts the mix, GROMACS does
    not. The minority residues are converted to the chain's majority type (a tie
    is left as is: pdb2gmx then reports it). Call after
    replaceNonstandardResidues and before findMissingAtoms. Returns a list of
    "chain:OLD->NEW resid" labels.
    """
    converted = []
    to_delete = []
    for chain in fixer.topology.chains():
        residues = list(chain.residues())
        dna = [r for r in residues if r.name in DNA_TO_RNA]
        rna = [r for r in residues if r.name in RNA_TO_DNA]
        if not dna or not rna or len(dna) == len(rna):
            continue
        if len(rna) > len(dna):
            minority, mapping, drop = dna, DNA_TO_RNA, DNA_TO_RNA_DROP
        else:
            minority, mapping, drop = rna, RNA_TO_DNA, RNA_ONLY_ATOMS
        for res in minority:
            new_name = mapping[res.name]
            converted.append(f"{chain.id}:{res.name}->{new_name}{res.id}{res.insertionCode.strip()}")
            res.name = new_name
            to_delete.extend(a for a in res.atoms() if a.name in drop)
    if to_delete:
        modeller = Modeller(fixer.topology, fixer.positions)
        modeller.delete(to_delete)
        fixer.topology = modeller.topology
        fixer.positions = modeller.positions
    return converted


def keep_ions_enabled():
    """Crystal ions are kept by default; MWDROP_KEEP_IONS=0 (or false/no/off) opts out."""
    return os.environ.get("MWDROP_KEEP_IONS", "1").strip().lower() not in ("0", "false", "no", "off")


def remove_heterogens_except_ions(fixer):
    """Remove waters and heterogens like PDBFixer.removeHeterogens(keepWater=False),
    but keep the force-field ions (FORCE_FIELD_IONS).

    Ions shape the hydration of nucleic acids in particular (a Mg2+ holds a shell of
    six waters, phosphates are screened), so dropping them leaves those crystal
    waters with nothing to reproduce them. Returns the kept ions as
    "chain:resname+resSeq" labels.
    """
    keep = set(proteinResidues) | set(dnaResidues) | set(rnaResidues) | {"N", "UNK"} | FORCE_FIELD_IONS
    to_delete = [r for r in fixer.topology.residues() if r.name not in keep]
    kept_ions = [f"{r.chain.id}:{r.name}{r.id}" for r in fixer.topology.residues()
                 if r.name in FORCE_FIELD_IONS]

    modeller = Modeller(fixer.topology, fixer.positions)
    modeller.delete(to_delete)
    fixer.topology = modeller.topology
    fixer.positions = modeller.positions
    return kept_ions


def move_ions_to_own_chain(fixer):
    """Move every ion into one new chain at the end of the structure.

    A crystal ion usually carries the chain ID of the molecule it binds and may be
    written inside that chain, so pdb2gmx would see a chain mixing e.g. 'DNA' and
    'Ion' residues and abort ("do not have a consistent type"). The topology is
    rebuilt in the same atom order, except that the ions move to a trailing chain
    with the first unused ID from ION_CHAIN_ID_CANDIDATES. Coordinates are copied
    unchanged. Returns that chain ID, or None when there are no ions.
    """
    ions = [r for r in fixer.topology.residues() if r.name in FORCE_FIELD_IONS]
    if not ions:
        return None
    used_ids = {chain.id for chain in fixer.topology.chains()}
    ion_chain_id = next(c for c in ION_CHAIN_ID_CANDIDATES if c not in used_ids)

    old_positions = fixer.positions.value_in_unit(nanometer)
    new_topology = Topology()
    new_topology.setPeriodicBoxVectors(fixer.topology.getPeriodicBoxVectors())
    new_positions = []
    new_atom_for = {}  # old atom -> new atom, to carry the bonds over

    def copy_residue(residue, chain):
        new_residue = new_topology.addResidue(residue.name, chain, residue.id, residue.insertionCode)
        for atom in residue.atoms():
            new_atom_for[atom] = new_topology.addAtom(atom.name, atom.element, new_residue, atom.id)
            new_positions.append(old_positions[atom.index])

    for chain in fixer.topology.chains():
        residues = [r for r in chain.residues() if r.name not in FORCE_FIELD_IONS]
        if residues:
            new_chain = new_topology.addChain(chain.id)
            for residue in residues:
                copy_residue(residue, new_chain)

    ion_chain = new_topology.addChain(ion_chain_id)
    for residue in ions:
        copy_residue(residue, ion_chain)

    for atom1, atom2 in fixer.topology.bonds():
        new_topology.addBond(new_atom_for[atom1], new_atom_for[atom2])

    fixer.topology = new_topology
    fixer.positions = new_positions * nanometer
    return ion_chain_id


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
        help="Target PDB mode: remove existing waters and other heterogens; force-field "
             "ions (Mg, Na, K, Ca, Zn, Cl) are kept unless MWDROP_KEEP_IONS=0.",
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

        # Target only (pdb2gmx needs single-type chains); before findMissingAtoms,
        # which then adds the atoms the converted residues' templates need.
        unified = [] if keep_water else unify_nucleic_chain_types(fixer)

        # Reference: keep waters, drop all other heterogens. Target: drop waters and
        # heterogens but keep the force-field ions, unless MWDROP_KEEP_IONS=0.
        keep_ions = not keep_water and keep_ions_enabled()
        kept_ions = []
        if keep_ions:
            kept_ions = remove_heterogens_except_ions(fixer)
        else:
            fixer.removeHeterogens(keepWater=keep_water)

        # Add missing heavy atoms and residues
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        # Target only: the reference never goes through pdb2gmx/pdbxyz, and its
        # 5' phosphates are not used (waters + the CA/C4' fit atoms are).
        stripped_5p = [] if keep_water else strip_five_prime_phosphates(fixer)

        # Kept ions get their own trailing chain (pdb2gmx needs single-type chains).
        ion_chain_id = move_ions_to_own_chain(fixer) if kept_ions else None

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
            print("[DONE] Waters and other heterogens removed"
                  f"{' (force-field ions kept)' if keep_ions else ''}.")
        print("[DONE] Missing heavy atoms and internal residues repaired where possible.")
        if stripped_5p:
            print(f"[DONE] Removed the 5'-terminal phosphate (no AMBER 5'-phosphate residue) "
                  f"from: {', '.join(stripped_5p)}")
        if kept_ions:
            print(f"[DONE] Kept {len(kept_ions)} crystal ion(s), moved to chain {ion_chain_id} "
                  f"(set MWDROP_KEEP_IONS=0 to remove them): {', '.join(kept_ions)}")
        if unified:
            print(f"[DONE] Converted nucleotides to their chain's majority type (DNA/RNA; "
                  f"pdb2gmx needs single-type chains): {', '.join(unified)}")
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
