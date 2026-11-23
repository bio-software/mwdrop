# Nucleic-acid compatibility of the MWDROP pipeline

Date: 2026-11-23 · Test system: **1VZK** (B-DNA dodecamer, 2 chains × 12 nt, 758 atoms with H; the
crystal also holds ligand `D1B` + one Mg²⁺, both stripped by the preprocessor as for proteins).
Commit `a218381` plus the uncommitted `tc-grps` fix described in §2 (step G5).
MobyWat source findings: `/home/hangb/mobywat/MOBYWAT_NUCLEIC_ACIDS_REPORT.md` (MobyWat 1.1,
source review + synthetic-DNA tests), cross-checked here on the real 1VZK trajectory (§4).

## Status after fixes (2026-11-23, later the same day)

**MWDROP now runs DNA and RNA end to end on both engines.** Fixed in the working tree (uncommitted):
`pipeline/pipeline_common.sh`, `gromacs.sh`, `tinker.sh`, `research.sh`, `CLAUDE.md`, plus the
patched MobyWat from https://github.com/bio-software/MobyWat (commit `c4c1d3b`, installed as
`~/bin/mobywat`; previous binary kept as `~/bin/mobywat.bak-20261123`).

Test runs (5 layers, final MD 0.01 ns, `MOBYWAT_DEBUG=1`; GROMACS on CPU, Tinker with `TINKER_GPU=1`):

| System | Engine | Result | MobyWat target / RMSD atom | SR MER / IDa / IDe / POS | Ref. waters | Per-frame SR mean (1000 frames) | RMSD frame 2 |
|---|---|---|---|---|---|---|---|
| 1VZK DNA, 2 chains | GROMACS | ✅ 300 s | `[AB]` / `C4'` | 69.388 / 63.265 / 59.184 / 69.388 | 49 | 51.441 | 0.026 Å |
| 1VZK DNA, 2 chains | Tinker GPU | ✅ 1081 s | `[AB]` / `C4'` | 79.592 / 85.714 / 83.673 / 77.551 | 49 | 63.145 | 0.052 Å |
| 483D RNA, 1 chain | GROMACS | ✅ 325 s | `[A]` / `C4'` | 69.048 / 69.841 / 69.048 / 69.048 | 126 | 56.961 | 0.028 Å |
| 483D RNA, 1 chain | Tinker GPU | ✅ 1187 s | `[A]` / `C4'` | 76.984 / 77.778 / 76.190 / 76.984 | 126 | 66.036 | 0.052 Å |

Runtimes are with all four runs in parallel; most of the Tinker time is MobyWat reading the
1000-model PDB trajectory twice (Prediction + Analysis).

- `research.sh` (Tinker GPU, 1VZK) writes a complete row: `rmsd_avg` 0.046, `rmsd_max` 0.058,
  all SRs, `sr_frame_*`, `matches_*` 39/42/41/38, `sr_denominator` 49, **`failure=0`**.
- Protein regression (1R6J, `testing/fast/gromacs.sh`): SRs 53.804 / 57.609 / 58.152 / 52.717 and the
  RMSD column identical to before the changes; target `[A]`, RMSD atom `CA`.
- Patched vs original MobyWat on proteins: all 28 output files byte-identical; CPU time unchanged.

Fix status per issue is in the §6 table. Still open: #7 (wdrop, separate repo) and the reference
B-factor problem described at the end of §6.

**Bottom line (original analysis, before the fixes):**

- **MobyWat itself supports nucleic acids as-is.** It has no residue tables, never superimposes,
  and selects the target by chain. Run on the real 1VZK Tinker trajectory with `-t "[AB]"` it gives
  SR MER 81.633 / IDa 81.633 / IDe 79.592 / POS 85.714 and a full per-frame Analysis (§4).
- **Everything up to and including MD works on DNA for both engines.**
- **What breaks is MWDROP's own glue before MobyWat**, which assumes a protein: the chain-ID
  detection keys on `CA` atoms (both engines die silently), and GROMACS fits the trajectory on
  numeric group `3` (= C-alpha for proteins, = water for DNA).
- Proteins are unaffected by everything listed here.

Legend: ✅ verified on 1VZK · ⚠️ works but with a caveat · ❌ fails on nucleic acids ·
⏳ not reached through the pipeline

---

## 1. Shared steps (`mwdrop.sh`, `pipeline_common.sh`)

| Step | What it does | DNA | Evidence / notes |
|---|---|---|---|
| 0B `pdb-preprocessor.py --target` | pdbfixer: missing atoms/residues, strip waters + heterogens | ✅ | Ran clean on 1VZK. Removes *all* heterogens, so drug/ion ligands bound to DNA (here `D1B`, Mg²⁺) are lost — same behaviour as for protein ligands. |
| wdrop | Deposits one HOH layer per iteration | ⚠️ | Works (134 waters in layer 1). Radii are element-based (`sigma_from_atom` in wdrop `src/wdrop.c`); **P has no entry and falls back to 1.70 Å** instead of the usual 1.80 Å. 22 P atoms in 1VZK — minor geometric effect on the backbone. |
| `format_pdb.py` | Fixed-width rewrite, water renaming | ✅ | Primed atom names (`O5'`, `C1'`) and `DA/DC/DG/DT` survive unchanged through both engines. |
| `protein_chain_ids()` / `mobywat_target_spec()` | Chain IDs of `CA`-bearing residues → MobyWat `-t [..]` | ❌ | DNA has 0 `CA` atoms → empty string. `mobywat_target_spec` then returns 1 (`[[ -n "" ]] && printf`), and under `set -e` the assignment `MOBYWAT_TARGET="$(mobywat_target_spec …)"` **kills the script with rc=1 and no log message**. Reproduced directly and in both engine runs. (Had it not died, MobyWat would have **segfaulted** on `-t ""` / `-t "[]"` / an empty target REMARK — MobyWat report Q5.) |

## 2. `gromacs.sh`

| # | Step | DNA | Evidence / notes |
|---|---|---|---|
| G0 | `format_pdb.py` | ✅ | |
| G1 | `pdb2gmx -ff amber99sb-ildn -water tip3p -ignh` | ✅ | Builds `DNA_chain_A/B` topologies with the force field's own DNA residue definitions. Total charge −22 e (−11 per strand) — expected for DNA; no counter-ions in the vacuum setup, worth knowing for physical interpretation. |
| G1 | Position restraints | ✅ | `posre_DNA_chain_{A,B}.itp`, included under `#ifdef POSRES`; 486 restrained heavy atoms = all DNA heavy atoms, 1000 kJ mol⁻¹ nm⁻² — identical to protein handling. |
| G2 | `editconf -c -box 1000` | ✅ | Note: this **translates the system by ~5000 Å** (solute centre 15.1/20.8/8.6 Å in the crystal → 4999.7/4999.9/4999.3 Å in `md.tpr`). Harmless for MD, but it is why the fit in G8 is mandatory (§4). |
| G3 | `grompp`/`mdrun` steepest descent | ✅ | |
| G4 | `grompp`/`mdrun` CG | ✅ | |
| G5 | `grompp` MD (`gromacs-md.mdp`) | ✅ (after fix) | **Was ❌**: `tc-grps`/`energygrps = Protein non-Protein` → *"Group Protein … not found"*. Changed to `non-Water Water` (also in `determinism_lib.sh`). 1R6J protein control gives identical SRs before/after (53.804 / 57.609 / 58.152 / 52.717). |
| G5 | `mdrun` MD | ✅ | 0.01 ns run completed. |
| G6 | `trjconv` pbc whole / compact / last frame (groups `0`, `1 0`) | ✅ | Group 1 is `DNA` here, `Protein` for proteins — both "the solute", so the numeric selection happens to work. |
| G8 | `pdb-preprocessor.py --reference` | ✅ | |
| G8 | `trjconv` pbc1/pbc2/pbc3 (groups `0`, `1 0`) | ✅ | |
| G8 | `confrms -one … <<< "3 3"` | ❌ | **First failure.** Group 3 is `C-alpha` for a protein but `SOL` for DNA (default groups: 0 System, 1 DNA, 2 Water, 3 SOL, 4 non-Water). Compares reference waters vs simulated waters → *"You selected groups with different number of atoms."* |
| G8 | `trjconv -fit progressive <<< "3 0"` (×2) | ❌ | Same group-3 problem. This fit is what puts the frames back into the crystal frame, and MobyWat does no fitting of its own, so it must be fixed, not skipped. |
| G8 | `ensure_chain_labels` | ❌ | Uses `protein_chain_ids` → always empty for DNA → stamps chain `A` on **every** atom (`editconf -label A`), merging both strands. Not needed for DNA: `md.tpr` already carries chains A/B (379 + 379 atoms, waters blank). |
| G8 | `mobywat_target_spec` | ❌ | Silent `set -e` exit (§1). |
| G8 | MobyWat Prediction / Analysis | ⏳ | Not reached through the pipeline. MobyWat's xtc + topology path is the same code as the Tinker path after conversion (MobyWat report, "Not verified"), chain IDs A/B are present, so it is expected to work once G8's fit and target are fixed. |

Prediction-only branch (no reference): same `trjconv -fit` group 3 + `ensure_chain_labels` +
`mobywat_target_spec` issues.

## 3. `tinker.sh`

| # | Step | DNA | Evidence / notes |
|---|---|---|---|
| T1 | `format_pdb.py` | ✅ | |
| T3 | `pdbxyz` with `amber99.prm` | ✅ | 1160 atoms, **0 untyped (type 0) atoms**; 758 solute atoms — same count as GROMACS. |
| T4 | Restraint atom count (atoms before first `OW`/`HW`) | ✅ | 758 — solute-agnostic logic, covers DNA incl. H. (Log text still says "protein atoms".) |
| T5 | `minimize` (`RESTRAIN-POSITION -1 N 2.0`) | ✅ | All 5 iterations: *"LBFGS — Normal Termination due to SmallGrad"*. |
| T6 | `dynamic` MD (`RESTRAIN-POSITION -1 N 300.0`) | ✅ | 0.01 ns, 1000 frames, rc=0 (6 min 40 s on CPU). Frames stay in the crystal coordinate frame (no recentring), which is what MobyWat needs. |
| T7–T11 | `arcedit`, `xyzpdb`, multi-model PDB, `format_pdb.py` | ✅ | All rc=0; `system_mdl.pdb` (1000 models, chains A/B, waters blank chain) written. |
| T12 | `pdb-preprocessor.py --reference` | ✅ | `system_ref.pdb`: chains A/B, 60 HOH, ligand and Mg²⁺ removed, no `CA`-prefixed atom names. |
| T12 | `mobywat_target_spec` (trajectory + reference) | ❌ | **Failure point, confirmed in the run:** the log ends right after `pdb-preprocessor.py --reference` with no error line; `tinker.sh` exits rc=1 (§1). `research.sh` records `failure=1`. |
| T12 | MobyWat Prediction + Analysis | ✅ (manual) | Run by hand on the run's own files with the target the pipeline *should* have built — see §4. |

Tinker needs no group selections or fitting, so its **only** blocker is the target selection.
End-to-end runtime for 1VZK (5 layers, 0.01 ns): 711 s.

## 4. MobyWat on nucleic acids

### What the source review found (summary of the MobyWat report)

- **No protein-specific logic.** No residue-name tables, no amino-acid filter. The only residue names
  it knows are waters (`WAT`, `SOL`, `H2O`, `HOH`). Primed and 4-character atom names are stored
  verbatim; hydrogens are recognised by a leading `H`, which is correct for nucleic acids.
- **Target selection** (`-t` and `REMARK mobywat_reference_target`): `[AB]` = every non-water atom
  of chains A and B, whatever the residue type. Note that ions/ligands carrying chain A/B would be
  included.
- **No superposition anywhere.** Frames and the reference are compared in raw coordinates, so the
  pipeline must deliver frames already in the crystal frame.
- **`CA` is hard-coded only in the target RMSD** (`O_system_rmst.txt`). With 0 CA atoms on both
  sides it does not abort; it writes `-nan`. Analysis aborts (exit 0, no `O_system_succ.txt`) if
  the reference and trajectory have *different, non-zero* counts of atoms whose names start with
  `CA` — e.g. a crystal Ca²⁺ in chain A/B that is missing from the MD system.
- **Exit status is unreliable**: many errors `exit(0)`; empty target specs **segfault** (exit 139).

### Confirmed on the real 1VZK Tinker output

Files: the run's own `system_mdl.pdb` (1000 models) and `system_ref.pdb`, with the REMARKs written by
`add-mobywat-analysis-params.sh system_ref.pdb "[AB]"` + `remove-ter-id.sh`, i.e. exactly what the
pipeline would produce with a nucleic-acid-aware target spec. Commands are the pipeline's own, with
`-t "[AB]"`.

| Pass | Exit | Result |
|---|---|---|
| Prediction (`-cls MER -v Diagnostic`) | 0, no `ERROR` line | SR **MER 81.633, IDa 81.633, IDe 79.592, POS 85.714**; matches 40/40/39/42 of **49** reference waters (100 × 40 / 49 = 81.633 ✓ — `research.sh`'s `matches_*`/`sr_denominator` logic holds). Log: `Target CA's : 0`. |
| Analysis | 0, no `ERROR` line | `O_system_succ.txt`, 1000 frames: min 51.020, max 81.633, mean 67.949, std 5.698 |
| `O_system_rmst.txt` | — | every row `-nan` |

The high SRs also confirm that the Tinker frames are in the crystal frame without any fitting.

### Consequence for `research.sh`

`research.sh` averages `O_system_rmst.txt` with awk, which turns all-`-nan` input into the string
`-nan`. So a DNA row would get `rmsd_avg=-nan`, `rmsd_max=-nan` — **non-empty, so `failure` stays 0**
— and the CSV gets a non-numeric value that will break later numeric import. Needs handling (§5).

## 5. `research.sh` / `add_sr_counts.py`

Otherwise engine- and molecule-independent. A DNA run currently records a row with `failure=1` and
empty metrics because the backend dies before MobyWat.

## 6. What needs to change (not implemented), by severity

Severity scale: **Critical** = DNA runs cannot complete · **High** = a run completes but records
wrong or misleading data without warning · **Medium** = harder to diagnose, or affects only some
inputs · **Low** = minor accuracy or cosmetic.

| # | Severity | Issue | Engines | Affects proteins? | Status |
|---|---|---|---|---|---|
| 1 | Critical | Target chains detected from `CA` atoms → empty for DNA → silent exit | both | no | ✅ fixed: `target_chain_ids()` (non-water atoms) |
| 2 | Critical | GROMACS fits on numeric group `3` (= water for DNA) → `confrms` fails | GROMACS | no | ✅ fixed: named `Fit` group (CA + C4′) via `write_fit_index()` |
| 3 | High | `-nan` RMSD recorded as a value, `failure` stays 0 | both | no | ✅ fixed: `MOBYWAT_RMS_ATOM=C4'` gives real RMSD; `research.sh` skips nan |
| 4 | High | MobyWat errors exit 0 → stale or missing results accepted as success | both | **yes** (pre-existing) | ✅ fixed in MobyWat (exit 1 on errors) |
| 5 | Medium | Empty target is never reported (silent exit; would segfault MobyWat) | both | no | ✅ fixed: logged `ERROR`, non-zero return |
| 6 | Low | A Ca²⁺ named `CA` in the reference would abort Analysis | both | yes (latent) | ✅ fixed in MobyWat (ions named like the atom are skipped) |
| 7 | Low | wdrop uses 1.70 Å for P instead of 1.80 Å | both | no | open (wdrop repo) |
| 8 | Low | Protein-only wording in logs and comments | both | — | ✅ fixed |

1. **Critical: target chain detection** (`protein_chain_ids` in `pipeline_common.sh`, used by
   `mobywat_target_spec` and `ensure_chain_labels`). Every DNA run stops before MobyWat, in both
   engines, and `ensure_chain_labels` would merge both strands into chain A. Fix: take the chain IDs
   of **all non-water ATOM/HETATM records** (or of residues with `CA` **or** `C4'`/`P`). For 1VZK this
   gives `[AB]`, verified to work (§4). Proteins resolve to the same chains as today, since the
   preprocessor already strips heterogens. **This alone unblocks Tinker.**
2. **Critical: GROMACS fit group** (`confrms`, `trjconv -fit` ×2 in each branch). Group `3` is water
   for DNA, so `confrms` fails. The fit cannot simply be skipped: `editconf -c` moves the system
   ~5000 Å (§2 G2) and MobyWat compares raw coordinates, so unfitted frames would score near zero.
   Fix: build an index with `gmx make_ndx`/`gmx select` containing protein `CA` **plus**
   nucleic-acid `C4'`, and pass `-n`. The MobyWat report recommends `C4'`/`C1'` over `P`: one per
   nucleotide and present at the 5′ terminus, where crystal structures often lack the phosphate.
   Numeric groups `0` and `1` work today only because of default group ordering; consider naming
   them too. **Together with #1, this unblocks GROMACS.**
3. **High: `-nan` RMSD in `research.sh`.** Certain to hit every DNA row once #1 and #2 are fixed:
   `rmsd_avg`/`rmsd_max` become the string `-nan`, so `failure` stays 0 and the CSV gets non-numeric
   values that will break the later database import. Fix: treat `nan`/`-nan` as missing (empty). Then
   decide whether a DNA run with valid SRs but no RMSD counts as `failure=1` (current definition) or
   whether `failure` should look only at the SRs. Alternatively apply the optional MobyWat patch (env
   `MOBYWAT_RMS_ATOM`, e.g. `C4'`) so the RMSD is real for DNA.
4. **High: MobyWat failures are not detected** (pre-existing, also affects proteins). MobyWat often
   exits 0 after an error, and an aborted pass leaves earlier output files in place, so the pipeline
   can record stale or missing results as a success. Ranked below #3 because it only bites when
   MobyWat errors, whereas #3 hits every DNA row. Fix: treat a stdout line starting with `ERROR`, a
   non-zero exit, or a missing/stale `O_system_mt*.lst` / `O_system_succ.txt` as a failed pass.
5. **Medium: empty target is never reported.** Today it causes the silent `set -e` exit, so the only
   trace is `failure=1` with no reason in the log. If the empty value ever reached MobyWat, it would
   segfault (exit 139). Mostly defused by #1, but still possible for odd inputs (e.g. no chain IDs at
   all). Fix: log `ERROR: no target chains found in <file>` and exit non-zero explicitly.
6. **Low: ions named `CA`.** If a future reference keeps a Ca²⁺ in chain A/B, Analysis silently
   aborts on MobyWat's CA-count check. Latent only: the preprocessor currently removes all
   heterogens.
7. **Low: wdrop P radius** (separate repo). Add `P` (1.80 Å) to `sigma_from_atom`; today P gets
   1.70 Å, a small geometric effect on backbone water placement.
8. **Low: wording.** "protein restraints / protein atom count" in `tinker.sh` logs; "C-alpha RMSD"
   in `research.sh` comments.

**Found while fixing, not yet addressed — reference B-factors are lost.** `pdb-preprocessor.py
--reference` rewrites `system_ref.pdb` through pdbfixer/OpenMM, which writes every B-factor as
0.00. MobyWat's B-factor limit (`-bmax`, 30 Å²) therefore never excludes a reference water: in
1VZK 40 of the 60 crystal waters have B ≥ 30, in 1R6J 7 of 237. This affects proteins as well and
every SR recorded so far. Fixing it (carrying the crystal B-factors through the preprocessor) will
change SRs and make new rows not comparable with the existing history, so it needs a decision.

**Cosmetic, MobyWat:** its log still labels the RMSD-atom count "Target CA's" (it shows 24 `C4'`
for 1VZK).

## 7. How this was tested

- `research.sh` on `nucleicacid/1123_test/input/1VZK` (Tinker combo, 5 layers, 0.01 ns final MD).
- `mwdrop.sh 1VZK.pdb gromacs 1VZK_r.pdb --layers 5 --final-md-ns 0.01` with `MOBYWAT_DEBUG=1`.
- `mobywat_target_spec` called directly on a Tinker-produced DNA PDB under `set -euo pipefail`.
- MobyWat Prediction + Analysis run by hand on copies of the Tinker run's `system_mdl.pdb` /
  `system_ref.pdb` with `-t "[AB]"` (§4).
- GROMACS chain IDs and coordinate offset read from `md.tpr` via `gmx editconf`.
- Protein regression control: `testing/fast/gromacs.sh` (1R6J) after the `tc-grps` change.
- All test outputs were written to a scratch directory; the repository's CSV history was not touched.
