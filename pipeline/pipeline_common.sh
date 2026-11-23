#!/bin/bash
set -euo pipefail

normalize_input_pdb() {
  local input="$1"

  if [[ ! -f "$input" ]]; then
    echo "Error: file not found: $input" >&2
    return 1
  fi

  if [[ "$input" != *.pdb ]]; then
    echo "Error: input must end with .pdb (got: $input)" >&2
    return 1
  fi

  printf '%s\n' "$(cd -P "$(dirname "$input")" && pwd)/$(basename "$input")"
}

make_output_dir() {
  local input_abs="$1"
  local suffix="$2"
  local base_name
  base_name="$(basename "${input_abs%.pdb}")"
  printf '%s/%s%s\n' "$(dirname "$input_abs")" "$base_name" "$suffix"
}

validate_script_dir_not_input_dir() {
  local input_path="$1"
  local input_dir_candidate
  local input_dir_real
  local SCRIPT_DIR_VAR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

  input_dir_candidate="$(dirname -- "$input_path")"

  if [[ -d "$input_dir_candidate" ]]; then
    input_dir_real="$(cd -P "$input_dir_candidate" && pwd)"

    if [[ "$SCRIPT_DIR_VAR" == "$input_dir_real" ]]; then
      log "Error: mwdrop.sh cannot be located in the same directory as INPUT_PDB." >&2
      log "Move the input PDB to a different folder!" >&2
      exit 1
    fi
  fi
}

# ============================================================================
# LOGGING FUNCTIONS
# ============================================================================
# These functions manage unified logging across all pipeline scripts.
# Global variable: LOGFILE (set by calling scripts)

setup_logging() {
  local logfile="$1"
  export LOGFILE="$logfile"

  # Initialize log file only if it doesn't exist (don't truncate existing files)
  if [[ ! -f "$LOGFILE" ]]; then
    : > "$LOGFILE"
  fi
}

log() {
  local line
  if [[ -z "${LOGFILE:-}" ]]; then
    while IFS= read -r line; do
      printf '[%(%F %T)T] %s\n' -1 "$line" >&2
    done <<< "$*"
  else
    while IFS= read -r line; do
      printf '[%(%F %T)T] %s\n' -1 "$line"
    done <<< "$*" | tee -a "$LOGFILE"
  fi
}

run_step() {
  log "RUN: $*"
  "$@"
  local rc=$?
  log "DONE (rc=$rc): $*"
  return "$rc"
}

# ============================================================================

# MD runs only when a positive duration is set. Shared by the MM backends
# (gromacs.sh, tinker.sh) to gate the MD simulation step.
md_enabled() { awk -v D="${1:-0}" 'BEGIN { exit !(D+0 > 0) }'; }

# MobyWat runs only on the final iteration, and only when MD produced a trajectory.
# $1 = md_duration for this cycle, $2 = run-mobywat flag (1/0; empty defaults to 1
# for backward compatibility with direct backend callers). MobyWat needs an MD
# trajectory, so a non-positive duration disables it regardless of the flag.
should_run_mobywat() {
  md_enabled "${1:-0}" || return 1
  case "${2:-1}" in
    ''|0|false|no|off) return 1 ;;
    *)                 return 0 ;;
  esac
}

# MobyWat picks its target molecule by chain ID (-t x-y/[xy...]), both on the
# command line (target inside the trajectory) and via REMARK mobywat_reference_target
# (target inside the reference). Both used to be hardcoded to [A], which silently
# validated only the first chain of every multi-chain protein.

# Sorted, unique chain IDs of the non-water atoms of a PDB, e.g. "AB". This is the
# set MobyWat's chain mode selects (every non-water atom of the listed chains,
# whatever the residue type), so proteins, DNA and RNA all resolve; the water names
# are MobyWat's own (HOH/WAT/SOL/H2O). Heterogens are stripped by
# pdb-preprocessor.py beforehand, so for a protein this is the same set as its
# CA-bearing chains. Only the first MODEL is read, so multi-model trajectory PDBs
# work unchanged.
target_chain_ids() {
  awk '
    /^ENDMDL/ { exit }
    /^(ATOM|HETATM)/ {
      res = substr($0, 18, 3); gsub(/ /, "", res)
      id = substr($0, 22, 1)
      if (res != "HOH" && res != "WAT" && res != "SOL" && res != "H2O" && id != " ") print id
    }' "$1" | sort -u | tr -d '\n'
}

# Target spec of a PDB, e.g. "[AB]". Derived per file: the trajectory and the
# reference may spell their chains differently (MobyWat manual 3.1.5.2).
# An empty spec would crash MobyWat, so it is a logged error (on stderr, because
# callers capture stdout).
mobywat_target_spec() {
  local ids
  ids="$(target_chain_ids "$1")"
  if [[ -z "$ids" ]]; then
    log "ERROR: no target chains (non-water atoms with a chain ID) found in $1" >&2
    return 1
  fi
  printf '[%s]' "$ids"
}

# Atom name MobyWat should use for the target RMSD (O_system_rmst.txt), passed via
# the MOBYWAT_RMS_ATOM environment variable (patched MobyWat, bio-software/MobyWat):
# CA when the target has C-alphas (proteins and complexes), else C4' (DNA/RNA: one
# per nucleotide, present at the 5' terminus where crystals often lack the
# phosphate). Defaults to CA, which MobyWat reports as nan with a warning.
mobywat_rms_atom() {
  awk '
    /^ENDMDL/ { exit }
    /^(ATOM|HETATM)/ {
      name = substr($0, 13, 4); gsub(/ /, "", name)
      res = substr($0, 18, 3); gsub(/ /, "", res)
      if (name == "CA" && res != "CA") ca = 1
      if (name == "C4'"'"'") c4 = 1
    }
    END { print ((ca || !c4) ? "CA" : "C4'"'"'") }' "$1"
}

run_mm_step() {
  local mode="$1"
  local input_pdb="$2"
  local script_dir="$3"
  local md_duration="${4:-1}" # in nanosec, default 1; >0 runs MD, 0 skips it
  local REFERENCE_PDB_COM="${5}"
  local run_mobywat="${6:-1}" # 1 on the final iteration (run MobyWat), 0 otherwise

  case "$mode" in
    gromacs)
      "$script_dir/gromacs.sh" "$input_pdb" "$md_duration" 1000 "$REFERENCE_PDB_COM" "$run_mobywat"
      ;;
    tinker)
      "$script_dir/tinker.sh" "$input_pdb" "$md_duration" 1000 "$REFERENCE_PDB_COM" "$run_mobywat"
      ;;
    *)
      echo "Error: unsupported mode: $mode (use gromacs or tinker)" >&2
      return 1
      ;;
  esac
}