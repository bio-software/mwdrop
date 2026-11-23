#!/usr/bin/env python3
"""Backfill MobyWat match counts and the SR denominator into a research.sh CSV history.

Jeszenői et al. 2016 (J. Chem. Inf. Model. 56, 148), Eq. 1:

    SR = 100 * (count of matches) / (count of reference water positions)

and its "overall" SR is pooled from the raw counts of all systems (sum of matches /
sum of reference waters), which cannot be recovered from the per-run percentages that
research.sh records. This script adds the raw counts as five columns:

  matches_mer, matches_ida, matches_ide, matches_pos
      Number of matched predicted waters per clustering method = rows flagged "x" in
      the last column of O_system_mt{MER,IDa,IDe,POS}.lst. Each reference water can be
      matched only once per prediction list (MobyWat manual 3.3.3), so this is also the
      number of reference waters that were found.
  sr_denominator
      Size of the reference pool MobyWat divided by: "Waters below B-factor limit" in
      O_system.log (reference waters within dmax of the target, then B <= bmax).

Each CSV row is mapped to the run directory research.sh produced for it
(<OUTPUT_ROOT>/<PDB>/<combo>/<PDB>_l<L>_int<INT>_fin<FIN>/<L>/), and a value is only
written once 100 * matches / denominator reproduces both the SR printed in the .lst
file and the sr_* value already stored in the CSV row. The output folder is only read;
the CSV is backed up before it is rewritten.

Usage:
    add_sr_counts.py <OUTPUT_ROOT> <CSV_FILE> [--dry-run]

Windows paths are accepted (C:\\..., \\\\wsl.localhost\\<distro>\\..., home\\user\\...).
"""

import argparse
import csv
import io
import os
import re
import shutil
import sys
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

NEW_COLUMNS = ["matches_mer", "matches_ida", "matches_ide", "matches_pos", "sr_denominator"]
METHODS = [("MER", "sr_mer", "matches_mer"),
           ("IDa", "sr_ida", "matches_ida"),
           ("IDe", "sr_ide", "matches_ide"),
           ("POS", "sr_pos", "matches_pos")]
KEY_COLUMNS = ["pdb_name", "engine", "layers", "intermediate_md_ns", "final_md_ns"]

SR_TOLERANCE = 0.0011           # SR values are printed with 3 decimals
TIME_TOLERANCE_S = 600          # log finish time vs CSV run_timestamp (after TZ shift)
MAX_TZ_OFFSET_H = 14            # real UTC offsets lie within +/-14 h

LOG_TS_RE = re.compile(r"^\[(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)\]")
MODE_RE = re.compile(r"\]\s*MODE=(\S+)")
RESULTS_RE = re.compile(r"Results are under:\s*(\S+)")
TAG_RE = re.compile(r"^(?P<pdb>.+)_l(?P<l>\d+)_int(?P<int>[0-9.eE+-]+)_fin(?P<fin>[0-9.eE+-]+)$")
MATCH_ROW_RE = re.compile(r"\|\s*x\s*$")
SR_LINE_RE = re.compile(r"^SR\s+([0-9.]+)")
DENOM_RE = re.compile(r"Waters below B-factor limit\s*:\s*(\d+)")


# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
def normalize_path(raw):
    """Turn a Windows / WSL-UNC / slash-less path into a Linux path that exists."""
    p = raw.strip().strip('"').strip("'")
    unc = re.match(r"^[\\/]{2}(?:wsl\.localhost|wsl\$)[\\/][^\\/]+(.*)$", p, re.IGNORECASE)
    if unc:
        p = unc.group(1) or "/"
    drive = re.match(r"^([A-Za-z]):[\\/]?(.*)$", p)
    if drive:
        p = f"/mnt/{drive.group(1).lower()}/{drive.group(2)}"
    p = p.replace("\\", "/")
    candidates = [Path(p).expanduser()]
    if not p.startswith("/"):
        candidates.append(Path("/" + p))
    for c in candidates:
        if c.exists():
            return c.resolve()
    return candidates[0]


def num_key(value):
    """Canonical numeric form so '0.1' == '0.10' and '5' == '5.0'."""
    try:
        return repr(float(value))
    except (TypeError, ValueError):
        return (value or "").strip()


# ----------------------------------------------------------------------------
# Run discovery
# ----------------------------------------------------------------------------
class Run:
    def __init__(self, work_dir):
        self.work_dir = work_dir
        self.pdb = work_dir.parent.name
        self.engine = None
        self.layers = self.int_ns = self.fin_ns = None
        self.analysis_dir = None
        self.finished = None        # naive local datetime from application.LOG
        self.claimed_by = None

    @property
    def key(self):
        return (self.pdb, self.engine, num_key(self.layers), num_key(self.int_ns), num_key(self.fin_ns))

    @property
    def run_dir(self):
        if self.analysis_dir is None:
            return None
        return self.analysis_dir / self.layers if int(self.layers) > 1 else self.analysis_dir

    def label(self):
        return str(self.work_dir)


def parse_application_log(run, log_path):
    last_ts = None
    runtime_ts = None
    results_tag = None
    with open(log_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            m = LOG_TS_RE.match(line)
            if m:
                last_ts = m.group(1)
            if run.engine is None:
                mm = MODE_RE.search(line)
                if mm:
                    run.engine = mm.group(1).strip().lower()
            if "Total runtime was" in line and m:
                runtime_ts = m.group(1)
            rm = RESULTS_RE.search(line)
            if rm:
                # e.g. /other/machine/.../1T2H_l5_int0_fin0.1/{1..5}/ -> keep the tag component
                for part in rm.group(1).split("/"):
                    if TAG_RE.match(part):
                        results_tag = part
    ts = runtime_ts or last_ts
    if ts:
        run.finished = datetime.strptime(ts, "%Y-%m-%d %H:%M:%S")
    return results_tag


def discover_runs(output_root, problems):
    runs = []
    for log_path in sorted(output_root.glob("*/*/application.LOG")):
        run = Run(log_path.parent)
        try:
            tag = parse_application_log(run, log_path)
            tag_dirs = []
            if tag and (run.work_dir / tag).is_dir():
                tag_dirs = [run.work_dir / tag]
            else:
                tag_dirs = [d for d in sorted(run.work_dir.iterdir())
                            if d.is_dir() and TAG_RE.match(d.name)
                            and TAG_RE.match(d.name).group("pdb") == run.pdb]
            if len(tag_dirs) != 1:
                problems.append(f"{run.label()}: expected one <PDB>_l*_int*_fin* dir, found {len(tag_dirs)}; skipped")
                continue
            m = TAG_RE.match(tag_dirs[0].name)
            run.analysis_dir = tag_dirs[0]
            run.layers, run.int_ns, run.fin_ns = m.group("l"), m.group("int"), m.group("fin")
            if not run.engine:
                engine_guess = run.work_dir.name.split("_")[0].lower()
                if engine_guess in ("gromacs", "tinker"):
                    run.engine = engine_guess
                else:
                    problems.append(f"{run.label()}: no MODE= line in application.LOG; skipped")
                    continue
            runs.append(run)
        except Exception as exc:  # one broken run must not stop the batch
            problems.append(f"{run.label()}: could not read run ({type(exc).__name__}: {exc}); skipped")
    return runs


# ----------------------------------------------------------------------------
# Count extraction
# ----------------------------------------------------------------------------
def read_match_list(path):
    """Return (matches, distinct_refs, sr) from an O_system_mt*.lst file."""
    matches = 0
    refs = set()
    sr = None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if MATCH_ROW_RE.search(line):
                matches += 1
                right = line.split("|", 1)[1].split()
                if right:
                    refs.add(right[0])
            else:
                m = SR_LINE_RE.match(line)
                if m:
                    sr = float(m.group(1))
    return matches, len(refs), sr


def extract_counts(run_dir):
    """Return (counts dict, lst SR dict, notes list). Raises ValueError when unusable."""
    if run_dir is None or not run_dir.is_dir():
        raise ValueError(f"final-layer dir missing: {run_dir}")
    counts, lst_sr, notes = {}, {}, []
    for method, _, col in METHODS:
        path = run_dir / f"O_system_mt{method}.lst"
        if not path.is_file():
            raise ValueError(f"no MobyWat result ({path.name} missing)")
        matches, distinct, sr = read_match_list(path)
        if distinct != matches:
            raise ValueError(f"{path.name}: {matches} match rows but {distinct} distinct reference waters")
        if sr is None:
            raise ValueError(f"{path.name}: no SR line")
        counts[col] = matches
        lst_sr[method] = sr

    denom = None
    log_path = run_dir / "O_system.log"
    if log_path.is_file():
        with open(log_path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                m = DENOM_RE.search(line)
                if m:
                    denom = int(m.group(1))
    if denom is None:
        # Prediction-only log: recover the pool size from matches and SR, if all lists agree.
        derived = {round(counts[col] * 100.0 / lst_sr[method])
                   for method, _, col in METHODS if lst_sr[method] > 0}
        if len(derived) != 1:
            raise ValueError("O_system.log has no 'Waters below B-factor limit' line and the "
                             f"denominator cannot be derived consistently from the SRs ({sorted(derived)})")
        denom = derived.pop()
        notes.append("denominator derived from SR (not in O_system.log)")
    if denom <= 0:
        raise ValueError(f"reference pool size is {denom}")
    counts["sr_denominator"] = denom
    return counts, lst_sr, notes


def verify(counts, lst_sr, row):
    """Every 100*m/d must reproduce the .lst SR and the CSV sr_* value."""
    d = counts["sr_denominator"]
    for method, sr_col, col in METHODS:
        calc = 100.0 * counts[col] / d
        if abs(calc - lst_sr[method]) > SR_TOLERANCE:
            return f"{method}: 100*{counts[col]}/{d} = {calc:.3f} != .lst SR {lst_sr[method]:.3f}"
        csv_val = (row.get(sr_col) or "").strip()
        if csv_val:
            try:
                if abs(float(csv_val) - lst_sr[method]) > SR_TOLERANCE:
                    return f"{method}: CSV {sr_col}={csv_val} != .lst SR {lst_sr[method]:.3f} (wrong run mapped?)"
            except ValueError:
                return f"{method}: CSV {sr_col}={csv_val!r} is not a number"
    return None


# ----------------------------------------------------------------------------
# Row <-> run mapping
# ----------------------------------------------------------------------------
def parse_row_time(value):
    try:
        return datetime.strptime((value or "").strip(), "%Y-%m-%dT%H:%M:%SZ")
    except ValueError:
        return None


def learn_offset_hours(rows_by_key, runs_by_key):
    """Whole-hour shift between CSV UTC timestamps and local log times (unambiguous pairs only)."""
    votes = Counter()
    for key, rows in rows_by_key.items():
        runs = runs_by_key.get(key, [])
        if len(rows) == 1 and len(runs) == 1:
            t_row = parse_row_time(rows[0][1].get("run_timestamp"))
            t_run = runs[0].finished
            if t_row and t_run:
                delta = (t_run - t_row).total_seconds()
                hours = round(delta / 3600)
                if abs(hours) <= MAX_TZ_OFFSET_H and abs(delta - hours * 3600) <= TIME_TOLERANCE_S:
                    votes[hours] += 1
    return votes.most_common(1)[0][0] if votes else None


def time_cost(row, run, offset_h):
    t_row = parse_row_time(row.get("run_timestamp"))
    if not t_row or not run.finished:
        return None
    delta = (run.finished - t_row).total_seconds()
    shift = offset_h if offset_h is not None else round(delta / 3600)
    if abs(shift) > MAX_TZ_OFFSET_H:
        return None
    return abs(delta - shift * 3600)


def map_rows(rows_by_key, runs_by_key, offset_h):
    """Return {row_index: (run, note)} and {row_index: reason} for unmapped rows."""
    mapped, unmapped = {}, {}
    for key, rows in rows_by_key.items():
        runs = runs_by_key.get(key, [])
        if not runs:
            for idx, _ in rows:
                unmapped[idx] = "no run directory found for (pdb, engine, layers, int, fin)"
            continue
        if len(rows) == 1 and len(runs) == 1:
            idx, row = rows[0]
            cost = time_cost(row, runs[0], offset_h)
            note = None
            if cost is None or cost > TIME_TOLERANCE_S:
                note = ("application.LOG finish time does not match run_timestamp (e.g. MobyWat re-run "
                        "later); accepted as the only candidate, SR check passed")
            mapped[idx] = (runs[0], note)
            continue
        pairs = []
        for idx, row in rows:
            for run in runs:
                cost = time_cost(row, run, offset_h)
                if cost is not None and cost <= TIME_TOLERANCE_S:
                    pairs.append((cost, idx, run))
        used_rows, used_runs = set(), set()
        for cost, idx, run in sorted(pairs, key=lambda p: p[0]):
            if idx in used_rows or id(run) in used_runs:
                continue
            mapped[idx] = (run, None)
            used_rows.add(idx)
            used_runs.add(id(run))
        for idx, _ in rows:
            if idx not in used_rows:
                unmapped[idx] = (f"{len(runs)} candidate run dirs, none finished at run_timestamp "
                                 f"(+/- {TIME_TOLERANCE_S}s); ambiguous, left blank")
    return mapped, unmapped


# ----------------------------------------------------------------------------
# CSV I/O
# ----------------------------------------------------------------------------
def read_csv(path):
    raw = path.read_bytes()
    bom = raw.startswith(b"\xef\xbb\xbf")
    text = raw.decode("utf-8-sig")
    newline = "\r\n" if "\r\n" in text else "\n"
    records = list(csv.reader(io.StringIO(text, newline="")))
    if not records or not records[0]:
        raise ValueError("CSV is empty or has no header")
    return records[0], records[1:], newline, bom


def write_csv(path, header, body, newline, bom):
    buf = io.StringIO(newline="")
    writer = csv.writer(buf, lineterminator=newline)
    writer.writerow(header)
    writer.writerows(body)
    data = buf.getvalue().encode("utf-8")
    if bom:
        data = b"\xef\xbb\xbf" + data
    fd, tmp = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent))
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(data)
        shutil.copymode(path, tmp)
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("output_root", help="research.sh OUTPUT_DIR (holds <PDB>/<combo>/ run dirs)")
    ap.add_argument("csv_file", help="research.sh CSV history to update")
    ap.add_argument("--dry-run", action="store_true", help="report only; do not back up or write")
    args = ap.parse_args(argv)

    output_root = normalize_path(args.output_root)
    csv_path = normalize_path(args.csv_file)
    fatal = []
    if not output_root.is_dir():
        fatal.append(f"output folder not found: {args.output_root} -> {output_root}")
    if not csv_path.is_file():
        fatal.append(f"CSV file not found: {args.csv_file} -> {csv_path}")
    if fatal:
        for f in fatal:
            print(f"ERROR: {f}", file=sys.stderr)
        return 2

    try:
        header, body, newline, bom = read_csv(csv_path)
    except Exception as exc:
        print(f"ERROR: cannot read {csv_path}: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2

    missing = [c for c in KEY_COLUMNS + [s for _, s, _ in METHODS] if c not in header]
    if missing:
        print(f"ERROR: CSV header lacks required column(s): {', '.join(missing)}", file=sys.stderr)
        return 2

    orig_width = len(header)
    out_header = header + [c for c in NEW_COLUMNS if c not in header]
    col_pos = {c: out_header.index(c) for c in NEW_COLUMNS}

    problems = []                 # run-discovery problems
    row_status = {}               # idx -> (status, detail)
    rows_by_key = defaultdict(list)
    for idx, rec in enumerate(body):
        if not rec or all(not f.strip() for f in rec):
            continue
        if len(rec) != orig_width:
            row_status[idx] = ("malformed", f"{len(rec)} fields, header has {orig_width}; left unchanged")
            continue
        row = dict(zip(header, rec))
        key = (row["pdb_name"].strip(), row["engine"].strip().lower(), num_key(row["layers"]),
               num_key(row["intermediate_md_ns"]), num_key(row["final_md_ns"]))
        if not all(key[:2]):
            row_status[idx] = ("malformed", "empty pdb_name/engine; left unchanged")
            continue
        rows_by_key[key].append((idx, row))

    runs = discover_runs(output_root, problems)
    runs_by_key = defaultdict(list)
    for run in runs:
        runs_by_key[run.key].append(run)

    offset_h = learn_offset_hours(rows_by_key, runs_by_key)
    mapped, unmapped = map_rows(rows_by_key, runs_by_key, offset_h)
    for idx, reason in unmapped.items():
        row_status[idx] = ("unmatched", reason)

    new_values = {}
    for idx, (run, note) in sorted(mapped.items()):
        run.claimed_by = idx
        row = dict(zip(header, body[idx]))
        csv_srs = [(row.get(s) or "").strip() for _, s, _ in METHODS]
        try:
            counts, lst_sr, notes = extract_counts(run.run_dir)
        except Exception as exc:
            if not any(csv_srs):
                row_status[idx] = ("no-result", f"failed run (no SR in CSV): {exc}; {run.label()}")
            else:
                row_status[idx] = ("error", f"{exc}; {run.label()}")
            continue
        if not any(csv_srs):
            row_status[idx] = ("verify-failed", "CSV row has no SR values but the run dir has MobyWat "
                               f"results; refusing to guess ({run.label()})")
            continue
        bad = verify(counts, lst_sr, row)
        if bad:
            row_status[idx] = ("verify-failed", f"{bad}; {run.label()}")
            continue
        new_values[idx] = counts
        extra = [n for n in [note] + notes if n]
        row_status[idx] = ("ok", "; ".join(extra))

    # Build output rows: new columns appended (or refilled in place on a re-run).
    out_body = []
    for idx, rec in enumerate(body):
        status = row_status.get(idx, (None, None))[0]
        if not rec or status == "malformed":
            out_body.append(rec)
            continue
        new = list(rec) + [""] * (len(out_header) - len(rec))
        if idx in new_values:
            for col in NEW_COLUMNS:
                new[col_pos[col]] = str(new_values[idx][col])
        out_body.append(new)

    # ---------------------------------------------------------------- report
    def row_label(idx):
        rec = body[idx]
        r = dict(zip(header, rec)) if len(rec) == orig_width else {}
        return (f"line {idx + 2} id={r.get('id', '?')} {r.get('pdb_name', '?')}/{r.get('engine', '?')} "
                f"l{r.get('layers', '?')} int{r.get('intermediate_md_ns', '?')} fin{r.get('final_md_ns', '?')}")

    tally = Counter(s for s, _ in row_status.values())
    print("=" * 78)
    print(f"add_sr_counts.py {'(DRY RUN) ' if args.dry_run else ''}report")
    print(f"  output folder : {output_root}")
    print(f"  CSV           : {csv_path}")
    print(f"  run dirs found: {len(runs)}   data rows: {sum(len(v) for v in rows_by_key.values()) + tally['malformed']}")
    print(f"  TZ offset     : {'+' if (offset_h or 0) >= 0 else ''}{offset_h} h (log local time vs CSV UTC)"
          if offset_h is not None else "  TZ offset     : not determined")
    print("-" * 78)
    print("WHAT WENT RIGHT")
    print(f"  {tally['ok']} row(s) filled and verified (100*matches/denominator == .lst SR == CSV sr_*)")
    for idx in sorted(i for i, (s, d) in row_status.items() if s == "ok"):
        v = new_values[idx]
        line = (f"    {row_label(idx)}: " + ",".join(str(v[c]) for c in NEW_COLUMNS))
        print(line + (f"  [{row_status[idx][1]}]" if row_status[idx][1] else ""))
    print("-" * 78)
    print("WHAT WENT WRONG / LEFT BLANK")
    labels = {"no-result": "failed run, no MobyWat output (blank, expected)",
              "unmatched": "no matching run directory (blank)",
              "verify-failed": "counts did not reproduce SR (blank)",
              "error": "could not extract counts (blank)",
              "malformed": "malformed CSV row (unchanged)"}
    any_wrong = False
    for status, title in labels.items():
        idxs = sorted(i for i, (s, _) in row_status.items() if s == status)
        if not idxs:
            continue
        any_wrong = True
        print(f"  {title}: {len(idxs)}")
        for idx in idxs:
            print(f"    {row_label(idx)}: {row_status[idx][1]}")
    unclaimed = [r for r in runs if r.claimed_by is None]
    if unclaimed:
        any_wrong = True
        print(f"  run dirs not claimed by any CSV row: {len(unclaimed)}")
        for r in unclaimed:
            print(f"    {r.label()} ({r.engine}, l{r.layers} int{r.int_ns} fin{r.fin_ns})")
    if problems:
        any_wrong = True
        print(f"  run dirs that could not be read: {len(problems)}")
        for p in problems:
            print(f"    {p}")
    if not any_wrong:
        print("  nothing")
    print("-" * 78)

    if args.dry_run:
        print("Dry run: no backup made, CSV not modified.")
        return 0

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    backup = csv_path.with_name(f"{csv_path.name}.bak-{stamp}")
    try:
        shutil.copy2(csv_path, backup)
    except Exception as exc:
        print(f"ERROR: backup to {backup} failed ({type(exc).__name__}: {exc}); CSV NOT modified.", file=sys.stderr)
        return 3
    print(f"Backup written : {backup}")
    try:
        write_csv(csv_path, out_header, out_body, newline, bom)
    except PermissionError as exc:
        print(f"ERROR: cannot write {csv_path} ({exc}). Is it open in Excel? Original left untouched; "
              f"backup at {backup}.", file=sys.stderr)
        return 3
    except Exception as exc:
        print(f"ERROR: writing {csv_path} failed ({type(exc).__name__}: {exc}); original left untouched; "
              f"backup at {backup}.", file=sys.stderr)
        return 3
    added = [c for c in NEW_COLUMNS if c not in header]
    print(f"CSV updated    : {csv_path}")
    print(f"  columns {'added: ' + ', '.join(added) if added else 'already present, refilled in place'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
