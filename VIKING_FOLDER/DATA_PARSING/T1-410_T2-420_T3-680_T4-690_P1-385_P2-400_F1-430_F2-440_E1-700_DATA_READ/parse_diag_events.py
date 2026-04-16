#!/usr/bin/env python3
"""
parse_diag_events.py
====================
Reads the stdout from KLong_save_vectors.C (with DEBUG_MODE=true and
DIAG_RECO_MODE=true) piped in on stdin, and prints a clean summary of
every event that passes all selection criteria.

Usage
-----
    root -l -b -q 'KLong_save_vectors.C' 2>&1 | python3 parse_diag_events.py
    root -l -b -q 'KLong_save_vectors.C' 2>&1 | python3 parse_diag_events.py --all
    root -l -b -q 'KLong_save_vectors.C' 2>&1 | python3 parse_diag_events.py --log diag_clean.txt

Options
-------
    --all           Print all events (including failures); failures show which
                    criteria they fail
    --log FILE      Tee the clean-event summary to FILE as well as stdout
    --poca  VAL     PoCA separation threshold in cm  (default 1.0)
    --dt    VAL     |time diff| threshold in ns       (default 0.1)
    --dz    VAL     |z residual| threshold in cm      (default 5.0)
    --dslope VAL    |slope error| threshold            (default 0.01)
"""

import sys
import re
import argparse
import math

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Filter and format KLong_save_vectors.C DIAG output")
parser.add_argument("--all",    action="store_true",
                    help="print all events, not just those passing all criteria")
parser.add_argument("--log",    metavar="FILE",
                    help="write clean-event summary to FILE")
parser.add_argument("--poca",   type=float, default=1.0,
                    help="PoCA separation threshold (cm)  [default: 1.0]")
parser.add_argument("--dt",     type=float, default=0.1,
                    help="|time diff| threshold (ns)       [default: 0.1]")
parser.add_argument("--dz",     type=float, default=5.0,
                    help="|z residual| threshold (cm)      [default: 5.0]")
parser.add_argument("--dslope", type=float, default=0.01,
                    help="|slope error| threshold           [default: 0.01]")
args = parser.parse_args()

# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------
log_fh = open(args.log, "w") if args.log else None

def out(text=""):
    print(text)
    if log_fh:
        log_fh.write(text + "\n")

PASS = "\033[32mPASS\033[0m"
FAIL = "\033[31mFAIL\033[0m"
WARN = "\033[33mWARN\033[0m"

def status(ok):
    return PASS if ok else FAIL

# ---------------------------------------------------------------------------
# Data structures for one event
# ---------------------------------------------------------------------------
class DiagEvent:
    def __init__(self):
        self.event_num  = None
        # [1]
        self.poca_sep   = None
        self.pt1_z      = None
        self.pt2_z      = None
        # [2]
        self.pip_dt     = None
        self.pim_dt     = None
        self.avg_dt     = None
        self.abs_diff   = None
        # [3]
        self.reco_xyz   = None   # (x, y, z)
        self.true_xyz   = None
        self.delta_xyz  = None
        # [4] pip
        self.pip_ax_fit = None; self.pip_ay_fit = None
        self.pip_ax_exp = None; self.pip_ay_exp = None
        self.pip_dax    = None; self.pip_day    = None
        # [4] pim
        self.pim_ax_fit = None; self.pim_ay_fit = None
        self.pim_ax_exp = None; self.pim_ay_exp = None
        self.pim_dax    = None; self.pim_day    = None
        # [5]
        self.pip_extrap = None   # (x, y)  at true_vz
        self.pim_extrap = None
        self.pip_extrap_delta = None
        self.pim_extrap_delta = None
        self.extrap_tvz = None
        # kinematics
        self.true_p     = None
        self.reco_p     = None
        self.ratio      = None
        # DEBUG tracker counts (from DEBUG block matching same event_num)
        self.t1_tot = None; self.t2_tot = None
        self.t3_tot = None; self.t4_tot = None

class DebugEvent:
    def __init__(self):
        self.event_num = None
        # Per pion: list of (pip_counts, pim_counts) where counts = dict T1..T4
        self.pip_counts = {}
        self.pim_counts = {}

# ---------------------------------------------------------------------------
# Regex patterns
# ---------------------------------------------------------------------------
RE_DIAG_HDR = re.compile(
    r"=== DIAG EVENT\s+(\d+)\s+\(\d+/\d+\)\s+===\s+"
    r"true_p=([\d.]+)\s+GeV/c\s+reco_p=([\d.]+)\s+GeV/c\s+ratio=([\d.]+)")
RE_DIAG_1 = re.compile(
    r"\[1\] PoCA separation:\s+([\d.]+)\s+cm.*point1_z=([-\d.]+).*point2_z=([-\d.]+)")
RE_DIAG_2 = re.compile(
    r"\[2\].*pip=([-\d.]+)\s+pim=([-\d.]+)\s+avg=([-\d.]+)\s+\|diff\|=([\d.]+)")
RE_DIAG_3_RECO = re.compile(
    r"reco:\s+\(([-\d.]+),\s*([-\d.]+),\s*([-\d.]+)\)")
RE_DIAG_3_TRUE = re.compile(
    r"true:\s+\(([-\d.]+),\s*([-\d.]+),\s*([-\d.]+)\)")
RE_DIAG_3_DELTA = re.compile(
    r"delta:\(([-\d.]+),\s*([-\d.]+),\s*([-\d.]+)\)")
RE_DIAG_4_PIP_FIT = re.compile(
    r"pip\s+fit:\s+ax=([-\d.]+)\s+ay=([-\d.]+)")
RE_DIAG_4_PIP_EXP = re.compile(
    r"pip\s+exp:\s+ax=([-\d.]+)\s+ay=([-\d.]+)")
RE_DIAG_4_PIP_D   = re.compile(
    r"pip\s+dax=([-\d.]+)\s+day=([-\d.]+)")
RE_DIAG_4_PIM_FIT = re.compile(
    r"pim\s+fit:\s+ax=([-\d.]+)\s+ay=([-\d.]+)")
RE_DIAG_4_PIM_EXP = re.compile(
    r"pim\s+exp:\s+ax=([-\d.]+)\s+ay=([-\d.]+)")
RE_DIAG_4_PIM_D   = re.compile(
    r"pim\s+dax=([-\d.]+)\s+day=([-\d.]+)")
RE_DIAG_5_PIP = re.compile(
    r"pip:\s+\(([-\d.]+),\s*([-\d.]+)\)\s+true vtx:\s+\(([-\d.]+),\s*([-\d.]+)\)\s+delta:\s+\(([-\d.]+),\s*([-\d.]+)\)")
RE_DIAG_5_PIM = re.compile(
    r"pim:\s+\(([-\d.]+),\s*([-\d.]+)\)\s+true vtx:\s+\(([-\d.]+),\s*([-\d.]+)\)\s+delta:\s+\(([-\d.]+),\s*([-\d.]+)\)")
RE_DIAG_5_TVZ = re.compile(
    r"\[5\].*true_vz=([-\d.]+)\s+cm")

RE_DEBUG_HDR = re.compile(
    r"=== DEBUG EVENT\s+(\d+)\s+\(\d+/\d+\)")
RE_DEBUG_COUNTS = re.compile(
    r"Hit counts:\s+T1=(\d+)\s+T2=(\d+)\s+T3=(\d+)\s+T4=(\d+)")
RE_DEBUG_PION_HDR = re.compile(r"\s+(pi[+-]):")

# ---------------------------------------------------------------------------
# Parse input line by line
# ---------------------------------------------------------------------------
diag_events   = {}   # event_num -> DiagEvent
debug_events  = {}   # event_num -> DebugEvent

current_diag  = None
current_debug = None
current_pion  = None   # "pi+" or "pi-" while inside a debug block

for raw_line in sys.stdin:
    line = raw_line.rstrip("\n")

    # ------------------------------------------------------------------
    # DEBUG block header
    # ------------------------------------------------------------------
    m = RE_DEBUG_HDR.search(line)
    if m:
        current_debug = DebugEvent()
        current_debug.event_num = int(m.group(1))
        current_pion = None
        debug_events[current_debug.event_num] = current_debug
        current_diag = None          # DEBUG and DIAG blocks don't overlap
        continue

    # ------------------------------------------------------------------
    # DIAG block header
    # ------------------------------------------------------------------
    m = RE_DIAG_HDR.search(line)
    if m:
        current_diag = DiagEvent()
        current_diag.event_num = int(m.group(1))
        current_diag.true_p    = float(m.group(2))
        current_diag.reco_p    = float(m.group(3))
        current_diag.ratio     = float(m.group(4))
        diag_events[current_diag.event_num] = current_diag
        current_debug = None
        continue

    # ------------------------------------------------------------------
    # Inside a DEBUG block
    # ------------------------------------------------------------------
    if current_debug is not None:
        m = RE_DEBUG_PION_HDR.match(line)
        if m:
            current_pion = m.group(1)
            continue
        m = RE_DEBUG_COUNTS.search(line)
        if m and current_pion:
            counts = {
                "T1": int(m.group(1)), "T2": int(m.group(2)),
                "T3": int(m.group(3)), "T4": int(m.group(4)),
            }
            if current_pion == "pi+":
                current_debug.pip_counts = counts
            else:
                current_debug.pim_counts = counts
        continue   # still in debug block until we see a non-matching line
        # (blocks end when a new header or blank line arrives, but we
        #  don't need to track the end explicitly)

    # ------------------------------------------------------------------
    # Inside a DIAG block
    # ------------------------------------------------------------------
    if current_diag is not None:
        if RE_DIAG_1.search(line):
            m = RE_DIAG_1.search(line)
            current_diag.poca_sep = float(m.group(1))
            current_diag.pt1_z   = float(m.group(2))
            current_diag.pt2_z   = float(m.group(3))
        elif RE_DIAG_2.search(line):
            m = RE_DIAG_2.search(line)
            current_diag.pip_dt   = float(m.group(1))
            current_diag.pim_dt   = float(m.group(2))
            current_diag.avg_dt   = float(m.group(3))
            current_diag.abs_diff = float(m.group(4))
        elif RE_DIAG_3_RECO.search(line):
            m = RE_DIAG_3_RECO.search(line)
            current_diag.reco_xyz = (float(m.group(1)), float(m.group(2)), float(m.group(3)))
        elif RE_DIAG_3_TRUE.search(line):
            m = RE_DIAG_3_TRUE.search(line)
            current_diag.true_xyz = (float(m.group(1)), float(m.group(2)), float(m.group(3)))
        elif RE_DIAG_3_DELTA.search(line):
            m = RE_DIAG_3_DELTA.search(line)
            current_diag.delta_xyz = (float(m.group(1)), float(m.group(2)), float(m.group(3)))
        elif RE_DIAG_4_PIP_FIT.search(line):
            m = RE_DIAG_4_PIP_FIT.search(line)
            current_diag.pip_ax_fit = float(m.group(1))
            current_diag.pip_ay_fit = float(m.group(2))
        elif RE_DIAG_4_PIP_EXP.search(line):
            m = RE_DIAG_4_PIP_EXP.search(line)
            current_diag.pip_ax_exp = float(m.group(1))
            current_diag.pip_ay_exp = float(m.group(2))
        elif RE_DIAG_4_PIP_D.search(line):
            m = RE_DIAG_4_PIP_D.search(line)
            current_diag.pip_dax = float(m.group(1))
            current_diag.pip_day = float(m.group(2))
        elif RE_DIAG_4_PIM_FIT.search(line):
            m = RE_DIAG_4_PIM_FIT.search(line)
            current_diag.pim_ax_fit = float(m.group(1))
            current_diag.pim_ay_fit = float(m.group(2))
        elif RE_DIAG_4_PIM_EXP.search(line):
            m = RE_DIAG_4_PIM_EXP.search(line)
            current_diag.pim_ax_exp = float(m.group(1))
            current_diag.pim_ay_exp = float(m.group(2))
        elif RE_DIAG_4_PIM_D.search(line):
            m = RE_DIAG_4_PIM_D.search(line)
            current_diag.pim_dax = float(m.group(1))
            current_diag.pim_day = float(m.group(2))
        elif RE_DIAG_5_TVZ.search(line):
            m = RE_DIAG_5_TVZ.search(line)
            current_diag.extrap_tvz = float(m.group(1))
        elif RE_DIAG_5_PIP.search(line):
            m = RE_DIAG_5_PIP.search(line)
            current_diag.pip_extrap       = (float(m.group(1)), float(m.group(2)))
            current_diag.pip_extrap_delta = (float(m.group(5)), float(m.group(6)))
        elif RE_DIAG_5_PIM.search(line):
            m = RE_DIAG_5_PIM.search(line)
            current_diag.pim_extrap       = (float(m.group(1)), float(m.group(2)))
            current_diag.pim_extrap_delta = (float(m.group(5)), float(m.group(6)))

# ---------------------------------------------------------------------------
# Merge tracker coverage from DEBUG blocks into DIAG events
# ---------------------------------------------------------------------------
for evnum, de in diag_events.items():
    db = debug_events.get(evnum)
    if db:
        pip_c = db.pip_counts
        pim_c = db.pim_counts
        de.t1_tot = pip_c.get("T1", 0) + pim_c.get("T1", 0)
        de.t2_tot = pip_c.get("T2", 0) + pim_c.get("T2", 0)
        de.t3_tot = pip_c.get("T3", 0) + pim_c.get("T3", 0)
        de.t4_tot = pip_c.get("T4", 0) + pim_c.get("T4", 0)

# ---------------------------------------------------------------------------
# Apply selection criteria, collect results
# ---------------------------------------------------------------------------
THR_POCA  = args.poca
THR_DT    = args.dt
THR_DZ    = args.dz
THR_DSLP  = args.dslope

def fmt(val, fmt_str):
    return fmt_str % val if val is not None else "n/a"

def apply_criteria(ev):
    """Returns (overall_pass, dict of criterion -> (pass_bool, value_str))."""
    results = {}

    # [1] PoCA separation
    ok1 = (ev.poca_sep is not None) and (ev.poca_sep < THR_POCA)
    results["[1] PoCA sep"] = (ok1, fmt(ev.poca_sep, "%.3f cm"))

    # [2] Time diff
    ok2 = (ev.abs_diff is not None) and (ev.abs_diff < THR_DT)
    results["[2] |dt|"] = (ok2, fmt(ev.abs_diff, "%.4f ns"))

    # [3] |delta z|
    dz = abs(ev.delta_xyz[2]) if ev.delta_xyz else None
    ok3 = (dz is not None) and (dz < THR_DZ)
    results["[3] |dz|"] = (ok3, fmt(dz, "%.3f cm"))

    # [4] Slope errors (all four components)
    slope_vals = [ev.pip_dax, ev.pip_day, ev.pim_dax, ev.pim_day]
    if all(v is not None for v in slope_vals):
        ok4 = all(abs(v) < THR_DSLP for v in slope_vals)
        slope_str = ("pip(%.4f,%.4f)  pim(%.4f,%.4f)" %
                     (ev.pip_dax, ev.pip_day, ev.pim_dax, ev.pim_day))
    else:
        ok4 = False
        slope_str = "n/a"
    results["[4] dslope"] = (ok4, slope_str)

    # [5] Tracker coverage
    if ev.t1_tot is not None:
        ok5 = (ev.t1_tot > 0 and ev.t2_tot > 0 and
               ev.t3_tot > 0 and ev.t4_tot > 0)
        cov_str = "T1=%d T2=%d T3=%d T4=%d" % (
            ev.t1_tot, ev.t2_tot, ev.t3_tot, ev.t4_tot)
    else:
        ok5 = None   # DEBUG block may not be present for this event
        cov_str = "(no DEBUG block – assuming OK)" if ok5 is None else ""
        ok5 = True   # default pass if we have no data
    results["[5] coverage"] = (ok5, cov_str)

    overall = all(ok for ok, _ in results.values())
    return overall, results

# ---------------------------------------------------------------------------
# Print output
# ---------------------------------------------------------------------------
divider       = "=" * 72
thin_divider  = "-" * 72

n_total = len(diag_events)
n_pass  = 0

out()
out(divider)
out("  DIAG EVENT SUMMARY   (thresholds: PoCA<%.2f cm  |dt|<%.3f ns  "
    "|dz|<%.1f cm  |dslope|<%.3f)" % (THR_POCA, THR_DT, THR_DZ, THR_DSLP))
out(divider)

for evnum in sorted(diag_events.keys()):
    ev = diag_events[evnum]
    overall, results = apply_criteria(ev)

    if not (overall or args.all):
        continue

    if overall:
        n_pass += 1

    tag = "CLEAN EVENT" if overall else "event (failed)"
    out()
    out(divider)
    out("  %s  |  Event %-6d  |  true_p=%s GeV/c  reco_p=%s GeV/c  ratio=%s" % (
        tag,
        evnum,
        fmt(ev.true_p, "%.3f"),
        fmt(ev.reco_p, "%.3f"),
        fmt(ev.ratio,  "%.4f"),
    ))
    out(divider)

    # --- Checklist -------------------------------------------------------
    out("  Criteria:")
    for label, (ok, val_str) in results.items():
        ok_str = status(ok)
        out("    %s  %-18s  %s" % (ok_str, label, val_str))

    # --- Vertex ----------------------------------------------------------
    if ev.reco_xyz and ev.true_xyz and ev.delta_xyz:
        out()
        out("  Vertex (cm):")
        out("    reco  x=%8.3f  y=%8.3f  z=%8.3f" % ev.reco_xyz)
        out("    true  x=%8.3f  y=%8.3f  z=%8.3f" % ev.true_xyz)
        out("    delta x=%+8.3f  y=%+8.3f  z=%+8.3f" % ev.delta_xyz)

    # --- Track slopes ----------------------------------------------------
    if ev.pip_ax_fit is not None:
        out()
        out("  Track slopes (dx/dz, dy/dz):")
        out("    pi+  fit: ax=%+.5f  ay=%+.5f" % (ev.pip_ax_fit, ev.pip_ay_fit))
        out("    pi+  exp: ax=%+.5f  ay=%+.5f" % (ev.pip_ax_exp, ev.pip_ay_exp))
        out("    pi+  err: dax=%+.5f  day=%+.5f" % (ev.pip_dax, ev.pip_day))
        out("    pi-  fit: ax=%+.5f  ay=%+.5f" % (ev.pim_ax_fit, ev.pim_ay_fit))
        out("    pi-  exp: ax=%+.5f  ay=%+.5f" % (ev.pim_ax_exp, ev.pim_ay_exp))
        out("    pi-  err: dax=%+.5f  day=%+.5f" % (ev.pim_dax, ev.pim_day))

    # --- Extrapolation at true_vz ----------------------------------------
    if ev.pip_extrap and ev.pim_extrap:
        out()
        out("  Track extrapolated to true_vz=%.3f cm:" % (ev.extrap_tvz or 0.))
        out("    pi+  at_vz=(%8.3f, %8.3f)  delta=(%+.4f, %+.4f)" % (
            ev.pip_extrap[0], ev.pip_extrap[1],
            ev.pip_extrap_delta[0], ev.pip_extrap_delta[1]))
        out("    pi-  at_vz=(%8.3f, %8.3f)  delta=(%+.4f, %+.4f)" % (
            ev.pim_extrap[0], ev.pim_extrap[1],
            ev.pim_extrap_delta[0], ev.pim_extrap_delta[1]))

    # --- PoCA detail & timing --------------------------------------------
    out()
    out("  PoCA detail:")
    out("    separation=%.4f cm  |  pt1_z=%.3f cm  pt2_z=%.3f cm" % (
        ev.poca_sep or 0., ev.pt1_z or 0., ev.pt2_z or 0.))
    if ev.pip_dt is not None:
        out("  Decay times (ns):")
        out("    pi+=%+.4f  pi-=%+.4f  avg=%+.4f  |diff|=%.4f" % (
            ev.pip_dt, ev.pim_dt, ev.avg_dt, ev.abs_diff))

# ---------------------------------------------------------------------------
# Final tally
# ---------------------------------------------------------------------------
out()
out(divider)
out("  TOTAL DIAG EVENTS : %d" % n_total)
out("  PASSED ALL CRITERIA: %d" % n_pass)
out("  FAILED             : %d" % (n_total - n_pass))
out(divider)
out()

if log_fh:
    log_fh.close()
    print("(summary also written to %s)" % args.log)
