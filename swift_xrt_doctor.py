#!/usr/bin/env python3
"""swift_xrt_doctor.py -- environment diagnostic for the Swift XRT pipeline.

Prints a green/red checklist of everything the pipeline assumes about the
environment (PATH, HEASoft, CALDB, Python, packages, Sherpa, disk) so a user
can verify readiness in one shot instead of watching the pipeline fail midway.

Exit status: 0 only if no check FAILs. Stdlib only.

Usage:
    swift_xrt_doctor.py            # full human-readable report (color on tty)
    swift_xrt_doctor.py --quiet    # only print failing checks
    swift_xrt_doctor.py --no-color # plain text, no ANSI
"""

import argparse
import importlib.metadata as _md
import importlib.util as _ilu
import os
import re
import shutil
import sys

# Directory this script lives in; used to confirm PATH points at the pipeline.
PIPELINE_DIR = os.path.dirname(os.path.realpath(__file__))

# HEASoft FTOOLS the pipeline shells out to (see extract/fit scripts).
HEASOFT_TOOLS = ["xrtpipeline", "xrtmkarf", "grppha", "xselect", "ftlist"]

# Required Python packages: (import_name, distribution_name).
REQUIRED_PKGS = [
    ("astropy", "astropy"),
    ("numpy", "numpy"),
    ("scipy", "scipy"),
    ("matplotlib", "matplotlib"),
    ("requests", "requests"),
]

OK, WARN, FAIL = "OK", "WARN", "FAIL"

_COLORS = {OK: "\033[32m", WARN: "\033[33m", FAIL: "\033[31m"}
_RESET = "\033[0m"

_results = []  # list of (status, message, [detail lines])
_use_color = False
_quiet = False


def emit(status, message, details=None):
    """Record a check result and print it (unless suppressed by --quiet)."""
    _results.append((status, message, details or []))
    if _quiet and status != FAIL:
        return
    tag = "[%s]" % status
    if _use_color:
        tag = _COLORS[status] + tag + _RESET
    print("%s %s" % (tag, message))
    for line in (details or []):
        print("       %s" % line)


def _pkg_version(dist):
    try:
        return _md.version(dist)
    except Exception:
        return "?"


# --- individual checks ------------------------------------------------------

def check_path():
    resolved = shutil.which("swift_xrt_summary.py")
    if resolved and os.path.realpath(resolved).startswith(PIPELINE_DIR):
        emit(OK, "Pipeline on PATH: swift_xrt_summary.py -> %s" % resolved)
    elif resolved:
        emit(WARN, "swift_xrt_summary.py resolves to %s, not %s "
             "(run setup_swiftxrt)" % (resolved, PIPELINE_DIR))
    else:
        emit(FAIL, "Pipeline not on PATH (no swift_xrt_summary.py). "
             "Run setup_swiftxrt.")


def check_heasoft():
    headas = os.environ.get("HEADAS")
    if not headas:
        emit(FAIL, "HEASoft not loaded: $HEADAS unset. Run 'heainit'.")
        return
    missing, details = [], []
    for tool in HEASOFT_TOOLS:
        path = shutil.which(tool)
        if path:
            details.append("%-12s %s" % (tool, path))
        else:
            missing.append(tool)
            details.append("%-12s NOT FOUND" % tool)
    if missing:
        emit(FAIL, "HEASoft tools missing from PATH: %s (run 'heainit')"
             % ", ".join(missing), details)
    else:
        emit(OK, "HEASoft loaded ($HEADAS set, all FTOOLS on PATH)", details)


def check_heasoft_version():
    # $HEADAS can be clobbered (e.g. CIAO repoints it at its own bundled
    # spectral tools), so also try the resolved xrtpipeline path.
    candidates = [os.environ.get("HEADAS", "")]
    tool = shutil.which("xrtpipeline")
    if tool:
        candidates.append(os.path.realpath(tool))
    for src in candidates:
        m = re.search(r"heasoft[-_]?(\d+)\.(\d+)(?:\.(\d+))?", src, re.I)
        if not m:
            continue
        major, minor = int(m.group(1)), int(m.group(2))
        ver = ".".join(p for p in m.groups() if p)
        if (major, minor) < (6, 30):
            emit(WARN, "HEASoft %s looks old (<6.30); consider upgrading" % ver)
        else:
            emit(OK, "HEASoft version %s" % ver)
        return
    emit(WARN, "Could not parse HEASoft version (HEADAS=%s)"
         % (os.environ.get("HEADAS") or "unset"))


def _looks_like_ciao(path):
    return "ciao" in path.lower()


def check_caldb():
    caldb = os.environ.get("CALDB")
    cfg = os.environ.get("CALDBCONFIG")
    if not caldb:
        emit(FAIL, "$CALDB unset. Run 'caldbinit' (HEASoft CALDB).")
        return
    if not os.path.isdir(caldb):
        emit(FAIL, "$CALDB=%s but directory does not exist" % caldb)
        return
    if _looks_like_ciao(caldb):
        emit(WARN, "$CALDB points inside a CIAO install (%s)" % caldb,
             ["This is the Chandra CALDB, not HEASoft's Swift CALDB.",
              "Pass --caldb /opt/CALDB to parallel_fit.py / "
              "swift_xrt_fit_spectra.py."])
        return
    details = ["$CALDB=%s" % caldb]
    if not cfg:
        emit(WARN, "$CALDB ok but $CALDBCONFIG unset (run 'caldbinit')",
             details)
    elif not os.path.isfile(cfg):
        emit(WARN, "$CALDBCONFIG=%s does not exist" % cfg, details)
    else:
        details.append("$CALDBCONFIG=%s" % cfg)
        emit(OK, "CALDB configured", details)


def check_swift_caldb():
    caldb = os.environ.get("CALDB", "")
    rel = "data/swift/xrt/cpf/rmf"
    if caldb and os.path.isdir(os.path.join(caldb, rel)):
        emit(OK, "Swift XRT response files present under $CALDB")
        return
    # Fall back to the conventional HEASoft CALDB on this host.
    if os.path.isdir(os.path.join("/opt/CALDB", rel)):
        emit(WARN, "Swift XRT RMFs found at /opt/CALDB but not under $CALDB",
             ["Point $CALDB at /opt/CALDB or pass --caldb /opt/CALDB."])
    else:
        emit(FAIL, "Swift XRT response files not found (%s under $CALDB)" % rel)


def check_python():
    v = sys.version_info
    msg = "Python %d.%d.%d (%s)" % (v.major, v.minor, v.micro, sys.executable)
    if (v.major, v.minor) < (3, 9):
        emit(FAIL, msg + " -- 3.9+ required")
    else:
        emit(OK, msg)


def check_required_pkgs():
    for imp, dist in REQUIRED_PKGS:
        if _ilu.find_spec(imp) is not None:
            emit(OK, "%-12s %s" % (imp, _pkg_version(dist)))
        else:
            emit(FAIL, "%-12s NOT IMPORTABLE (pip install %s)" % (imp, dist))


def check_optional_pkgs():
    if _ilu.find_spec("astroquery") is not None:
        emit(OK, "astroquery   %s (optional)" % _pkg_version("astroquery"))
    else:
        emit(WARN, "astroquery not installed (optional; download script falls "
             "back to SIMBAD/NED/Sesame)")


def check_sherpa():
    if _ilu.find_spec("sherpa") is not None:
        emit(OK, "sherpa       %s" % _pkg_version("sherpa"))
    else:
        emit(FAIL, "sherpa not importable -- needed for fitting. This usually "
             "means you are outside a CIAO environment (source ciao.bash).")


def check_disk():
    try:
        usage = shutil.disk_usage("/opt")
    except OSError as exc:
        emit(WARN, "Could not stat /opt: %s" % exc)
        return
    free_gb = usage.free / (1024 ** 3)
    msg = "Disk free at /opt: %.1f GB" % free_gb
    if free_gb < 10:
        emit(WARN, msg + " (<10 GB; reductions may run out of space)")
    else:
        emit(OK, msg)


CHECKS = [
    check_path,
    check_heasoft,
    check_heasoft_version,
    check_caldb,
    check_swift_caldb,
    check_python,
    check_required_pkgs,
    check_optional_pkgs,
    check_sherpa,
    check_disk,
]


def main(argv=None):
    global _use_color, _quiet
    parser = argparse.ArgumentParser(
        description="Verify the Swift XRT pipeline environment.")
    parser.add_argument("--quiet", action="store_true",
                        help="only print failing checks")
    parser.add_argument("--no-color", action="store_true",
                        help="disable ANSI color output")
    args = parser.parse_args(argv)

    _quiet = args.quiet
    _use_color = sys.stdout.isatty() and not args.no_color \
        and os.environ.get("NO_COLOR") is None

    for check in CHECKS:
        try:
            check()
        except Exception as exc:  # a check must never crash the doctor
            emit(FAIL, "%s crashed: %s" % (check.__name__, exc))

    n_ok = sum(1 for s, _, _ in _results if s == OK)
    n_warn = sum(1 for s, _, _ in _results if s == WARN)
    n_fail = sum(1 for s, _, _ in _results if s == FAIL)

    summary = "%d checks: %d ok, %d warn, %d fail" % (
        len(_results), n_ok, n_warn, n_fail)
    if _use_color:
        color = _COLORS[FAIL] if n_fail else (_COLORS[WARN] if n_warn
                                              else _COLORS[OK])
        summary = color + summary + _RESET
    print("\n" + summary)
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
