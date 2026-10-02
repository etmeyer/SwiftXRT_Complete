"""
swift_xrt_env.py

Shared environment checks for the Swift XRT pipeline scripts.

HEASoft and CIAO cannot share one shell for this pipeline:

  * Once CIAO is set up (`ciao`, i.e. sourcing ciao.sh), `python3`
    resolves to CIAO's wrapper script. That wrapper re-runs
    ciao_setup.sh, which sets HEADAS to CIAO's bundled spectral
    directory inside every Python process it starts. Every HEASoft
    tool a pipeline script then launches (xrtpipeline, xselect,
    xrtmkarf, grppha) inherits that HEADAS and fails, e.g.
        Failed to open .../ciao-4.16/spectral/bin//xselect.mdb
    CIAO also points $CALDB at its Chandra-only calibration tree.

  * The fit step needs Sherpa, which ships only with CIAO's Python.

So the pipeline runs in two terminals:

    HEASoft terminal (no CIAO): download, xrtpipeline, survey,
                                inspection, extraction (docs steps 2-7)
    CIAO terminal:              fitting, with --caldb   (docs step 8)

The pipeline scripts import this module from their own directory.
"""

import os
import shutil
import sys

# Where the Swift XRT response matrices live inside a HEASoft CALDB.
SWIFT_RMF_SUBDIR = os.path.join('data', 'swift', 'xrt', 'cpf', 'rmf')

# Conventional HEASoft CALDB location, suggested in error messages.
DEFAULT_HEASOFT_CALDB = '/opt/CALDB'

HEASOFT_SHELL_HINT = (
    "Run this step from a terminal with HEASoft set up and CIAO NOT set up:\n"
    "    setup_swiftxrt; heainit\n"
    "(elsewhere: source headas-init.sh and caldbinit.sh, but not ciao.sh).\n"
    "Running heainit again in a terminal where `ciao` was already run does\n"
    "not undo CIAO -- open a new terminal.")

CIAO_SHELL_HINT = (
    "Run the fit step from a terminal with CIAO set up:\n"
    "    setup_swiftxrt; ciao\n"
    "and pass the HEASoft CALDB, e.g. --caldb %s" % DEFAULT_HEASOFT_CALDB)


def _is_under(path, root):
    """True if path is root or lies inside it (symlinks resolved)."""
    if not path or not root:
        return False
    path = os.path.realpath(path)
    root = os.path.realpath(root)
    return path == root or path.startswith(root + os.sep)


def ciao_install():
    """CIAO install root if CIAO is set up in this environment, else ''."""
    return os.environ.get('ASCDS_INSTALL', '')


def headas_is_ciao():
    """True if $HEADAS, as this process sees it, points into CIAO."""
    return _is_under(os.environ.get('HEADAS', ''), ciao_install())


def caldb_has_swift(caldb):
    """
    True/False for a local CALDB directory; None when it cannot be
    checked (unset, or a remote CALDB URL).
    """
    if not caldb or '://' in caldb:
        return None
    return os.path.isdir(os.path.join(caldb, SWIFT_RMF_SUBDIR))


def _caldb_problem(caldb):
    """Describe why caldb can't serve Swift XRT, or None if it can."""
    if caldb_has_swift(caldb) is not False:
        return None
    if _is_under(caldb, ciao_install()):
        return ("CALDB=%s is CIAO's Chandra calibration tree; it has no "
                "Swift XRT files." % caldb)
    return "No Swift XRT calibration files under CALDB=%s" % caldb


def require_heasoft_shell(tools):
    """
    Exit with an actionable message unless HEASoft (and not CIAO) is
    set up for this process. Returns (headas, caldb) on success.
    """
    headas = os.environ.get('HEADAS', '')
    caldb = os.environ.get('CALDB', '')
    problems = []

    if headas_is_ciao():
        problems.append(
            "CIAO is set up in this terminal, and CIAO's python has "
            "replaced $HEADAS with %s, so the HEASoft tools this step "
            "runs cannot find their files." % headas)
    elif not headas:
        problems.append("HEADAS is not set (HEASoft is not initialized).")

    if not caldb:
        problems.append("CALDB is not set.")
    else:
        caldb_problem = _caldb_problem(caldb)
        if caldb_problem:
            problems.append(caldb_problem)

    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        problems.append("Not on PATH: %s" % ', '.join(missing))

    if problems:
        print("ERROR: this step needs a HEASoft environment without CIAO.",
              file=sys.stderr)
        for problem in problems:
            print("  - " + problem, file=sys.stderr)
        print(HEASOFT_SHELL_HINT, file=sys.stderr)
        sys.exit(1)

    return headas, caldb


def sherpa_missing_message():
    """Explain how to get Sherpa when it is not importable."""
    return ("ERROR: sherpa is not importable from this Python (%s).\n"
            "Sherpa ships with CIAO.\n%s\n"
            "(Alternatively, pip install sherpa into this Python.)"
            % (sys.executable, CIAO_SHELL_HINT))


def require_fit_caldb(caldb_override):
    """
    Exit with an actionable message unless the CALDB the fit will use
    (--caldb, else $CALDB) contains the Swift XRT response files.
    Returns the CALDB path used.
    """
    caldb = caldb_override or os.environ.get('CALDB', '')
    source = '--caldb' if caldb_override else '$CALDB'
    if not caldb:
        print("ERROR: no CALDB to resolve Swift XRT responses from "
              "(neither --caldb nor $CALDB is set).", file=sys.stderr)
        print("  Pass --caldb %s" % DEFAULT_HEASOFT_CALDB, file=sys.stderr)
        sys.exit(1)
    problem = _caldb_problem(caldb)
    if problem:
        print("ERROR: %s (from %s)" % (problem, source), file=sys.stderr)
        print("  Without the Swift RMFs every fit fails. Pass the HEASoft "
              "CALDB, e.g. --caldb %s" % DEFAULT_HEASOFT_CALDB,
              file=sys.stderr)
        sys.exit(1)
    return caldb
