"""
textlog.py -- tee stdout into a text file in output/.

Each runnable script wraps its __main__ body in `tee_stdout("<name>.txt")` so
the printed summary (rates, angles, sanity numbers) is persisted next to the
figures instead of living only in the terminal.  No physics here.
"""

import sys
from contextlib import contextmanager
from pathlib import Path

OUTDIR = Path(__file__).resolve().parent.parent / "output"


class _Tee:
    """Write-through to several streams (terminal + log file)."""

    def __init__(self, *streams):
        self._streams = streams

    def write(self, text):
        for s in self._streams:
            s.write(text)

    def flush(self):
        for s in self._streams:
            s.flush()


@contextmanager
def tee_stdout(filename):
    """Duplicate everything printed inside the block into output/<filename>."""
    path = OUTDIR / filename
    with open(path, "w", encoding="utf-8") as f:
        orig = sys.stdout
        sys.stdout = _Tee(orig, f)
        try:
            yield path
        finally:
            sys.stdout = orig
