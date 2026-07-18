"""
plot_detector_scan.py -- fixed-detector / rotating-crystal scan, two panels.

A detector sits at a FIXED lab (external) angle theta_s = 3 deg; the crystal
is rotated and the wavelength-integrated count is recorded vs the external
rotation angle alpha:

    N_det(theta_m) = SUM over lam_s of  Ns(lam_s, theta_s = 3 deg, theta_m)

Two panels:
  1) BROADBAND  -- all signal wavelengths (with the idler < 3500 nm physical
                   cutoff, see detector_scan_common.LAM_I_MAX).
  2) NARROWBAND -- hard 805-815 nm top-hat band (e.g. a bandpass filter);
                   for the Gaussian-filter / monochromatic versions see
                   plot_detector_scan_bandpass.py.

Key physical result: the broadband curve is strongly ASYMMETRIC (dropping two
orders of magnitude takes only ~ +0.8 deg on the high-alpha side -- the sharp
degenerate-caustic edge -- but ~ -7 deg on the low side, where the tuning-curve
arms keep crossing 3 deg).  The narrowband curve is ~symmetric (+/- 0.3 deg):
that is the "intuitive"/filtered case (paper Fig. 3).

Physics from spdc_physics.py; scan machinery from detector_scan_common.py.
Writes detector_scan.png.
"""

import sys
from pathlib import Path

# archived script: the shared modules live one level up, in src/
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib.pyplot as plt

import detector_scan_common as com

OUTDIR = Path(__file__).resolve().parent.parent.parent / "output" / "archive"


def run(theta_det_deg=3.0, dth_det=2e-3, theta_m_cut=None,
        tm_lo=24.0, tm_hi=30.5, dtm=0.01,
        lam_lo=430e-9, lam_hi=1600e-9, dlam=0.5e-9,
        narrow_band=(805e-9, 815e-9)):
    th_det = np.radians(theta_det_deg)
    lam    = np.arange(lam_lo, lam_hi + 1e-12, dlam)
    tms    = np.arange(tm_lo, tm_hi + 1e-12, dtm)      # theta_m grid, deg

    # Full per-wavelength spectrum for every theta_m (one matrix, both panels
    # derived from it).  Counts depend only on theta_m -- the alpha axis is
    # just a relabelling set by the cut angle.
    spectra = com.spectra_vs_thetam(th_det, tms, lam, dlam, dth_det)

    in_band = (lam >= narrow_band[0]) & (lam <= narrow_band[1])
    broad   = spectra.sum(axis=1)                      # counts/s
    narrow  = spectra[:, in_band].sum(axis=1)          # counts/s

    # Cut the crystal at the optimal phase-matching angle so the peak sits at
    # alpha = 0 (normal incidence = optimal).  Default: centre on the
    # broadband peak.
    if theta_m_cut is None:
        theta_m_cut = tms[int(np.argmax(broad))]
    alphas, a2t, t2a = com.rotation_axes(tms, theta_m_cut)

    nb_lo, nb_hi = narrow_band[0] * 1e9, narrow_band[1] * 1e9
    print(f"Crystal cut at theta_m_cut = {theta_m_cut:.2f} deg "
          f"(= optimal, so alpha=0 is the peak)")
    ib,  bth = com.summarize("BROADBAND (all wavelengths, idler<3500nm)",
                             alphas, tms, broad)
    inr, nth = com.summarize(f"NARROWBAND ({nb_lo:.0f}-{nb_hi:.0f} nm, "
                             f"e.g. with a bandpass filter)",
                             alphas, tms, narrow)

    # ---- plot ----
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.5, 9.5))

    com.draw_panel(ax1, alphas, broad, ib, bth, "detector count (1/s)",
                   "All wavelengths (broadband)",
                   "always some non-degenerate pair phase-matched at 3deg -> "
                   "plateau + sharp caustic edge", a2t, t2a)
    com.draw_panel(ax2, alphas, narrow, inr, nth, "detector count (1/s)",
                   f"Bandpass {nb_lo:.0f}-{nb_hi:.0f} nm (near-degenerate only)",
                   "the bright 810 nm ring sweeps through 3deg -> sharp, "
                   "~symmetric peak (the 'intuitive' picture)", a2t, t2a)

    fig.suptitle(f"SPDC count at a fixed {theta_det_deg:.0f}$^\\circ$ detector "
                 f"vs. crystal rotation\n(crystal cut at the optimal "
                 f"$\\theta_m$ = {theta_m_cut:.2f}$^\\circ$, so $\\alpha$=0 is the peak)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUTDIR / "detector_scan.png", dpi=140)
    print(f"saved {OUTDIR / 'detector_scan.png'}")


if __name__ == "__main__":
    from textlog import tee_stdout
    with tee_stdout("archive/detector_scan.txt"):
        print(f"indices: {com.phys.sellmeier_label()}")
        run()
