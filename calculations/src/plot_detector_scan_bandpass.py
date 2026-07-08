"""
plot_detector_scan_bandpass.py -- fixed-detector / rotating-crystal scan,
three wavelength-selection variants.

Same geometry as plot_detector_scan.py (detector fixed at external
theta_s = 3 deg, crystal rotated), but with reworked wavelength selection:

  1) BROADBAND     -- integrate over ALL signal wavelengths landing on the
                      detector (idler < 3500 nm physical cutoff).
  2) BANDPASS      -- weight the spectrum by a GAUSSIAN filter centred at
                      810 nm with FWHM = 10 nm:
                          T(lam) = exp(-4 ln2 (lam - 810nm)^2 / FWHM^2),
                      peak transmission T = 1.
  3) MONOCHROMATIC -- only "exactly" 810 nm.  An infinitely narrow filter
                      passes zero total flux, so the meaningful quantity is
                      the SPECTRAL DENSITY dN/dlam at 810 nm, in counts/s/nm.

Physics from spdc_physics.py; scan machinery from detector_scan_common.py.
Writes detector_scan_bandpass.png.
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

import detector_scan_common as com

OUTDIR = Path(__file__).resolve().parent.parent / "output"


def gaussian_bandpass(lam, lam0, fwhm):
    """Gaussian filter transmission T(lam), normalised to peak T = 1 at lam0.
    exp(-4 ln2 x^2 / fwhm^2) has T = 0.5 at x = +/- fwhm/2, i.e. exactly the
    requested full width at half maximum."""
    return np.exp(-4.0 * np.log(2.0) * (lam - lam0)**2 / fwhm**2)


def run(theta_det_deg=3.0, dth_det=2e-3, theta_m_cut=None,
        tm_lo=24.0, tm_hi=30.5, dtm=0.01,
        lam_lo=430e-9, lam_hi=1600e-9, dlam=0.5e-9,
        lam0=810e-9, fwhm=10e-9):
    th_det = np.radians(theta_det_deg)
    lam    = np.arange(lam_lo, lam_hi + 1e-12, dlam)
    tms    = np.arange(tm_lo, tm_hi + 1e-12, dtm)      # theta_m grid, deg

    # Gaussian bandpass (panel 2) and the index of the 810 nm bin (panel 3).
    T       = gaussian_bandpass(lam, lam0, fwhm)
    i0      = int(np.argmin(np.abs(lam - lam0)))
    dlam_nm = dlam * 1e9

    # Full per-wavelength spectrum for every theta_m, computed ONCE; all three
    # panels are derived from this matrix so they are guaranteed consistent.
    spectra = com.spectra_vs_thetam(th_det, tms, lam, dlam, dth_det)

    broad = spectra.sum(axis=1)                 # all wavelengths,   counts/s
    band  = (spectra * T[None, :]).sum(axis=1)  # Gaussian bandpass, counts/s
    mono  = spectra[:, i0] / dlam_nm            # dN/dlam at 810 nm, counts/s/nm

    # Cut the crystal at the optimal phase-matching angle so the peak sits at
    # alpha = 0.  Default: centre on the broadband peak.
    if theta_m_cut is None:
        theta_m_cut = tms[int(np.argmax(broad))]
    alphas, a2t, t2a = com.rotation_axes(tms, theta_m_cut)

    fwhm_nm, l0_nm = fwhm * 1e9, lam0 * 1e9
    print(f"Crystal cut at theta_m_cut = {theta_m_cut:.2f} deg "
          f"(= optimal, so alpha=0 is the peak)")
    ib,  bth = com.summarize("BROADBAND (all wavelengths, idler<3500nm)",
                             alphas, tms, broad)
    inr, nth = com.summarize(f"BANDPASS (Gaussian, {l0_nm:.0f} nm centre, "
                             f"{fwhm_nm:.0f} nm FWHM)", alphas, tms, band)
    imo, mth = com.summarize(f"MONOCHROMATIC (exactly {l0_nm:.0f} nm, "
                             f"spectral density)", alphas, tms, mono, unit="/s/nm")

    # ---- plot ----
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8.5, 13.5))

    com.draw_panel(ax1, alphas, broad, ib, bth, "detector count (1/s)",
                   "All wavelengths (broadband)",
                   "always some non-degenerate pair phase-matched at 3deg -> "
                   "plateau + sharp caustic edge", a2t, t2a)
    com.draw_panel(ax2, alphas, band, inr, nth, "detector count (1/s)",
                   f"Gaussian bandpass: {l0_nm:.0f} nm centre, {fwhm_nm:.0f} nm FWHM",
                   "near-degenerate ring weighted by the filter T(lam) -> "
                   "sharp, ~symmetric peak", a2t, t2a)
    com.draw_panel(ax3, alphas, mono, imo, mth, "spectral count (1/s/nm)",
                   f"Monochromatic: exactly {l0_nm:.0f} nm",
                   "zero-bandwidth limit (spectral density) -> the cleanest, "
                   "sharpest peak", a2t, t2a)

    fig.suptitle(f"SPDC count at a fixed {theta_det_deg:.0f}$^\\circ$ detector "
                 f"vs. crystal rotation\n(crystal cut at the optimal "
                 f"$\\theta_m$ = {theta_m_cut:.2f}$^\\circ$, so $\\alpha$=0 is the peak)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.975])
    fig.savefig(OUTDIR / "detector_scan_bandpass.png", dpi=140)
    print(f"saved {OUTDIR / 'detector_scan_bandpass.png'}")


if __name__ == "__main__":
    run()
