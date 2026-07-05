"""
detector_scan_common.py -- shared machinery for the fixed-detector /
rotating-crystal scans (plot_detector_scan.py, plot_detector_scan_bandpass.py).

Setup common to both diagrams: a detector sits at a FIXED lab (external) angle
theta_s (3 deg by default).  The crystal is rotated (varying the internal
phase-matching angle theta_m) and the photon count is recorded vs rotation.

Detector model: count into the full azimuthal ring at theta_s within an
angular acceptance dth_det (default 2 mrad ~ the paper's angular resolution).
A real point detector subtending azimuth dphi sees a fixed fraction dphi/2pi
of this ring -- that only rescales the y-axis, leaving the shape and the
orders-of-magnitude analysis unchanged.

All physics comes from spdc_physics; this module only orchestrates the scan,
locates the peak/threshold crossings, and draws the standard panel.
"""

import numpy as np
import spdc_physics as phys

# Signal wavelengths whose idler exceeds the BBO transparency edge (~3500 nm;
# also the limit of the refractive-index data) are unphysical: such pairs
# cannot actually be generated (e.g. lam_s = 433 nm -> lam_i = 6370 nm).
# Without this cutoff the short-lam_s tail of the broadband scan is propped up
# by them.
LAM_I_MAX = 3500e-9


def spectra_vs_thetam(theta_det, tms_deg, lam, dlam, dth_det, Nxi=400):
    """Per-wavelength detector count spectrum for each crystal angle.

    theta_det : fixed external detector angle, rad.
    tms_deg   : array of internal phase-matching angles theta_m, DEGREES.
    lam, dlam : signal wavelength grid (m) and its spacing.
    dth_det   : detector polar acceptance, rad.

    Returns array of shape (len(tms_deg), len(lam)) in counts/s per lam-bin.
    Every panel quantity (broadband sum, filtered sum, spectral density) is
    derived from this one matrix, so the panels are guaranteed consistent."""
    return np.array([phys.Ns(lam, theta_det, np.radians(t), dlam, dth_det,
                             Nxi=Nxi, lam_i_max=LAM_I_MAX)
                     for t in tms_deg])


def rotation_axes(tms_deg, theta_m_cut):
    """Map the theta_m grid to the external crystal-rotation angle alpha and
    build the interpolators used for the secondary (top) axis.

    Returns (alphas_deg, alpha_to_thetam, thetam_to_alpha)."""
    alphas = phys.crystal_rotation(tms_deg, theta_m_cut)
    a2t = lambda a: np.interp(a, alphas, tms_deg)
    t2a = lambda t: np.interp(t, tms_deg, alphas)
    return alphas, a2t, t2a


def two_order_crossings(alphas, counts, ipk):
    """Interpolated alpha where counts fall to peak/100 (two orders of
    magnitude) on each side of the peak index ipk.  Interpolation is done in
    log10(counts) since the data spans decades.  If a side never crosses the
    threshold inside the scan range, the end of the range is returned.

    Returns (alpha_low, alpha_high, threshold)."""
    thresh = counts[ipk] / 100.0

    def cross(direction):
        i = ipk
        while 0 <= i < len(counts) and counts[i] >= thresh:
            i += direction
        if not (0 <= i < len(counts)):          # ran off the scan range
            return alphas[i - direction]
        j = i - direction                       # last point still above thresh
        y0, y1 = np.log10(max(counts[j], 1e-300)), np.log10(max(counts[i], 1e-300))
        f = (np.log10(thresh) - y0) / (y1 - y0)
        return alphas[j] + f * (alphas[i] - alphas[j])

    return cross(-1), cross(+1), thresh


def summarize(header, alphas, tms_deg, data, unit="/s"):
    """Print the peak position and the two-orders-down rotation range for one
    scan curve.  Returns (ipk, thresh) for use by draw_panel."""
    ipk = int(np.argmax(data))
    lo, hi, thresh = two_order_crossings(alphas, data, ipk)
    print(f"=== {header} ===")
    print(f"  peak at alpha={alphas[ipk]:+.2f} deg (thetam={tms_deg[ipk]:.2f}), "
          f"value={data[ipk]:.2e}{unit}")
    print(f"  2 orders down at alpha = {lo:+.2f} and {hi:+.2f} deg "
          f"(rotate {lo - alphas[ipk]:+.2f} / {hi - alphas[ipk]:+.2f} from peak)")
    return ipk, thresh


def draw_panel(ax, alphas, data, ipk, thresh, ylabel, title, subtitle, a2t, t2a):
    """Standard scan panel: semilog count vs alpha, dashed line at the peak,
    dotted line at the 2-orders-down threshold, secondary theta_m axis on top."""
    ax.semilogy(alphas, data, color="C0", lw=1.6)
    ax.axhline(thresh, color="0.5", ls=":", lw=1)
    ax.axvline(alphas[ipk], color="C3", ls="--", lw=1)
    ax.set_ylabel(ylabel)
    ax.grid(True, which="both", alpha=0.25)
    ax.set_xlabel(r"crystal rotation angle $\alpha$ (deg, 0 = normal incidence)")
    ax.set_title(title, pad=26)
    ax.text(0.5, 0.02, subtitle, transform=ax.transAxes, ha="center",
            va="bottom", fontsize=7.5, color="0.35")
    ax.set_ylim(data.max() / 3e3, data.max() * 2)
    sec = ax.secondary_xaxis("top", functions=(a2t, t2a))
    sec.set_xlabel(r"internal phase-matching angle $\theta_m$ (deg)")
