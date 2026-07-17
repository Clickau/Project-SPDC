"""
detector_counts.py -- photons/s on a real, finite-aperture detector at a fixed
lab angle, for THIS setup (not the paper's), plus the crystal-rotation scan.

Every experiment knob lives in the PARAMETERS block below -- including the
detector lab angle.  spdc_physics keeps the paper values as its defaults; this
script overrides them per call, so the existing scripts are untouched.

Setup (requested 2026-07-05):
  pump 405 nm, 40 mW, collimated, 2 mm beam diameter (W = 1 mm)
  BBO, L = 0.1 mm, deff = 1.75 pm/V, cut at the IDEAL angle for the detector
    position (degenerate cone exactly through the detector at normal incidence)
  detector at theta_s = 3 deg (lab/external), 625 mm from the crystal,
    6 mm circular opening
  Gaussian bandpass, 10 nm FWHM, centred on the degenerate wavelength
  optical-component efficiencies (EFFICIENCIES block): pump chain
    laser -> F01 mirror -> HWP -> F01 mirror -> crystal face scales the pump
    power; collection chain crystal face -> P01 mirror -> bandpass peak ->
    lens/fiber coupling -> fiber -> APD scales the detected rate.
    Datasheet values filled in 2026-07-05.

Detector model
--------------
phys.Ns gives counts/s into the FULL 2*pi azimuthal ring per (dlam x dtheta)
bin.  A circular aperture (radius r, distance D) centred on the ring at
R0 = D tan(theta_s) intersects the cone of polar angle theta (ring radius
rho = D tan(theta)) over the azimuthal arc

    dphi(theta) = 2 arccos[(rho^2 + R0^2 - r^2) / (2 rho R0)] ,  |rho - R0| < r,

so the rate on the detector is

    N = INT dtheta [dN/dtheta](theta) * dphi(theta) / (2*pi) ,

with [dN/dtheta] the filter-weighted, wavelength-integrated ring density.
This is the count on ONE detector (the twin at -theta_s sees the same rate).

Outputs: prints the detected rate at normal incidence (ideal cut) and writes
detector_scan_mysetup.png -- the rate vs crystal rotation angle alpha.

Run:  python detector_counts.py
"""

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import spdc_physics as phys
import detector_scan_common as com

OUTDIR = Path(__file__).resolve().parent.parent / "output"

# ----------------------------------------------------------------------------
# PARAMETERS (this setup; paper values remain the defaults inside spdc_physics)
# ----------------------------------------------------------------------------
LAM_P  = 405e-9        # pump wavelength, m
PP     = 40e-3         # pump power, W
W_PUMP = 1e-3          # pump 1/e^2 intensity RADIUS, m (collimated, 2 mm dia.)
L_CRY  = 0.1e-3        # crystal length, m
DEFF   = 1.75e-12      # effective nonlinearity, m/V (paper value)

THETA_S_DEG = 3.0      # detector lab (external) angle, deg
DIST_DET    = 625e-3   # crystal -> detector distance, m
APERTURE    = 6e-3     # detector opening DIAMETER, m

FILTER_FWHM   = 10e-9  # Gaussian bandpass FWHM, m
FILTER_CENTER = None   # centre wavelength, m; None -> degenerate 2*LAM_P

LAM_I_MAX = com.LAM_I_MAX          # idler transparency cutoff (3500 nm)

# ----------------------------------------------------------------------------
# EFFICIENCIES -- datasheet/measured values (2026-07-05).  PP above is the
# LASER OUTPUT power; the pump power reaching the crystal is PP * ETA_PUMP.
# ----------------------------------------------------------------------------
# pump path: laser -> mirror F01 -> half-wave plate -> mirror F01 -> crystal
ETA_MIRROR_F01 = 0.920  # F01 mirror reflectance @ 405 nm
N_MIRROR_F01   = 2      # number of F01 mirrors in the pump path
ETA_HWP        = 0.999  # half-wave plate transmission @ 405 nm
T_CRYSTAL_IN   = 1.0    # BBO entrance face: AR coated, no loss

ETA_PUMP = ETA_MIRROR_F01 ** N_MIRROR_F01 * ETA_HWP * T_CRYSTAL_IN

# collection path, per detector arm (both arms identical):
#   crystal exit -> mirror P01 -> bandpass -> lens+fiber mount -> fiber -> APD
T_CRYSTAL_OUT      = 1.0   # BBO exit face: AR coated, no loss
ETA_MIRROR_P01     = 0.958 # P01 mirror reflectance @ 810 nm
T_FILTER_PEAK      = 0.983 # bandpass PEAK transmission (the Gaussian shape in
                           #   filter_grid stays normalised to 1; this scales it)
ETA_FIBER_COUPLING = 0.80  # lens -> fiber -> detector box, all-in.  Assumes an
                           #   AR-coated lens; budget: 15% loss coupling from
                           #   the lens into the fiber core, 3.5% Fresnel loss
                           #   at the fiber input, negligible loss in the fiber
                           #   itself, 5% loss at the mating with the detector
                           #   box  (0.85 * 0.965 * 0.95 = 0.78, rounded to 0.80)
ETA_FIBER          = 1.0   # propagation loss negligible; connector/mating
                           #   losses are already in ETA_FIBER_COUPLING
ETA_APD            = 0.30  # APD detection efficiency @ 810 nm

# collection chain WITHOUT the filter peak (so an open/unfiltered arm, as
# allowed in coincidence_counts.py, is handled correctly there)
ETA_COLLECT_BASE = (T_CRYSTAL_OUT * ETA_MIRROR_P01 * ETA_FIBER_COUPLING
                    * ETA_FIBER * ETA_APD)
ETA_COLLECT = ETA_COLLECT_BASE * T_FILTER_PEAK   # full chain incl. filter peak

PHYS_KW = dict(lam_p=LAM_P, Pp=PP * ETA_PUMP, L=L_CRY, W=W_PUMP, deff=DEFF)

OUTFILE = OUTDIR / "detector_scan_mysetup.png"


# ----------------------------------------------------------------------------
# Ideal cut angle: degenerate cone exactly through the detector
# ----------------------------------------------------------------------------
def ideal_cut_deg(theta_target_deg, lam_signal):
    """Internal phase-matching angle theta_m (deg) whose perfect-phase-matching
    cone at lam_signal sits exactly at the external angle theta_target_deg.
    Bisection: the cone opening angle grows monotonically with theta_m."""
    lam   = np.array([lam_signal])
    lam_i = lam_signal * LAM_P / (lam_signal - LAM_P)
    ks    = phys.n_o(lam_signal) * 2 * np.pi / lam_signal
    ki    = phys.n_o(lam_i) * 2 * np.pi / lam_i

    def ext(tm_deg):
        kp = phys.n_eff(LAM_P, np.radians(tm_deg)) * 2 * np.pi / LAM_P
        if ks + ki < kp:                  # below the collinear threshold
            return -1.0
        v = phys.tuning_curve(lam, np.radians(tm_deg), lam_p=LAM_P)[0]
        return 91.0 if np.isnan(v) else v # NaN here = beyond max opening angle

    lo, hi = 20.0, 40.0
    if not (ext(lo) < theta_target_deg < ext(hi)):
        raise ValueError("detector angle not reachable by any cut in 20-40 deg")
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if ext(mid) < theta_target_deg:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ----------------------------------------------------------------------------
# Rate through the circular aperture
# ----------------------------------------------------------------------------
def detected_rate(theta_m_deg, lam, T, dlam, Nth=81, Nxi=200):
    """(aperture, ring): filter-weighted photons/s through the 6 mm opening of
    one detector, and the same integrated over the full azimuthal ring within
    the aperture's polar span (reference)."""
    r, D = APERTURE / 2, DIST_DET
    R0 = D * np.tan(np.radians(THETA_S_DEG))
    if R0 <= r:
        raise ValueError("aperture contains the pump axis; azimuthal-arc "
                         "model needs D*tan(theta_s) > aperture radius")
    th  = np.linspace(np.arctan((R0 - r) / D), np.arctan((R0 + r) / D), Nth)
    rho = D * np.tan(th)
    cosarg = (rho**2 + R0**2 - r**2) / (2 * rho * R0)
    dphi   = 2 * np.arccos(np.clip(cosarg, -1.0, 1.0))

    # counts/s per polar radian into the full ring, filter-weighted
    dNdth = (phys.Ns(lam[:, None], th[None, :], np.radians(theta_m_deg),
                     dlam, 1.0, Nxi=Nxi, lam_i_max=LAM_I_MAX, **PHYS_KW)
             * T[:, None]).sum(axis=0)

    aperture = np.trapezoid(dNdth * dphi / (2 * np.pi), x=th)
    ring     = np.trapezoid(dNdth, x=th)
    return ETA_COLLECT * aperture, ETA_COLLECT * ring


def filter_grid():
    """Wavelength grid spanning the bandpass and its Gaussian transmission."""
    center = FILTER_CENTER if FILTER_CENTER is not None else 2 * LAM_P
    lam = np.arange(center - 2.5 * FILTER_FWHM, center + 2.5 * FILTER_FWHM,
                    FILTER_FWHM / 100)
    T = np.exp(-4 * np.log(2) * ((lam - center) / FILTER_FWHM) ** 2)
    return lam, T, center


# ----------------------------------------------------------------------------
def main():
    lam, T, center = filter_grid()
    dlam = lam[1] - lam[0]

    cut = ideal_cut_deg(THETA_S_DEG, center)
    n_det, n_ring = detected_rate(cut, lam, T, dlam)

    r, D = APERTURE / 2, DIST_DET
    print("--- Setup ---")
    print(f"  indices: {phys.sellmeier_label()}")
    print(f"  pump {LAM_P*1e9:.0f} nm, {PP*1e3:.0f} mW, W = {W_PUMP*1e3:.1f} mm"
          f" (1/e^2 radius);  BBO L = {L_CRY*1e3:.2f} mm, deff = {DEFF*1e12:.2f} pm/V")
    print(f"  detector: theta_s = {THETA_S_DEG:.2f} deg, D = {D*1e3:.0f} mm, "
          f"opening {APERTURE*1e3:.0f} mm  (polar acceptance "
          f"{APERTURE/D*1e3:.1f} mrad, ring fraction "
          f"{2*r/(2*np.pi*D*np.tan(np.radians(THETA_S_DEG))):.3f})")
    print(f"  filter: Gaussian {FILTER_FWHM*1e9:.0f} nm FWHM at {center*1e9:.1f} nm")
    print(f"  ideal cut angle theta_m = {cut:.3f} deg "
          f"(degenerate cone through {THETA_S_DEG:.1f} deg at normal incidence)")
    print("--- Efficiencies ---")
    print(f"  pump chain (2x F01 mirror, HWP, crystal face): "
          f"ETA_PUMP = {ETA_PUMP:.3f}  ->  "
          f"{PP*ETA_PUMP*1e3:.1f} mW at the crystal")
    print(f"  collection chain (crystal face, P01 mirror, filter peak, "
          f"fiber coupling, fiber, APD): ETA_COLLECT = {ETA_COLLECT:.3f}")
    print("--- Result (crystal at normal incidence) ---")
    print(f"  photons on ONE detector : {n_det:.3e} /s")
    print(f"  full azimuthal ring     : {n_ring:.3e} /s "
          f"(aperture/ring = {n_det/n_ring:.4f})")

    # --- crystal-rotation scan: same detected rate vs rotation angle alpha ---
    tms    = np.arange(cut - 6.0, cut + 6.0 + 1e-9, 0.05)
    alphas = phys.crystal_rotation(tms, cut, lam_p=LAM_P)
    a2t = lambda a: np.interp(a, alphas, tms)
    t2a = lambda t: np.interp(t, tms, alphas)

    # coarser grids for the sweep (converged to <1% vs the headline number)
    lam_c = lam[::4]
    T_c   = T[::4]
    counts = np.array([detected_rate(t, lam_c, T_c, dlam * 4,
                                     Nth=41, Nxi=200)[0] for t in tms])

    ipk, thresh = com.summarize(
        f"my setup: {FILTER_FWHM*1e9:.0f} nm FWHM bandpass, one detector",
        alphas, tms, counts)

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    com.draw_panel(
        ax, alphas, counts, ipk, thresh,
        ylabel="photons/s on detector",
        title=(f"Detector at {THETA_S_DEG:.0f}$^\\circ$, "
               f"⌀{APERTURE*1e3:.0f} mm @ {DIST_DET*1e3:.0f} mm  --  "
               f"L = {L_CRY*1e3:.1f} mm, {PP*1e3:.0f} mW, "
               f"W = {W_PUMP*1e3:.0f} mm, "
               f"{FILTER_FWHM*1e9:.0f} nm FWHM @ {center*1e9:.0f} nm"),
        subtitle=(f"cut at ideal angle $\\theta_m$ = {cut:.2f}$^\\circ$ for "
                  f"{THETA_S_DEG:.0f}$^\\circ$; circular-aperture azimuthal "
                  f"arc included"),
        a2t=a2t, t2a=t2a, subtitle_va="top")
    ax.title.set_fontsize(10)            # long parameter line; keep it inside
    fig.tight_layout()
    fig.savefig(OUTFILE, dpi=200)
    print(f"wrote {OUTFILE}")


if __name__ == "__main__":
    from textlog import tee_stdout
    with tee_stdout("detector_scan_mysetup.txt"):
        main()
