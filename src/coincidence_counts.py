"""
coincidence_counts.py -- signal-idler COINCIDENCE rate between two real,
finite-aperture detectors at +theta_s and -theta_s, for THIS setup (the one in
detector_counts.py, whose pump/crystal parameters and ideal-cut solver are
reused; detector_counts.py itself is untouched).

Physics
-------
phys.Ns is a MARGINAL: its xi_i quadrature integrates the idler over ALL
transverse momenta.  The integrand exp(-(xi_s - xi_i)^2/2) * sinc^2(L dkz / 2)
is the JOINT signal-idler transverse-momentum density, so the coincidence rate
is the same integral with

  (a) the xi_i quadrature restricted to detector B's polar acceptance
      (phys.Ns(..., idler_theta_window=...), added for this script;
      default None keeps every existing caller bit-identical),
  (b) the wavelength grid weighted by BOTH filters, T_A(lam_s) * T_B(lam_i)
      with lam_i = lam_s lam_p / (lam_s - lam_p),
  (c) the azimuthal arc factor dphi_A/2pi of the singles replaced by the
      overlap of arm A's arc with arm B's arc mirrored through the pump axis,
      smeared by the azimuthal momentum correlation width
      sigma_phi = (1/W) / (k0 sin theta)  (~0.14 deg here, vs a ~10.5 deg
      full arc -- so geometry costs only a percent-level rim loss).

With identical near-degenerate Gaussian filters in both arms the spectral
factor is T(lam)^2 (dlam_i/dlam_s = -1 at degeneracy), so the heralding ratio
is ~ integral(T^2)/integral(T) = 1/sqrt(2) ~ 0.71 -- the filters, not the
geometry, set the coincidence/singles ratio in this loss-free model.  Real
detector quantum efficiencies and optical losses multiply on top.

Validation (printed on every run): opening arm B's polar window and filter and
using arm A's plain arc must reproduce the singles rate of detector_counts.py.

Outputs: prints singles for each arm, coincidences, heralding ratios and the
accidental rate for the coincidence window TAU_C, and writes
coincidence_scan.png (singles + coincidences + heralding ratio vs crystal
rotation angle alpha).

Run:  python coincidence_counts.py
"""

from pathlib import Path

import numpy as np
from math import erf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import spdc_physics as phys
import detector_scan_common as com
import detector_counts as dc          # pump/crystal parameters, ideal_cut_deg

OUTDIR = Path(__file__).resolve().parent.parent / "output"

# ----------------------------------------------------------------------------
# PARAMETERS -- the two detector arms.  theta_deg is the magnitude of the lab
# angle; the model places arm B diametrically opposite arm A (i.e. at -theta,
# azimuth phi + pi), which is where momentum conservation sends the idler.
# filter_fwhm=None means an open (unfiltered) arm; filter_center=None means
# centred on the degenerate wavelength 2*lam_p.
# ----------------------------------------------------------------------------
ARM_A = dict(theta_deg=3.0, dist=625e-3, aperture=6e-3,
             filter_fwhm=10e-9, filter_center=None)
ARM_B = dict(theta_deg=3.0, dist=625e-3, aperture=6e-3,
             filter_fwhm=10e-9, filter_center=None)

TAU_C = 1e-9            # coincidence window (s), for the accidentals estimate

LAM_P     = dc.LAM_P
PHYS_KW   = dc.PHYS_KW
LAM_I_MAX = dc.LAM_I_MAX
SIGMA_KAPPA = 1.0 / dc.W_PUMP   # transverse-momentum correlation width, 1/m

OUTFILE = OUTDIR / "coincidence_scan.png"


# ----------------------------------------------------------------------------
# Geometry helpers (same circular-aperture model as detector_counts.py)
# ----------------------------------------------------------------------------
def arc_geometry(arm, Nth):
    """Polar grid across the aperture and the azimuthal arc dphi(theta)."""
    r, D = arm["aperture"] / 2, arm["dist"]
    R0 = D * np.tan(np.radians(arm["theta_deg"]))
    if R0 <= r:
        raise ValueError("aperture contains the pump axis; azimuthal-arc "
                         "model needs D*tan(theta) > aperture radius")
    th  = np.linspace(np.arctan((R0 - r) / D), np.arctan((R0 + r) / D), Nth)
    rho = D * np.tan(th)
    cosarg = (rho**2 + R0**2 - r**2) / (2 * rho * R0)
    dphi   = 2 * np.arccos(np.clip(cosarg, -1.0, 1.0))
    return th, dphi


def polar_window(arm):
    """(lo, hi) external polar angles subtended by the aperture, rad."""
    r, D = arm["aperture"] / 2, arm["dist"]
    R0 = D * np.tan(np.radians(arm["theta_deg"]))
    return (np.arctan((R0 - r) / D), np.arctan((R0 + r) / D))


def transmission(arm, lam):
    """Filter transmission of one arm at wavelength(s) lam."""
    if arm["filter_fwhm"] is None:
        return np.ones_like(lam)
    c0 = arm["filter_center"] if arm["filter_center"] is not None else 2 * LAM_P
    return np.exp(-4 * np.log(2) * ((lam - c0) / arm["filter_fwhm"]) ** 2)


def azimuthal_overlap(th, dphi_a, arm_b, Nphi=201):
    """Coincidence azimuthal weight replacing the singles' dphi_a/(2 pi).

    A signal at azimuth phi inside arm A's arc sends its twin to phi + pi
    + delta, delta ~ N(0, sigma_phi) with sigma_phi = SIGMA_KAPPA/(k0 sin th).
    Weight = (1/2pi) INT_{A arc} dphi P(twin inside B's mirrored arc).
    B's arc is evaluated at the same polar angle (the polar mirror is the
    identity for near-degenerate pairs; the <=0.04 deg wavelength walk-off is
    handled exactly in the polar window, and is negligible azimuthally)."""
    r, D = arm_b["aperture"] / 2, arm_b["dist"]
    R0  = D * np.tan(np.radians(arm_b["theta_deg"]))
    rho = D * np.tan(th)
    cosarg = (rho**2 + R0**2 - r**2) / (2 * rho * R0)
    dphi_b = 2 * np.arccos(np.clip(cosarg, -1.0, 1.0))

    kappa_ring = (2 * np.pi / (2 * LAM_P)) * np.sin(th)   # degenerate ring
    sig_phi = SIGMA_KAPPA / kappa_ring
    verf = np.vectorize(erf)

    out = np.zeros_like(th)
    for j in range(len(th)):
        if dphi_a[j] == 0.0 or dphi_b[j] == 0.0:
            continue
        phi = np.linspace(-dphi_a[j] / 2, dphi_a[j] / 2, Nphi)
        s   = np.sqrt(2) * sig_phi[j]
        P   = 0.5 * (verf((dphi_b[j] / 2 - phi) / s)
                     + verf((dphi_b[j] / 2 + phi) / s))
        out[j] = np.trapezoid(P, x=phi) / (2 * np.pi)
    return out


# ----------------------------------------------------------------------------
# Rates
# ----------------------------------------------------------------------------
def singles_rate(arm, theta_m_deg, lam, w, dlam, Nth=81, Nxi=200):
    """Filter-weighted photons/s through one arm's aperture (as in
    detector_counts.detected_rate, but for an arbitrary arm)."""
    th, dphi = arc_geometry(arm, Nth)
    dNdth = (phys.Ns(lam[:, None], th[None, :], np.radians(theta_m_deg),
                     dlam, 1.0, Nxi=Nxi, lam_i_max=LAM_I_MAX, **PHYS_KW)
             * w[:, None]).sum(axis=0)
    return np.trapezoid(dNdth * dphi / (2 * np.pi), x=th)


def coincidence_geometry(Nth=81):
    """Precompute the theta grid, azimuthal weight and idler polar window
    (geometry only -- independent of the crystal angle)."""
    th, dphi_a = arc_geometry(ARM_A, Nth)
    return th, azimuthal_overlap(th, dphi_a, ARM_B), polar_window(ARM_B)


def coincidence_rate(theta_m_deg, lam, w_joint, dlam, geom, Nxi=200):
    """Coincidences/s: joint (signal in A) x (idler in B) rate.
    w_joint = T_A(lam_s) * T_B(lam_i(lam_s))."""
    th, wphi, win = geom
    dNdth = (phys.Ns(lam[:, None], th[None, :], np.radians(theta_m_deg),
                     dlam, 1.0, Nxi=Nxi, lam_i_max=LAM_I_MAX,
                     idler_theta_window=win, **PHYS_KW)
             * w_joint[:, None]).sum(axis=0)
    return np.trapezoid(dNdth * wphi, x=th)


def wavelength_grids():
    """Signal grid + per-arm weights.  Grid spans the tightest filter present
    (fallback: broadband 480-1300 nm; only meaningful with lam_i_max)."""
    fwhms = [a["filter_fwhm"] for a in (ARM_A, ARM_B)
             if a["filter_fwhm"] is not None]
    center = 2 * LAM_P
    if fwhms:
        fw  = min(fwhms)
        lam = np.arange(center - 2.5 * fw, center + 2.5 * fw, fw / 100)
    else:
        lam = np.arange(480e-9, 1300e-9, 0.5e-9)
    lam_i  = lam * LAM_P / (lam - LAM_P)          # conjugate idler wavelength
    w_a    = transmission(ARM_A, lam)             # A sees the signal photon
    w_b    = transmission(ARM_B, lam)             # B's own singles weight
    w_conj = transmission(ARM_B, lam_i)           # B sees the idler twin
    return lam, lam[1] - lam[0], w_a, w_b, w_conj


# ----------------------------------------------------------------------------
def main():
    lam, dlam, w_a, w_b, w_conj = wavelength_grids()
    center = (ARM_A["filter_center"] if ARM_A["filter_center"] is not None
              else 2 * LAM_P)
    cut = dc.ideal_cut_deg(ARM_A["theta_deg"], center)

    geom = coincidence_geometry()
    s_a = singles_rate(ARM_A, cut, lam, w_a, dlam)
    s_b = singles_rate(ARM_B, cut, lam, w_b, dlam)
    r_c = coincidence_rate(cut, lam, w_a * w_conj, dlam, geom)

    # Validation: open arm B completely (full polar window, no filter, plain
    # arc-A azimuthal factor) -> must recover the singles marginal of arm A.
    th, _, _ = geom
    _, dphi_a = arc_geometry(ARM_A, len(th))
    geom_open = (th, dphi_a / (2 * np.pi), (1e-4, np.radians(89.0)))
    s_check = coincidence_rate(cut, lam, w_a, dlam, geom_open)

    print("--- Two-detector setup (arms at +/-%.2f deg) ---" % ARM_A["theta_deg"])
    print(f"  pump {LAM_P*1e9:.0f} nm, {dc.PP*1e3:.0f} mW, W = {dc.W_PUMP*1e3:.1f} mm; "
          f"BBO L = {dc.L_CRY*1e3:.2f} mm; ideal cut theta_m = {cut:.3f} deg")
    for name, arm in (("A", ARM_A), ("B", ARM_B)):
        filt = ("open" if arm["filter_fwhm"] is None else
                f"{arm['filter_fwhm']*1e9:.0f} nm FWHM Gaussian")
        print(f"  arm {name}: {arm['theta_deg']:.2f} deg, "
              f"D = {arm['dist']*1e3:.0f} mm, "
              f"opening {arm['aperture']*1e3:.0f} mm, filter {filt}")
    print("--- Rates (crystal at normal incidence, unit detector efficiency) ---")
    print(f"  singles arm A          : {s_a:.3e} /s")
    print(f"  singles arm B          : {s_b:.3e} /s")
    print(f"  coincidences           : {r_c:.3e} /s")
    print(f"  heralding  Rc/Rs_A     : {r_c/s_a:.3f}   (Rc/Rs_B: {r_c/s_b:.3f})")
    print(f"  accidentals Ra*Rb*tau  : {s_a*s_b*TAU_C:.2e} /s  (tau = {TAU_C*1e9:.1f} ns)")
    print(f"  [check] open-B recovers singles: {s_check:.3e} /s "
          f"(rel. diff {abs(s_check - s_a)/s_a:.1e})")

    # --- crystal-rotation scan --------------------------------------------
    tms = np.arange(cut - 6.0, cut + 6.0 + 1e-9, 0.05)
    alphas, a2t, t2a = com.rotation_axes(tms, cut)

    lam_c, dlam_c = lam[::4], dlam * 4       # converged coarse grids, as in
    wa_c, wb_c, wj_c = w_a[::4], w_b[::4], (w_a * w_conj)[::4]   # detector_counts
    geom_c = coincidence_geometry(Nth=41)

    s_scan = np.array([singles_rate(ARM_A, t, lam_c, wa_c, dlam_c, Nth=41)
                       for t in tms])
    c_scan = np.array([coincidence_rate(t, lam_c, wj_c, dlam_c, geom_c)
                       for t in tms])

    ipk_s, _ = com.summarize("singles, arm A", alphas, tms, s_scan)
    ipk_c, thresh_c = com.summarize("coincidences", alphas, tms, c_scan)

    fig, (ax, axr) = plt.subplots(2, 1, figsize=(7.2, 6.0), sharex=True,
                                  height_ratios=[2.2, 1.0])
    ax.semilogy(alphas, s_scan, color="C0", lw=1.6, label="singles (one arm)")
    ax.semilogy(alphas, c_scan, color="C1", lw=1.6, label="coincidences")
    ax.axhline(thresh_c, color="0.5", ls=":", lw=1)
    ax.axvline(alphas[ipk_c], color="C3", ls="--", lw=1)
    ax.set_ylabel("counts/s")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="lower left", fontsize=8)
    ax.set_ylim(c_scan.max() / 3e3, s_scan.max() * 2)
    ax.set_title(f"Two detectors at $\\pm${ARM_A['theta_deg']:.0f}$^\\circ$, "
                 f"⌀{ARM_A['aperture']*1e3:.0f} mm @ {ARM_A['dist']*1e3:.0f} mm"
                 f"  --  L = {dc.L_CRY*1e3:.1f} mm, {dc.PP*1e3:.0f} mW, "
                 f"{ARM_A['filter_fwhm']*1e9:.0f} nm FWHM both arms",
                 fontsize=10, pad=26)
    sec = ax.secondary_xaxis("top", functions=(a2t, t2a))
    sec.set_xlabel(r"internal phase-matching angle $\theta_m$ (deg)")

    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(s_scan > 0, c_scan / s_scan, np.nan)
    axr.plot(alphas, ratio, color="C2", lw=1.6)
    axr.axhline(1 / np.sqrt(2), color="0.5", ls=":", lw=1,
                label=r"$1/\sqrt{2}$ (matched-filter limit)")
    axr.set_ylabel(r"heralding $R_c/R_s$")
    axr.set_ylim(0, 1)
    axr.grid(True, alpha=0.25)
    axr.legend(loc="lower left", fontsize=8)
    axr.set_xlabel(r"crystal rotation angle $\alpha$ (deg, 0 = normal incidence)")

    fig.tight_layout()
    fig.savefig(OUTFILE, dpi=200)
    print(f"wrote {OUTFILE}")


if __name__ == "__main__":
    main()
