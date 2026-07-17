"""
spdc_physics.py -- physics library for Type-I SPDC in BBO.

Numerical recreation of the angle-resolved parametric-fluorescence flux from
Hsu & Lai, "Angle-Resolved Spectroscopy of Parametric Fluorescence"
(arXiv:1302.6226).  Type-I (e -> o + o), pump 405 nm, degenerate signal 810 nm.

This module is the SINGLE source of the physics: constants, experiment
parameters, the BBO Sellmeier indices, the flux master equation Ns, the
perfect-phase-matching tuning curve, and the crystal-rotation refraction
relation.  The plotting scripts (plot_fig6.py, plot_detector_scan*.py) import
it and contain no physics of their own.

Running this file directly performs the sanity check against the paper's
Table I / Fig. 6 numbers:  python spdc_physics.py


MASTER EQUATION -- Eq. (9), angle-resolved signal flux
------------------------------------------------------
Ns() below evaluates the paper's Appendix intermediate equation directly:

  dNs = [hbar deff^2 w_s w_i w_p L^2 Np sqrt(2pi) / (8 pi^4 c^3 eps0 ns ni np)]
        * dw_s * (kappa_s dkappa_s) * dphi * INT dxi_i (...)

changing variables w_s -> lam_s and kappa_s -> theta and integrating the
signal azimuth phi over 2*pi, with the Jacobians

  |dw_s/dlam_s| = 2 pi c / lam_s^2 ,
  kappa_s dkappa_s = (k_perp^2/2) sin(2 theta) dtheta .

In the paper's own variables (theta_s = EXTERNAL angle, see below) this equals

  Ns(lam_s, th_s) = 2*pi * [ 2*pi*sqrt(2*pi) * hbar*deff^2*omega_p*L^2*Np
                    / (eps0 * ns*ni*np * lam_s^5 * lam_i)
                    * sin(2 th_s) * dlam * dth
                    * INT dxi_i exp(-1/2 (xi_s - xi_i)^2) * sinc^2(1/2 L dkz) ]

i.e. exactly the printed Eq. (9) times the 2*pi azimuth ring, with
sinc(x) = sin(x)/x, the longitudinal phase mismatch (Appendix form)

  dkz  = sqrt(ks^2 - kappa_s^2) + sqrt(ki^2 - (xi_i/W)^2) - kp ,
  k    = 2*pi*n/lam  (in-medium),   kappa_s = k0*sin(th_s),  k0 = 2*pi/lam_s ,
  xi_s = kappa_s*W ,   xi_i = kappa_i*W    (dimensionless transverse momenta),

and energy conservation fixing the idler:  lam_i = lam_s*lam_p/(lam_s - lam_p).

Indices: signal/idler ordinary ns = no(lam_s), ni = no(lam_i); pump effective
extraordinary np = n_eff(theta_m, lam_p) (Eq. 4).


PAPER ERRATA STATUS (re-checked 2026-07-13; derivation in spdc_eq9_note.tex)
----------------------------------------------------------------------------
The printed Eq. (9) is CORRECT: theta_s is the EXTERNAL (lab) angle -- Fig. 2
defines theta_s' (internal) vs theta_s (external), and unprimed ks in the
appendix's "kappa_s = ks sin(theta_s)" is the VACUUM wavenumber (the paper
primes in-medium quantities: ks' = ns ws/c).  Then kappa_s dkappa_s =
(k0^2/2) sin(2 th) dth carries no ns^2, and the printed prefactor
2*pi*sqrt(2*pi) is exact -- as a density per unit azimuth (Eq. 9 merely drops
the dphi symbol; Sec. III.D integrates phi = 0..2*pi explicitly).  An earlier
version of this code's docs claimed a "missing 2*pi*ns^2 ~ 17x" typo in
Eq. (9); that came from misreading theta_s as internal and is RETRACTED.
The genuine misprints are confined to the main-text Eq. (7):
1. Factor 4: Eq.(7)'s prefactor is 2 pi^4 in the main text but 8 pi^4 in the
   Appendix; 8 pi^4 is correct (reproduces Table I's simulated eta ~ 2.6e-10,
   the independent closed-form Eq. (8) = 2.63e-10, and the measured
   2.7-2.8e-10; this code gives 2.72e-10).
2. Sign: Eq.(7)'s exponent prints exp(+xi^2/2); must be -xi^2/2 (as in the
   Appendix).


ANGLE CONVENTION
----------------
The paper's Eq. (9) (and Fig. 6) use the EXTERNAL (air/lab) angle:
kappa_s = k0 sin theta_ext, k0 = 2pi/lam (Fig. 2 convention; kappa is
conserved across the exit face, so this equals ks sin theta_int with
ks = ns k0).  external=True (default) implements exactly that.
external=False parametrises the same kappa_s by the INTERNAL angle
(kappa_s = ks sin theta_int) -- an equivalent change of variables whose
Jacobian carries ns^2.  Both give identical TOTAL flux -- only the y-axis
scale of a per-angle plot differs by ~ns; this is checked in sanity_check().
"""

import numpy as np

pi = np.pi

# ----------------------------------------------------------------------------
# Physical constants (SI)
# ----------------------------------------------------------------------------
hbar = 1.054571817e-34      # J s
eps0 = 8.8541878188e-12     # F/m
c    = 2.99792458e8         # m/s

# ----------------------------------------------------------------------------
# Experiment / crystal parameters (paper Fig. 6)
# ----------------------------------------------------------------------------
deff  = 1.75e-12            # effective nonlinearity, m/V (1.75 pm/V)
lam_p = 405e-9              # pump wavelength, m
L     = 3e-3                # crystal length, 3 mm
Pp    = 80e-3               # pump power, 80 mW
W     = np.sqrt(60e-6 * 30e-6)   # pump 1/e^2 radius: sqrt(60um x 30um) ~ 42.4 um

omega_p = 2 * pi * c / lam_p
Np      = Pp / (hbar * omega_p)  # pump photon flux, /s (paper Fig. 6: 1.63e17)


# ----------------------------------------------------------------------------
# BBO refractive indices -- Sellmeier equations (lambda in MICROMETRES inside)
#
# Two selectable coefficient sets ("references"):
#   ref 1, "kato"   -- Kato 1986 (IEEE JQE 22, 1013).  DEFAULT.  More accurate
#                      at 405/810 nm: agrees with the modern measurement
#                      (Tamosauskas et al. 2018, Opt. Mater. Express 8, 1410)
#                      to ~1-2e-4 there, vs 5-10e-4 for Eimerl.
#   ref 2, "eimerl" -- Eimerl et al. 1987 (J. Appl. Phys. 62, 1968).  What the
#                      paper used, via the NIST noncollinear phase-matching
#                      program (paper refs [17]/[25]); select it to reproduce
#                      the paper's figures/labels (Fig. 6) exactly.
# n^2 = A + B/(lam^2 - C) - D*lam^2, lam in um.
# ----------------------------------------------------------------------------
SELLMEIER_REFS = {
    "kato":   {"no": (2.7359, 0.01878, 0.01822, 0.01354),
               "ne": (2.3753, 0.01224, 0.01667, 0.01516),
               "label": "Sellmeier ref 1: Kato 1986"},
    "eimerl": {"no": (2.7405, 0.0184, 0.0179, 0.0155),
               "ne": (2.3730, 0.0128, 0.0156, 0.0044),
               "label": "Sellmeier ref 2: Eimerl 1987 (paper/NIST)"},
}
_SELLMEIER_ALIASES = {"kato": "kato", "eimerl": "eimerl", 1: "kato", 2: "eimerl",
                      "1": "kato", "2": "eimerl"}
_sellmeier = "kato"


def use_sellmeier(ref):
    """Select the Sellmeier coefficient set: 1/"kato" (default) or
    2/"eimerl" (the paper's).  Returns the human-readable label."""
    global _sellmeier
    key = _SELLMEIER_ALIASES.get(ref if isinstance(ref, int) else str(ref).lower())
    if key is None:
        raise ValueError(f"unknown Sellmeier ref {ref!r}; use 1/'kato' or 2/'eimerl'")
    _sellmeier = key
    return sellmeier_label()


def sellmeier_label():
    """Label of the active Sellmeier set (for printing into script outputs)."""
    return SELLMEIER_REFS[_sellmeier]["label"]


def _sellmeier_n(lam, coeffs):
    l2 = (lam * 1e6) ** 2
    A, B, C, D = coeffs
    return np.sqrt(A + B / (l2 - C) - D * l2)


def n_o(lam):
    """Ordinary index of BBO (active Sellmeier set).  lam in metres."""
    return _sellmeier_n(lam, SELLMEIER_REFS[_sellmeier]["no"])


def n_e_principal(lam):
    """Principal extraordinary index of BBO (propagation at 90 deg to the
    optic axis, active Sellmeier set).  lam in metres."""
    return _sellmeier_n(lam, SELLMEIER_REFS[_sellmeier]["ne"])


def n_eff(lam, theta_m):
    """Angle-dependent extraordinary index (paper Eq. 4) seen by the pump
    e-wave propagating at angle theta_m (rad) to the optic axis:
        n_eff^-2 = cos^2(theta_m)/n_o^2 + sin^2(theta_m)/n_e^2 ."""
    return 1.0 / np.sqrt(np.cos(theta_m)**2 / n_o(lam)**2
                         + np.sin(theta_m)**2 / n_e_principal(lam)**2)


def sinc(x):
    """sin(x)/x.  (NumPy's np.sinc is the engineering sinc sin(pi x)/(pi x),
    so we divide the argument by pi.)"""
    return np.sinc(x / pi)


# ----------------------------------------------------------------------------
# Master equation: photon flux Ns per (dlam, dtheta) bin
# ----------------------------------------------------------------------------
def Ns(lam_s, theta_s, theta_m, dlam, dtheta, *, external=True,
       Nxi=200, half_window=8.0, lam_i_max=None, idler_theta_window=None,
       lam_p=lam_p, Pp=Pp, L=L, W=W, deff=deff):
    """Signal photon flux (counts/s) into one (dlam x dtheta) bin, integrated
    over the full 2*pi azimuth of the emission cone.

    Evaluates the Appendix intermediate equation (= the printed Eq. 9, which
    is a per-unit-azimuth density in the external angle, times the 2*pi
    ring -- see the module docstring for the errata discussion).

    The experiment parameters (lam_p, Pp, L, W, deff) can be overridden per
    call; they default to the module-level paper values, so existing callers
    are unaffected.

    Parameters
    ----------
    lam_s   : signal wavelength(s), m.  Scalar or array.
    theta_s : signal polar angle(s), rad.  Scalar or array; may be signed
              (the cone is symmetric, only |theta_s| enters the physics).
              lam_s and theta_s are numpy-broadcast against each other, so you
              can pass e.g. a lam vector at fixed theta (detector scan) or a
              theta vector at fixed lam (one column of a Fig. 6 map).
    theta_m : pump angle to the optic axis (phase-matching angle), rad. Scalar.
    dlam    : wavelength bin width, m.
    dtheta  : polar-angle bin width, rad.
    external: True  -> theta_s is the EXTERNAL (air) angle,
                       kappa_s = k0 sin(theta_s) with k0 = 2 pi/lam_s
                       (the paper's Eq. 9 / Fig. 6 convention, per Fig. 2);
              False -> theta_s is the INTERNAL angle,
                       kappa_s = ks sin(theta_s) with ks = ns k0
                       (equivalent change of variables).
              Both give the same total flux; see sanity_check().
    Nxi     : number of points for the xi_i (idler transverse momentum)
              quadrature.
    half_window : xi_i integration window, +/- half_window around xi_s.  The
              integrand carries exp(-(xi_s - xi_i)^2 / 2), so +/-8 sigma-units
              is effectively exact.
    lam_i_max : if given, zero the flux where the idler wavelength exceeds this
              (e.g. 3500e-9 = BBO transparency edge / limit of the index data).
              Without it, short-lam_s results are propped up by photon pairs
              that cannot actually be generated.
    idler_theta_window : optional (lo, hi) in rad.  If given, the xi_i
              quadrature only counts idlers whose EXTERNAL polar angle
              arcsin(kappa_i / k0_i) lies inside [lo, hi], with
              kappa_i = xi_i / W and k0_i = 2 pi / lam_i (Snell: kappa is
              conserved across the exit face).  Angles are measured on the
              OPPOSITE side of the pump from the signal (the model's Gaussian
              exp(-(xi_s - xi_i)^2 / 2) peaks at xi_i = xi_s, i.e. at the
              momentum-conserving mirror direction), so a detector at -3 deg
              facing a signal at +3 deg is window ~ (2.7, 3.3) deg.  This
              turns Ns into the JOINT rate "signal in this (dlam x dtheta)
              bin AND idler inside the window" -- the coincidence kernel.
              Idlers past total internal reflection (kappa_i > k0_i) are
              excluded.  Default None integrates over everything and
              reproduces the singles marginal exactly (unchanged code path).

    Returns
    -------
    counts/s per bin, shaped like broadcast(lam_s, theta_s).
    """
    lam_s = np.asarray(lam_s, dtype=float)
    th    = np.abs(np.asarray(theta_s, dtype=float))     # cone symmetry
    lam_s, th = np.broadcast_arrays(lam_s, th)

    # Pump frequency and photon flux for the (possibly overridden) pump
    # parameters; identical to the module-level omega_p, Np at the defaults.
    omega_p = 2 * pi * c / lam_p
    Np      = Pp / (hbar * omega_p)

    # Energy conservation: 1/lam_p = 1/lam_s + 1/lam_i.
    lam_i = lam_s * lam_p / (lam_s - lam_p)

    # Refractive indices: signal/idler ordinary, pump effective extraordinary.
    ns   = n_o(lam_s)
    ni   = n_o(lam_i)
    npmp = n_eff(lam_p, theta_m)

    omega_s = 2 * pi * c / lam_s
    omega_i = 2 * pi * c / lam_i

    k0 = 2 * pi / lam_s                 # vacuum wavenumber of the signal
    ks = ns * k0                        # in-medium wavenumbers
    ki = ni * 2 * pi / lam_i
    kp = npmp * 2 * pi / lam_p

    # Transverse momentum of the detected signal, in the chosen angle
    # convention, and its dimensionless version xi_s = kappa_s * W.
    k_perp  = k0 if external else ks
    kappa_s = k_perp * np.sin(th)
    xi_s    = kappa_s * W

    # --- xi_i quadrature -------------------------------------------------
    # The Gaussian exp(-(xi_s - xi_i)^2 / 2) confines the integrand to a
    # +/- half_window neighbourhood of xi_s, so integrate on a window that
    # tracks xi_s.  A trailing axis of length Nxi is added for the quadrature.
    offs = np.linspace(-half_window, half_window, Nxi)   # (Nxi,)
    xi_i = xi_s[..., None] + offs                        # (..., Nxi)

    # Longitudinal phase mismatch dkz = kz_s + kz_i - kp with
    # kz_s = sqrt(ks^2 - kappa_s^2) (= ks cos th for the internal convention)
    # and kz_i = sqrt(ki^2 - kappa_i^2), kappa_i = xi_i / W.
    long_s = np.sqrt(np.clip(ks**2 - kappa_s**2, 0.0, None))
    arg_i  = ki[..., None]**2 - (xi_i / W)**2            # kz_i^2; <=0 -> evanescent
    dkz    = long_s[..., None] + np.sqrt(np.clip(arg_i, 0.0, None)) - kp

    integrand = np.exp(-0.5 * (xi_s[..., None] - xi_i)**2) * sinc(0.5 * L * dkz)**2
    integrand[arg_i <= 0.0] = 0.0                        # no propagating idler

    if idler_theta_window is not None:
        # Keep only idlers exiting into the external polar window: Snell
        # conserves kappa_i across the face, so sin(theta_i_ext) =
        # kappa_i / k0_i; sin is monotonic on (-pi/2, pi/2) so compare sines.
        lo, hi = idler_theta_window
        sin_ti = (xi_i / W) * (lam_i[..., None] / (2 * pi))
        integrand = np.where((sin_ti >= np.sin(lo)) & (sin_ti <= np.sin(hi)),
                             integrand, 0.0)

    integral = np.trapezoid(integrand, x=offs, axis=-1)

    # --- prefactor ---------------------------------------------------------
    # Appendix intermediate equation ("8 pi^4" form), then change variables:
    #   |domega_s/dlam_s| = 2 pi c / lam_s^2
    #   kappa_s dkappa_s  = (k_perp^2 / 2) sin(2 th) dtheta
    #   integral over azimuth phi -> 2 pi
    A = (hbar * deff**2 * omega_s * omega_i * omega_p * L**2 * Np
         * np.sqrt(2 * pi) / (8 * pi**4 * c**3 * eps0 * ns * ni * npmp))
    pref = A * (2 * pi * c / lam_s**2) * (0.5 * k_perp**2 * np.sin(2 * th)) * (2 * pi)

    out = pref * integral * dlam * dtheta                # counts/s per bin

    if lam_i_max is not None:
        out = np.where(lam_i <= lam_i_max, out, 0.0)
    # Guard against out-of-domain evaluations (lam_s <= lam_p, or lam_i beyond
    # the Sellmeier validity range) leaking NaN/inf into sums.
    return np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)


# ----------------------------------------------------------------------------
# Perfect phase-matching tuning curve  theta_ext(lam_s)
# ----------------------------------------------------------------------------
def tuning_curve(lam_array, theta_m, lam_p=lam_p):
    """External emission angle (deg) of exact phase matching vs signal
    wavelength, for pump angle theta_m (rad).  lam_p overridable per call
    (defaults to the paper pump).  Solves

        sqrt(ks^2 - kappa^2) + sqrt(ki^2 - kappa^2) = kp

    for the common transverse momentum kappa (signal and idler anti-parallel
    transverse, kappa_s = kappa_i = kappa), then kappa = k0 sin(theta_ext).
    Returns NaN where no phase-matching angle exists."""
    kp = n_eff(lam_p, theta_m) * 2 * pi / lam_p
    th_ext = np.full(lam_array.shape, np.nan)
    for idx, lam_s in enumerate(lam_array):
        if lam_s <= lam_p:
            continue
        lam_i = lam_s * lam_p / (lam_s - lam_p)
        k0 = 2 * pi / lam_s
        ks = n_o(lam_s) * k0
        ki = n_o(lam_i) * 2 * pi / lam_i
        if ks + ki < kp:                  # not phase-matchable at any angle
            continue
        lo, hi = 0.0, min(ks, ki) * (1 - 1e-12)
        f = lambda k: np.sqrt(ks**2 - k**2) + np.sqrt(ki**2 - k**2) - kp
        if f(hi) > 0:                     # PM angle exceeds available aperture
            continue
        for _ in range(60):               # bisection (f is decreasing in kappa)
            mid = 0.5 * (lo + hi)
            if f(mid) > 0:
                lo = mid
            else:
                hi = mid
        s = 0.5 * (lo + hi) / k0          # sin(theta_ext)
        if s <= 1.0:
            th_ext[idx] = np.degrees(np.arcsin(s))
    return th_ext


# ----------------------------------------------------------------------------
# Crystal rotation <-> internal phase-matching angle (entrance-face refraction)
# ----------------------------------------------------------------------------
def crystal_rotation(theta_m, theta_m_cut, lam_p=lam_p):
    """External crystal rotation angle alpha (deg, signed) that produces the
    internal phase-matching angle theta_m (deg), for a crystal cut at
    theta_m_cut (deg) (= theta_m at normal incidence, alpha = 0).  lam_p
    overridable per call (defaults to the paper pump).

    The pump refracts at the entrance face; tangential-k matching for the
    extraordinary pump wave gives
        sin(alpha) = n_eff(lam_p, theta_m) * sin(theta_m - theta_m_cut) ,
    so a given d(theta_m) needs a ~n_eff (~1.66x) larger external rotation.

    Note: rotating the crystal also tilts the internal pump (and with it the
    emission-cone axis), but a 2D ray-trace shows pump-in + signal-out
    refraction leaves the internal signal-pump angle for a fixed external
    detector essentially constant (n_o(signal) ~ n_eff(pump) ~ 1.66 in BBO),
    so the fixed-external-angle detector model stays geometrically correct."""
    beta = np.radians(theta_m - theta_m_cut)   # internal deviation from face normal
    npmp = n_eff(lam_p, np.radians(theta_m))
    s = np.clip(npmp * np.sin(beta), -1.0, 1.0)
    return np.degrees(np.arcsin(s))


# ----------------------------------------------------------------------------
# Sanity check against the paper (run:  python spdc_physics.py)
#
# Must hold:  eta = 2 Ns / Np ~ 2.6-2.7e-10 at 810 nm, total degenerate flux
# 2 Ns ~ 4.2e7 /s for theta_m = 29.12 deg, Pp = 80 mW.  Integration region:
# lam_s in [809.5, 810.5] nm, full azimuth (the 2*pi is baked into Ns).
# The external- and internal-angle conventions must give the SAME total.
# ----------------------------------------------------------------------------
def sanity_check():
    theta_m  = np.radians(29.12)
    lam_grid = np.arange(809.5, 810.5 + 1e-6, 0.05) * 1e-9      # m
    dlam = lam_grid[1] - lam_grid[0]

    def total(external, th_max_deg):
        th  = np.arange(0.0, th_max_deg + 1e-9, 0.01) * pi / 180
        dth = th[1] - th[0]
        return Ns(lam_grid[:, None], th[None, :], theta_m, dlam, dth,
                  external=external).sum()

    Ns_ext = total(external=True,  th_max_deg=12.0)   # external angle (Fig. 6 axis)
    Ns_int = total(external=False, th_max_deg=7.5)    # internal angle (Eq. 9 variable)

    print("--- Sanity check (theta_m = 29.12 deg, Pp = 80 mW) ---")
    print(f"  indices          : {sellmeier_label()}")
    print(f"  Np               = {Np:.3e} /s   (paper: 1.63e17)")
    print(f"  n_o(810 nm)      = {n_o(810e-9):.4f}   (paper: ~1.66)")
    print(f"  n_eff(405, th_m) = {n_eff(lam_p, theta_m):.4f} (paper Table I: ~1.660)")
    print(f"  2*Ns (external)  = {2*Ns_ext:.3e} /s   (paper: ~4.2e7)")
    print(f"  2*Ns (internal)  = {2*Ns_int:.3e} /s   (should match external)")
    print(f"  eta = 2Ns/Np     = {2*Ns_ext/Np:.3e}        (paper: ~2.6-2.7e-10)")
    return 2 * Ns_ext / Np


if __name__ == "__main__":
    from textlog import tee_stdout
    with tee_stdout("sanity_check.txt"):
        sanity_check()
