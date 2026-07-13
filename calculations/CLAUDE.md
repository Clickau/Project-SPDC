# SPDC-in-BBO flux calculations

Numerical recreation of the angle-resolved parametric-fluorescence flux from
**Hsu & Lai, "Angle-Resolved Spectroscopy of Parametric Fluorescence"
(arXiv:1302.6226)**. Type-I (e → o + o) SPDC in BBO, pump λp = 405 nm,
degenerate signal/idler at 810 nm. The user is a physicist; treat the physics
as the source of truth and keep the scripts numerically consistent with each
other (they share `spdc_physics.py` — **edit physics there, nowhere else**).

`paper.pdf` is the reference paper (equations verified against the rendered PDF,
not just text extraction). `code.txt` is the user's original Mathematica
implementation of Eq.(9) that was debugged/rewritten in Python (see bug list
below); keep it only for reference. Both live at the repo root (reference
material, not source or generated output).

## Layout (post 2026-07-05 reorg: `src/` + `output/` folders)

All Python/LaTeX source lives in `src/`; all generated PNGs and the LaTeX log
live in `output/`. Every plot script resolves its output path off `__file__`
(`OUTDIR = Path(__file__).resolve().parent.parent / "output"`), so they write
to `output/` correctly whether run as `python src/plot_fig6.py` from the repo
root or `python plot_fig6.py` from inside `src/`. Run
`python src/spdc_physics.py` to execute `sanity_check()` against the paper.

## Files

- **`src/spdc_physics.py`** — the physics library; all constants and physics
  live here, the plot scripts contain none. Defines the SI constants, experiment
  parameters, BBO Sellmeier (`n_o`, `n_e_principal`, `n_eff`), `sinc`, `Np`,
  the unified flux `Ns(lam_s, theta_s, ...)` (numpy-broadcasts λs against θs,
  so one function serves both the Fig. 6 maps and the detector scans; keyword
  options: `external` angle convention, `Nxi` quadrature points, `lam_i_max`
  idler cutoff, `idler_theta_window` — restricts the ξi quadrature to idlers
  exiting into a given external polar-angle window, turning `Ns` into the
  joint signal-AND-idler coincidence kernel; default None is the unchanged
  singles marginal), `tuning_curve`, and `crystal_rotation`. Running
  `python src/spdc_physics.py` executes `sanity_check()` against the paper.
  `Ns` also takes per-call overrides of the experiment parameters
  (`lam_p, Pp, L, W, deff`), and `tuning_curve`/`crystal_rotation` of `lam_p`;
  all default to the paper constants, so existing callers are bit-identical
  (sanity check re-verified after adding them).
- **`src/detector_scan_common.py`** — shared machinery for the two
  detector-scan diagrams: `LAM_I_MAX` (3500 nm idler cutoff),
  `spectra_vs_thetam` (the θm × λ count matrix every panel derives from),
  `rotation_axes` (α axis + θm secondary-axis interpolators),
  `two_order_crossings`, `summarize`, `draw_panel`.
- **`src/plot_fig6.py`** — Fig. 6 recreation. Writes `output/spdc_fig6.png`:
  the four stacked log-scale λs–θs maps Fig. 6(e)–(h) with white-dashed
  perfect-phase-matching tuning curves. Also runs the sanity check first.
- **`src/plot_detector_scan.py`** — fixed detector at external angle θs = 3°,
  rotate the crystal, plot wavelength-integrated count vs crystal rotation
  angle α. Two panels: broadband + hard 805–815 nm top-hat. Writes
  `output/detector_scan.png`.
- **`src/plot_detector_scan_bandpass.py`** — same scan with reworked
  wavelength selection, three panels: (1) broadband (all λ), (2) Gaussian
  bandpass FWHM 10 nm @ 810 nm, (3) monochromatic spectral density dN/dλ at
  810 nm (counts/s/nm). Writes `output/detector_scan_bandpass.png`.
- **`src/detector_counts.py`** — photons/s on a real finite-aperture detector
  for the USER'S OWN setup (not the paper's), all knobs in a PARAMETERS block
  at the top (incl. the detector lab angle): pump 405 nm / 40 mW, collimated
  W = 1 mm (2 mm beam dia.), BBO L = 0.1 mm, detector θs = 3° at 625 mm with
  6 mm circular opening, Gaussian bandpass 10 nm FWHM @ 810 nm. Solves the
  ideal cut angle (degenerate cone through the detector at normal incidence,
  θm = 29.24° for 3°), models the circular aperture exactly (polar span +
  azimuthal arc `2·arccos[(ρ²+R0²−r²)/(2ρR0)]` of the ring), prints the rate
  (1.24e3 /s detected per detector; full ring 5.40e4 /s; the loss-free values
  were 6.46e3 / 2.83e5), and writes `output/detector_scan_mysetup.png` (rate
  vs crystal rotation α — a broad symmetric thin-crystal sinc² pattern,
  2 orders down at α ≈ ±3.1°).
  EFFICIENCIES block (added 2026-07-05, datasheet values same day): pump
  chain `ETA_PUMP = ETA_MIRROR_F01^N_MIRROR_F01 · ETA_HWP · T_CRYSTAL_IN`
  = 0.92²·0.999·1 = 0.846 — PP is the LASER OUTPUT power and `PHYS_KW` feeds
  `Pp = PP·ETA_PUMP` (33.8 mW) to `Ns`; collection chain `ETA_COLLECT =
  ETA_COLLECT_BASE · T_FILTER_PEAK` = (1·0.958·0.80·1·0.30)·0.983 = 0.226,
  with `ETA_COLLECT_BASE = T_CRYSTAL_OUT · ETA_MIRROR_P01 ·
  ETA_FIBER_COUPLING · ETA_FIBER · ETA_APD`, multiplies `detected_rate`'s
  outputs (the Gaussian in `filter_grid` stays peak-normalised;
  `T_FILTER_PEAK` = 0.983 carries the peak). The 0.80 fiber coupling is
  all-in (AR-coated lens assumed; 15% lens→core, 3.5% Fresnel at input,
  negligible fiber propagation, 5% detector-box mating; 0.85·0.965·0.95 =
  0.78 ≈ 0.80), so `ETA_FIBER` stays 1.0. Crystal AR-coated both sides ⟹
  both `T_CRYSTAL_*` = 1. APD 0.30 @ 810 nm. Verified: all-1.0 reproduced
  the loss-free numbers exactly, and a 0.5 collection × 0.8 pump override
  scaled the rate by exactly 0.4.
- **`src/coincidence_counts.py`** — signal–idler coincidences between TWO
  detectors at ±3° for the user's setup (reuses `detector_counts`'
  parameters and `ideal_cut_deg`; that script is untouched). Per-arm
  PARAMETERS (angle, distance, aperture, filter; `filter_fwhm=None` = open
  arm). Coincidence = the `Ns` integral with (a) `idler_theta_window` = arm
  B's polar span, (b) λ weight `T_A(λs)·T_B(λi(λs))`, (c) the singles' arc
  factor dφ_A/2π replaced by the A-arc × mirrored-B-arc overlap smeared by
  the azimuthal momentum-correlation width σφ = (1/W)/(k0·sinθ) ≈ 0.14° (vs
  ~10.5° full arc). With the real efficiencies prints singles 1.24e3 /s per
  arm, coincidences 192 /s, heralding 0.156 = η_B·0.689 (loss-free values:
  6.46e3, 4.45e3, 0.689 ≈ (1/√2)·(geometry ≈ 0.97)), accidentals ~1.5e−3 /s
  @ τ = 1 ns, and a built-in validation (open B ⟹ recovers the singles,
  rel. diff 0).
  Writes `output/coincidence_scan.png` (singles + coincidences + heralding vs
  α; heralding is flat ≈ 0.69 across the whole scan — near-degenerate filters
  keep the twins mirrored, so only the filters set Rc/Rs, not geometry).
  Component efficiencies come from `detector_counts`' EFFICIENCIES block
  (ONE shared set, per the user — both arms have identical chains):
  `collection_eta(arm)` = `ETA_COLLECT_BASE` × (`T_FILTER_PEAK` only if the
  arm has a filter, so the open-B validation stays exact); singles ×η_arm,
  coincidences ×η_A·η_B ⟹ heralding Rc/Rs_A picks up a factor η_B; pump
  chain enters via the shared `PHYS_KW`. Accidentals use the η-scaled
  singles.
- **`src/spdc_eq9_note.tex`** (built log at `output/spdc_eq9_note.log`) —
  write-up of the paper's angle conventions and errata: Eq.(9) is correct as
  printed (external angle, per-unit-azimuth); the genuine misprints are in
  Eq.(7) (see below). Rewritten 2026-07-13, retracting the earlier
  "2π·ns² typo" claim.

(History: the pre-refactor layout had the physics inside `spdc_fig6.py` and a
duplicated flux kernel in each detector script. The unified `Ns` was verified
**bitwise identical** to both old kernels — and `tuning_curve` /
`crystal_rotation` / the sanity totals likewise — before the old files were
deleted. On 2026-07-05 the flat file layout was split into `src/` + `output/`;
all scripts re-verified to reproduce the same numbers after the move.)

## Master equation — Eq.(9), angle-resolved signal flux

As implemented in `spdc_physics.py` (`Ns` evaluates the paper's Appendix
intermediate equation and integrates the signal azimuth over 2π; in the
paper's variables that is exactly **2π × the printed Eq.(9)**):

```
Ns(λs,θs) = 2π · [ 2π·√(2π) · ħ·deff²·ωp·L²·Np / (ε0·ns·ni·np·λs⁵·λi)
            · sin(2θs) · dλs · dθs
            · ∫dξi  exp(−½(ξs−ξi)²) · sinc²(½·L·Δkz) ]
```

- `sinc(x) = sin(x)/x`, so `sinc²` is squared (NumPy's `np.sinc` is
  `sin(πx)/(πx)` — divide the argument by π).
- θs is the **external (lab) angle** (paper Fig. 2: θ′s internal, θs external);
  κs = k0·sin(θs) with k0 = 2π/λs (κ is conserved across the exit face, so
  this equals ks·sinθ′s with ks = ns·k0).
- ξs = κs·W, ξi = κi·W (dimensionless transverse momenta);
  W = √(60 µm × 30 µm) ≈ 42.4 µm.
- Longitudinal mismatch (Appendix form):
  `Δkz = √(ks² − κs²) + √(ki² − (ξi/W)²) − kp`, with `k = 2π·n/λ` (in-medium)
  for each wave.
- Energy conservation → idler: `λi = λs·λp/(λs − λp)`.

### Paper errata status — re-checked 2026-07-13 (derivation in `spdc_eq9_note.tex`)

**The printed Eq.(9) is CORRECT** — an earlier claim here that it was
"~2π·ns² ≈ 17× too low" is **retracted**. That claim misread θs as the
internal angle; Fig. 2 defines θs as external, and the unprimed ks in the
appendix's "κs = ks·sinθs" is the vacuum wavenumber (the paper primes
in-medium quantities, `ks' = ns·ωs/c`). With κs = k0·sinθs the Jacobian
`κs·dκs = ½·k0²·sin2θs·dθs` carries no ns², and the printed prefactor
`2π·√2π` is exact — as a density **per unit azimuth** (Eq.(9) merely drops
the dφ symbol; Sec. III.D integrates φ = 0→2π explicitly, which is the outer
2π in the code's normalisation above). The genuine misprints are confined to
the main-text Eq.(7):

1. **Factor 4:** Eq.(7)'s prefactor is `2π⁴` in the main text but `8π⁴` in the
   Appendix. `8π⁴` is correct — it reproduces Table I's simulated η (2.6e-10;
   code gives 2.72e-10), the independent closed-form Eq.(8) from ref. [18]
   (Koch et al. 1995; evaluates to 2.63e-10), and the measured 2.7–2.8e-10.
   The relabel ∫d²ξ→∫d²ξi is a unit-Jacobian shift and cannot produce the 4.
2. **Sign:** Eq.(7)'s exponent prints `exp(+½ξ²)`; must be `−½ξ²`, as the
   Appendix has it.

## Indices, Sellmeier, constants

- Signal & idler ordinary: `ns = no(λs)`, `ni = no(λi)`. Pump effective
  extraordinary (Eq. 4): `np = n_eff(θm,λp) = [cos²θm/no² + sin²θm/ne²]^(−½)`.
- BBO Sellmeier (λ in µm):
  - `no² = 2.7359 + 0.01878/(λ²−0.01822) − 0.01354·λ²`
  - `ne² = 2.3753 + 0.01224/(λ²−0.01667) − 0.01516·λ²`
- Constants/params: ħ = 1.054571817e−34 J·s; ε0 = 8.8541878188e−12 F/m;
  c = 2.99792458e8 m/s; deff = 1.75 pm/V; L = 3 mm; λp = 405 nm; Pp = 80 mW ⟹
  Np = Pp/(ħ·ωp) = 1.63e17 /s.
- **Angle convention:** the paper's Eq.(9) and Fig. 6 both use the *external*
  (air/lab) angle (κs = k0·sinθ_ext, k0 = 2π/λs; Fig. 2 convention) —
  `external=True`, the default. `external=False` parametrises the same κs by
  the *internal* angle (κs = ks·sinθ, ks = ns·2π/λs), an equivalent change of
  variables. Both give identical total flux — only the y-axis scale of a
  per-angle plot differs by ~ns. Use external for figures.

## Sanity checks (must hold)

- η = 2Ns/Np ≈ 2.6–2.7e−10 at 810 nm for θm ≈ 28.8–29.4°.
- Total degenerate flux 2Ns ≈ 4.2e7 /s at Pp = 80 mW (code gives 4.44e7).
- Fig. 6 phase-matching angles θm = 28.6°, 28.8°, 29.1°, 29.4°;
  λs ≈ 430–1000 nm, θs ≈ −6°…+6°, color on log scale.
- Note: collinear-degenerate (810 nm at θ=0) lands at θm ≈ 28.82° with this
  Sellmeier vs the paper's 28.6° (paper used a NIST index program) — a ~0.2°
  label offset; the physics is correct.

## Detector-scan analysis — key physical findings

Detector fixed at external θs = 3°, crystal rotated. Because the pump refracts
at the entrance face, the external crystal rotation α ≠ internal Δθm:
`sin α = n_eff(λp,θm)·sin(θm − θ_cut)` — a given Δθm needs ~n_eff ≈ 1.66× larger
external α. Plot x-axis = α, secondary top axis = θm.

- **Cone-axis tilt cancels (verified by 2D ray-trace):** rotating the crystal
  tilts the internal pump, but pump-in + signal-out refraction gives an internal
  signal–pump angle ≈ 1.80° for *all* α (1.79–1.81° over α = +0.4…−8.3°),
  because n_o(signal) ≈ n_eff(pump) ≈ 1.66 in BBO. So the simple
  κs = k0·sin(3°) model is geometrically correct.
- **Idler cap:** wavelength integration zeroes signals whose idler exceeds
  ~3500 nm (BBO transparency / index-data limit). Without it the short-λ tail is
  propped up by pairs that cannot be generated (e.g. 433 nm → 6370 nm idler).
- **Broadband asymmetry is real (broadband vs narrowband):** broadband peak at
  α ≈ +0.37° (θm ≈ 29.22°); dropping 2 orders takes only +0.8° on the high side
  (sharp degenerate-caustic edge) but −7° on the low side (the tuning-curve arms
  keep crossing 3°). A 10 nm bandpass @810 nm makes it **symmetric** (±0.3° for
  2 orders) — that is the "intuitive"/filtered case (paper Fig. 3).
- The detector-acceptance factor (dθ, azimuth) only rescales the y-axis; it does
  not change the shape or the orders-of-magnitude conclusions.
- **Coincidences ≠ singles:** `Ns` is a marginal (its ∫dξi integrates the idler
  over everything); the integrand `exp(−½(ξs−ξi)²)·sinc²` IS the joint
  transverse-momentum density, so coincidences = same integral with the idler
  restricted to detector B's acceptance. For the ±3° twin-detector setup the
  geometry is near-lossless (polar blur 1/(W·k0) ≈ 0.007°, azimuthal 0.14°,
  λ-walk-off ≤ 0.04° — all ≪ the 0.275°/±5.2° aperture) and the heralding
  ratio is set by the idler-arm filter at the conjugate wavelength: with
  identical 10 nm Gaussians in both arms, T(λs)·T(1620 nm−λs) = T², giving
  Rc/Rs → 1/√2 ≈ 0.71 (model gives 0.689 incl. geometry). Component/detector
  losses are now modeled (EFFICIENCIES block, 2026-07-05): they multiply the
  heralding by η_B = 0.226, giving the observed-rate prediction 0.156.
