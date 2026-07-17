# Pump focusing vs. coincidence rate — model gap and MMF collection

**Status:** findings memo, 2026-07-17. Deferred future work; not yet implemented.
See `CLAUDE.md` and the memory note `mmf-collection-focusing-model`.

## Question

What changes in the computation if the pump is focused onto the crystal
(collimated 2 mm dia → tighter spot), and what does it do to the single and
coincidence counts for the user's own setup (`detector_counts.py` /
`coincidence_counts.py`)?

## What the current model says

`W` (pump 1/e² radius) enters `spdc_physics.Ns` only inside the transverse
phase-matching Gaussian `exp(−½(ξs−ξi)²)`, ξ = κ·W — **not** the prefactor.
Consequences, confirmed by running the code with `W` overridden:

- **Total / singles flux is W-independent** (total 2·Ns = 4.44e7 /s at any W;
  detector singles 1235 /s at both W = 1 mm and 0.1 mm). This is the standard
  result: pair rate is set by pump *power*, not spot size, in the loose-focusing
  regime.
- **Coincidences fall** with focusing: heralding 0.156 → 0.129 (≈ −17% for
  W 1 mm → 0.1 mm), because the azimuthal correlation width
  σφ = (1/W)/(k₀ sinθ) grows 10× (0.14° → 1.4°) and the idler cone spills over
  the detector's far-field aperture.
- Focusing stays within the model's validity (Rayleigh range z_R ≈ 78 mm at
  0.1 mm radius ≫ the 0.1 mm crystal), so no converging-beam correction needed.

## What the experiment shows (the model is wrong)

Measured rates (per second; dark counts 400/800 subtracted):

| | Collimated | Focused | Δ |
|---|---|---|---|
| Singles A | 950 | 1100 | +16% |
| Singles B | 1150 | 1300 | +13% |
| Coincidences | 95 | 125 | +32% |
| Heralding C/S_A | 0.100 | 0.114 | +14% |

Focusing **raised** coincidences ~32% and heralding ~14% — opposite sign to the
model. Both singles and coincidences up, coincidences ~2× more: the fingerprint
of improved **pump–collection spatial-mode overlap**, not far-field aperture
geometry.

## Root cause

The current coincidence model treats collection as a **bare 6 mm far-field
angular bucket** at 625 mm with an implicit **plane-wave pump**. In that picture
singles can't depend on W and coincidences can only *lose* from correlation-width
spill — so it structurally cannot produce a focusing gain.

The real collection is **multimode fiber**: Thorlabs MR17L01 (50 µm core,
0.22 NA) + CFC11P-B collimator (f = 11 mm, output waist 2.38 mm), crystal at
625 mm. That is a **finite-étendue spatial + angular mode**, sensitive to
mode-matching the momentum-only model ignores.

### Geometry (from datasheets)

- **Collection footprint at the crystal ≈ 2.38 mm dia** (collimated mode,
  z_R ≈ 5.5 m ≫ 0.625 m so it barely diverges) — ~matched to the 2 mm pump, so
  the setup already sits near the mode-matching regime.
- **Angular acceptance** = core/f = 50 µm / 11 mm = **±0.13°** about the 3° axis
  (a few hundred collected modes).
- **Pump divergence** at 0.2 mm dia = ±λ/πW ≈ **±0.074°** — inside ±0.13°, so
  focusing concentrates pairs spatially with no angular penalty yet ⇒ singles up
  a little, coincidences (both twins) up more, heralding up. Matches the data.

### Prediction

An **optimum near ~0.1 mm pump diameter** (pump divergence ≈ collection
acceptance). Below that the idler angle spills past the fiber acceptance and
coincidences should start to fall again — a testable prediction.

## To implement (future)

Replace the far-field hard-aperture kernel in `coincidence_counts.py` with a
**finite-étendue collection mode** per arm:

1. Pump as a real Gaussian of waist W (add its transverse profile — currently
   plane-wave).
2. Each arm a collection mode: footprint ≈ 1.19 mm waist + angular acceptance
   ±0.13°, few-hundred modes (MMF étendue).
3. Singles and coincidences as mode-overlap integrals that depend on W.
4. Validate against the two data points above; then scan W to locate the
   optimal spot size.

**Still needed:** the actual focused pump waist at the crystal, and the pump
beam quality (M² / divergence).

**References:** Bennink, *Optimal collinear Gaussian beams for SPDC*, PRA 81,
053805 (2010); Ling, Lamas-Linares, Kurtsiefer, PRA 77, 043834 (2008).
