# Momentum Reconstruction Issues

**File**: `VIKING_FOLDER/DATA_PARSING/KLong_save_vectors.C`  
**Acceptance mirror**: `VIKING_FOLDER/DATA_PARSING/KLong_save_momentum_acceptance.C`  
**Date diagnosed**: March 2026  
**Status**: Fixes applied to both macros ✅

---

## Table of Contents

1. [Background: The Reconstruction Chain](#1-background-the-reconstruction-chain)
2. [Issue 1: Systematic Momentum Underestimation](#2-issue-1-systematic-momentum-underestimation)
   - 2.1 Root Cause: Time–Path Inconsistency
   - 2.2 Jensen Bias from the TOF Y Smear
   - 2.3 Relativistic Amplification at High Momentum
   - 2.4 Near-Parallel Track PoCA Noise
   - 2.5 Proposed Fix
3. [Issue 2: Resolution Independent of PIZZA–TOF Gap](#3-issue-2-resolution-independent-of-pizzatof-gap)
   - 3.1 Expected Behaviour
   - 3.2 Root Cause: PIZZA Position Smear Dominates L_PV Noise
   - 3.3 Configuration-Dependence of Noise Terms (Quantified)
   - 3.4 Why the TOF Baseline Benefit is Masked
   - 3.5 Secondary Effect: Tracking Lever Arm
   - 3.6 Proposed Fix
4. [Summary of Changes Required](#4-summary-of-changes-required)
5. [Validation Tests](#5-validation-tests)

---

## 1. Background: The Reconstruction Chain

The full K_L momentum reconstruction proceeds in the following steps:

```
1. Track fit (TGraph)
       ↓  fit_pip, fit_pim  (TrackFit: orig + dir)
2. PoCA vertex
       ↓  decay_vertex  (TVector3)
3. Pion path: L_track = |pip_tof_pos − pip_pizza_pos|
              L_PV    = |pip_pizza_pos − decay_vertex|
4. Pion TOF velocity: v_π = L_track / (tof_time − pizza_time)
       ↓  β_π, γ_π
5. Pion momentum: p_π = γ_π m_π v_π
6. Kaon decay time: t_decay = t_vertex + L_PV / (β_π c)
       (propagating pion backwards from PIZZA to vertex)
7. Kaon path: L_K = |decay_vertex − production_vertex|
              t_flight = t_decay − t_production
8. Kaon velocity: v_K = L_K / t_flight
       ↓  β_K, γ_K
9. Kaon momentum: p_K = γ_K m_K v_K
```

Steps 3–6 depend critically on the consistency of positions and times used. The two issues below arise from inconsistencies introduced by the original use of smeared/fabricated positions.

---

## 2. Issue 1: Systematic Momentum Underestimation

**Symptom**: `p_reco < p_true` across all configurations, with the bias growing worse at high momentum (> 1 GeV/c). The bias is not a constant offset — it increases roughly as `p²` or faster.

### 2.1 Root Cause: Time–Path Inconsistency

In the original code, `pip_tof_time` is the **real Geant4 hit time** at the true pion hit position on the TOF bar. However, `pip_tof_pos` was constructed as:

```cpp
// ORIGINAL (incorrect):
TVector3 pip_tof_pos = bar_centre;
pip_tof_pos.SetY(pip_tof_pos.Y() + gRandom->Gaus(0, 10.0));  // ±10 cm Y smear
```

This gives the **bar centre** (a fixed geometry point) plus a random Y offset — it is **not the same point** at which `pip_tof_time` was recorded. The velocity computed from this pair is therefore:

```
v_fabricated = L_fabricated / t_real  ≠  v_true
```

Because `L_fabricated` (bar-centre distance) is typically **shorter** than `L_real` (distance to true hit position, which is displaced in Y from the bar centre by the pion trajectory), the computed velocity is **biased low**, making the pion appear slower. A slower pion implies a shorter path back to the vertex (step 6), which then biases the kaon velocity and momentum downward.

The effect is amplified at high momentum because high-momentum pions are more forward-boosted and hit the TOF bar well away from its centre in Y, maximising `L_fabricated − L_real`.

### 2.2 Jensen Bias from the TOF Y Smear

Even ignoring the time inconsistency, the `Gaus(0, 10 cm)` Y smear introduces a **Jensen bias** in the path length:

```
E[ sqrt(ΔY² + ΔZ²) ]  >  sqrt( E[ΔY²] + ΔZ² )
```

because `sqrt` is a concave function. The expectation of the path is always larger than the path computed from the expectation of ΔY. This means that on average:

```
<L_track>  >  L_track(ΔY=0)
```

producing a systematic overestimate of the pion velocity and, consequently, of the pion momentum. The net effect on kaon momentum is more complex (it feeds through step 6 in both numerator and denominator), but the combined effect of 2.1 and 2.2 yields a net **underestimate** because the time inconsistency (2.1) dominates.

### 2.3 Relativistic Amplification at High Momentum

Both biases are amplified at high momentum through the relativistic factor. The kaon momentum is:

```
p_K = γ_K m_K v_K = m_K β_K / sqrt(1 − β_K²)
```

The fractional momentum error is related to the fractional velocity error by:

```
δp_K / p_K ≈ γ_K² × δβ_K / β_K
```

At:
- `p_K = 1 GeV/c`: γ_K ≈ 2.0,  amplification ≈ ×4
- `p_K = 2 GeV/c`: γ_K ≈ 4.1,  amplification ≈ ×17
- `p_K = 3 GeV/c`: γ_K ≈ 6.1,  amplification ≈ ×37

This is why even a small fractional bias in `v_K` produces a large fractional bias in `p_K` at high momentum, and why the momentum underestimation grows steeply with `p_true`.

### 2.4 Near-Parallel Track PoCA Noise

The Point of Closest Approach (PoCA) between two tracks is computed using:

```cpp
TVector3 n = v1.Cross(v2);          // normal to both track directions
double dist = (p2 - p1).Dot(n) / n.Mag();
```

The denominator `|v1 × v2|` = `|v1||v2| sin(θ)` where θ is the opening angle between the tracks. For high-momentum K_L decays, the daughter pions are strongly forward-boosted (small opening angles), so:

```
sin(θ) → 0   ⟹   |v1 × v2| → 0
```

The PoCA becomes numerically unstable: any noise in the track directions is divided by a very small number, producing a wildly displaced `decay_vertex`. This bad vertex then corrupts `L_PV` (step 3) and `t_decay` (step 6), cascading into a bad `p_K`.

The effect is event-by-event catastrophic (not a smooth bias), but for a sample dominated by high-momentum decays it appears as a heavy high-side or low-side tail in `Δp`, and shifts the mean.

### 2.5 Proposed Fix

Replace all fabricated/smeared positions with track-extrapolated positions, and reject near-parallel events:

```cpp
// (A) Near-parallel cut — before PoCA
double cos_opening = fit_pip.dir.Dot(fit_pim.dir);
if (cos_opening > 0.9998) continue;   // reject θ < ~1.1 degrees

// (B) Helper lambda — evaluate fitted track at a given z
auto eval_track_at_z = [](const TrackFit& f, double z) -> TVector3 {
    double t_param = (z - f.orig.Z()) / f.dir.Z();
    return f.orig + t_param * f.dir;
};

// (C) PIZZA position from track extrapolation
TVector3 pip_pizza_pos = eval_track_at_z(fit_pip, pip_pizza_z);
TVector3 pim_pizza_pos = eval_track_at_z(fit_pim, pim_pizza_z);

// (D) TOF position from track extrapolation
TVector3 pip_tof_pos = eval_track_at_z(fit_pip, z_tof);
TVector3 pim_tof_pos = eval_track_at_z(fit_pim, z_tof);
```

This ensures that `pip_tof_pos` and `pip_tof_time` correspond to the **same physical point** (the track crossing of the TOF plane at `z_tof`), eliminating the Jensen bias and the time–path inconsistency simultaneously.

---

## 3. Issue 2: Resolution Independent of PIZZA–TOF Gap

**Symptom**: The momentum resolution (FWHM of `Δp = p_reco − p_true` per momentum bin) is nearly identical across all four detector configurations (T1-240 through T1-480), even though the PIZZA–TOF baseline varies by a factor of ~2 (245 cm → 485 cm). This makes it impossible to assess the benefit of the longer TOF arm.

### 3.1 Expected Behaviour

The fractional time-of-flight resolution scales as:

```
σ(p)/p  ∝  σ(t) / L_TOF
```

where `L_TOF` is the pion flight path from PIZZA to TOF. Since Config A (T1-240) has `L_TOF ~ 485 cm` and Config D (T1-480) has `L_TOF ~ 245 cm`, we expect Config A to have roughly **twice the momentum resolution** of Config D — a clearly visible separation in the FWHM plots.

### 3.2 Root Cause: PIZZA Position Smear Dominates L_PV Noise

In the original code:

```cpp
double smear_sigma = 5.0;   // cm — PIZZA position smear
// ...
TVector3 pip_pizza_pos;
pip_pizza_pos.SetX(gRandom->Gaus(pizza_x, smear_sigma));
pip_pizza_pos.SetY(pip_pizza_pos.Y());
pip_pizza_pos.SetZ(pizza_z);
```

The 5 cm Gaussian smear on the PIZZA x-position enters `L_PV = |pip_pizza_pos − decay_vertex|`. For a typical decay, `L_PV ∼ 50–200 cm`, so:

```
σ(L_PV) / L_PV  ≈  5 cm / 100 cm  =  5%
```

This 5% fractional noise on `L_PV` is **identical for all four configurations** because `L_PV` depends on the decay vertex position (driven by the K_L kinematics and decay physics), not on where the TOF wall is placed. The smear therefore dominates the momentum resolution and masks the factor-of-2 variation from the TOF baseline.

### 3.3 Configuration-Dependence of Noise Terms (Quantified)

| Noise source | Config A (T1-240) | Config D (T1-480) | Config-dependent? |
|---|---|---|---|
| TOF time resolution (smear_time = 200 ps) | σ(p)/p ∝ 1/L_TOF ≈ 1/485 cm | σ(p)/p ∝ 1/245 cm | **YES** — factor 2 |
| PIZZA x smear (5 cm) → L_PV noise | σ(L_PV)/L_PV ≈ 5% | σ(L_PV)/L_PV ≈ 5% | **NO** — identical |
| Track direction noise → L_track | Small (many tracker hits) | Small | Weakly config-dep. |
| PoCA vertex noise | Event-dependent (θ) | Event-dependent (θ) | Weakly config-dep. |

The PIZZA smear (row 2) overwhelms the TOF baseline effect (row 1) because:

```
σ(p)_PIZZA / σ(p)_TOF  ≈  (5% × p × γ²)  /  (σ_t c / L_TOF × p × γ²)
                        =  5 cm × L_TOF / (c × σ_t)
                        ≈  5 cm × 485 cm / (30 cm/ns × 0.2 ns)
                        ≈  5 / 6  ≈ 0.83   [Config A]
                        ≈  5 cm × 245 cm / (30 cm/ns × 0.2 ns)
                        ≈  5 / 12 ≈ 0.42   [Config D]
```

So the PIZZA smear contributes ~83% of the TOF-baseline noise in Config A and ~42% in Config D. The absolute noise from PIZZA smear is the same in both — it simply represents a larger *fraction* of Config D's TOF noise, pulling the two configs' total noise closer together.

### 3.4 Why the TOF Baseline Benefit is Masked

Combining in quadrature (schematically):

```
σ_total² = σ_TOF_baseline² + σ_PIZZA_smear² + (other terms)

Config A:  σ_total² = σ_TOF_A² + (5cm)²
Config D:  σ_total² = σ_TOF_D² + (5cm)²  where σ_TOF_D ≈ 2×σ_TOF_A
```

If `σ_TOF_A ≈ 4 cm` equivalent and `σ_TOF_D ≈ 8 cm` equivalent, the 5 cm PIZZA smear pushes both totals close together:

```
Config A:  sqrt(4² + 5²) = 6.4 cm  →  only 20% better than Config D
Config D:  sqrt(8² + 5²) = 9.4 cm  (not 2× worse as expected without smear)
```

Setting `smear_sigma = 0` restores the factor-of-2 separation between Config A and Config D.

### 3.5 Secondary Effect: Tracking Lever Arm

The PIZZA position smear also affects the tracking reconstruction indirectly. The PIZZA detector sits at a different z for each configuration, providing a different lever arm for the TGraph fit. With a 5 cm artificial smear applied uniformly regardless of the PIZZA z position, the smear partially negates the geometric advantage of larger tracker separation (configs T1-240 vs T1-480 differ in tracker arm by ~200 cm). This further compresses the apparent difference between configurations.

### 3.6 Proposed Fix

```cpp
double smear_sigma = 0.0;   // cm — use track extrapolation instead (see eval_track_at_z)
```

With the PIZZA position taken from track extrapolation (fix (C) in §2.5 above), the `smear_sigma` variable becomes unused and should be set to zero to avoid confusion.

---

## 4. Summary of Changes Required

| # | Location | Change | Reason |
|---|---|---|---|
| 1 | Both macros, `smear_sigma` declaration | `5.0` → `0.0` cm | Remove artificial PIZZA smear |
| 2 | Both macros, before PoCA calculation | Add `cos_opening > 0.9998` cut | Reject numerically unstable near-parallel events |
| 3 | Both macros, before PIZZA position | Add `eval_track_at_z` lambda | Helper for consistent track extrapolation |
| 4 | Both macros, PIZZA position computation | Replace `Gaus(pizza_x, smear)` with `eval_track_at_z(fit_pip, pip_pizza_z)` | Consistent position for L_PV |
| 5 | Both macros, TOF position computation | Replace `bar_centre + Gaus(0,10)` with `eval_track_at_z(fit_pip, z_tof)` | Make TOF position consistent with TOF time |

All 5 changes have been applied to both `KLong_save_vectors.C` and `KLong_save_momentum_acceptance.C`. ✅

### Code applied (both macros)

```cpp
// Change 1: smear disabled (line ~448 in vectors, ~358 in acceptance)
double smear_sigma = 0.0;   // cm — PIZZA position smear: set 0 (track extrapolation used instead)

// Change 2: near-parallel cut (before PoCA)
double cos_opening = fit_pip.dir.Dot(fit_pim.dir);
if (cos_opening > 0.9998) continue;   // cos(1.1 deg) ~ 0.9998
// (use `break` instead of `continue` in KLong_save_momentum_acceptance.C
//  because the logic is inside a do { } while(false) block)

// Change 3: track extrapolation helper
auto eval_track_at_z = [](const TrackFit& f, double z) -> TVector3 {
    double t_param = (z - f.orig.Z()) / f.dir.Z();
    return f.orig + t_param * f.dir;
};

// Change 4: PIZZA position from track
TVector3 pip_pizza_pos = eval_track_at_z(fit_pip, pip_pizza_z);
TVector3 pim_pizza_pos = eval_track_at_z(fit_pim, pim_pizza_z);

// Change 5: TOF position from track (consistent with tof_time which is also at this z plane)
TVector3 pip_tof_pos = eval_track_at_z(fit_pip, z_tof);
TVector3 pim_tof_pos = eval_track_at_z(fit_pim, z_tof);
// Both path and time now correspond to the same track-extrapolated point
```

---

## 5. Validation Tests

After re-running the reconstruction with the fixes applied, perform the following checks:

### Test 1: Momentum bias removal

Plot `p_reco / p_true` vs `p_true` (truth-matched events, all configs). 

- **Before fix**: ratio drops below 1.0 with a slope, reaching ~0.7–0.8 at 3 GeV/c
- **After fix**: ratio should be flat at ~1.0 across the full momentum range
- If still biased: check that `pip_tof_time` unit is seconds and `pip_tof_pos` is in cm (velocity computed as cm/s before dividing by c)

### Test 2: Configuration separation in FWHM

Run `KLong_plot_gaussian_spreads_truncated.C` with all four configurations.

- **Before fix**: four FWHM curves nearly coincide
- **After fix**: curves should separate with Config A (T1-240, longest TOF arm) giving the best (smallest) FWHM, and Config D (T1-480, shortest TOF arm) giving the worst
- The separation should grow with momentum (because of the γ² amplification)

### Test 3: Δp centred at zero

Plot the `Δp = p_reco − p_true` distribution integrated over all momenta.

- **Before fix**: distribution mean < 0 (systematic underestimate)
- **After fix**: distribution should be centred at 0 within statistical uncertainty
- A residual small offset (< 1%) is acceptable and may arise from the near-parallel cut changing the accepted phase space

### Test 4: Ideal baseline (optional)

Temporarily set `smear_time_sigma = 0` (if such a parameter exists) and rerun to obtain the ideal (time-resolution-only-limited) momentum resolution. This provides a lower bound against which detector timing improvements can be benchmarked.

---

*Document created March 2026. Fixes applied to both macros March 17, 2026.*
