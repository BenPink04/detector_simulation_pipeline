# Tracking Overhaul: Module-Based Hit Reconstruction for Decay Vertex Reconstruction

**Status**: Implementation complete — TGraph approach is the current default; DSSSD pixelation preserved as `_pixelation.C` backup files  
**Affects**: `KLong_save_vectors.C`, `KLong_save_momentum_acceptance.C`, `KLong_save_vectors_pixelation.C`, `KLong_save_momentum_acceptance_pixelation.C`  
**Optional**: `JLabKSteppingAction.cc` (see Implementation Strategy below)

---

## Table of Contents

1. [Overview and Motivation]
(#1-overview-and-motivation)
2. [Implementation Strategy: Where to Apply the Mapping](#2-implementation-strategy-where-to-apply-the-mapping)
3. [Complete Geometry Reference]
(#3-complete-geometry-reference)
4. [Step-by-Step Implementation Plan]
(#4-step-by-step-implementation-plan)
   - [Step 1: Add Device ID Classification Functions]
   (#step-1-add-device-id-classification-functions)
   - [Step 2: Build the Module Geometry Lookup Tables]
   (#step-2-build-the-module-geometry-lookup-tables)
   - [Step 3: Module-Centre Reconstruction Function]
   (#step-3-module-centre-reconstruction-function)
   - [Step 4: Hit Grouping, Intra-Component Intersection, and Coincidence Requirement]
   (#step-4-hit-grouping-intra-component-intersection-and-coincidence-requirement)
   - [Step 5: 2D Overlap Area Construction from Intersection Results]
   (#step-5-2d-overlap-area-construction-from-intersection-results)
   - [Step 6: Line-of-Best-Fit Track Reconstruction]
   (#step-6-line-of-best-fit-track-reconstruction)
   - [Step 7: TOF Wall Treatment] 
   (#step-7-tof-wall-treatment)
   - [Step 8: Mirror Changes to Acceptance Macro](#step-8-mirror-changes-to-acceptance-macro)
5. [Pipeline Integration Notes](#5-pipeline-integration-notes)
6. [Validation Checklist](#6-validation-checklist)
7. [DSSSD Simplification and Arrival-Angle Correction (March 2026)](#7-dsssd-simplification-and-arrival-angle-correction-march-2026)
8. [TGraph-Based Track Fitting — Raw-Position Approach (March 2026)](#8-tgraph-based-track-fitting--raw-position-approach-march-2026)
   - [8.1 Motivation and Differences from the Pixelation Approach](#81-motivation-and-differences-from-the-pixelation-approach)
   - [8.2 The `TrackFit` Struct](#82-the-trackfit-struct)
   - [8.3 `fit_3dgraph_track()`: Implementation Details](#83-fit_3dgraph_track-implementation-details)
   - [8.4 Per-Detector Uncertainty Table](#84-per-detector-uncertainty-table)
   - [8.5 T3/T4 Stereo Uncertainty Propagation](#85-t3t4-stereo-uncertainty-propagation)
   - [8.6 TOF Bar: Geometry Lookup vs Raw Position](#86-tof-bar-geometry-lookup-vs-raw-position)
   - [8.7 Fit Validity Conditions](#87-fit-validity-conditions)
   - [8.8 All Sub-Layers Enabled](#88-all-sub-layers-enabled)
   - [8.9 File Structure After Overhaul](#89-file-structure-after-overhaul)

---

## 1. Overview and Motivation

### Current Behaviour

In `JLabKSteppingAction.cc`, every step in a sensitive volume records `postStepPoint->GetPosition()/cm` — the **true Geant4 hit position** inside the physical material. In `KLong_save_vectors.C`, a 5 cm Gaussian smear is then applied independently to X and Y, and the smeared front/back hits define the track direction via PoCA.

**Problems with this approach:**

- A Gaussian smear in both X and Y is not physically motivated. A tracker straw running along Y can only tell you X (the perpendicular coordinate). It gives no information about Y (the coordinate along the tube).
- A scintillator strip running along Y encodes the same constraint. Two perpendicular strips are needed to locate a point in 2D.
- Using the raw first/last hit as track direction is sensitive to noise and does not exploit all hit information.
- The smear magnitude (5 cm σ) is not tied to any module dimension, making it arbitrary.

### Target Behaviour

Replace raw smearing with a **geometry-aware module reconstruction**:

1. Each device ID maps to a known module type, centre position, and orientation.
2. A hit in a module constrains only the coordinate **perpendicular to the module's long axis** (the "measurement axis"). The coordinate along the long axis is unconstrained (smear by the full half-length).
3. When multiple elements within the same detector plane and axis (e.g. adjacent FRI strips from Wall 1, or straws from different sub-layers of the same tracker station) both register a hit, the particle must lie within the **geometric intersection** of their individual constraint regions. This intersection is narrower than any single element's region — it is the tightest X (or Y) interval consistent with all fired elements on that axis.
4. Combining the intersected X constraint with the intersected Y constraint at the same detector plane gives a **2D overlap area** (a rectangle in XY). The centre and half-width of this rectangle, not the centre of a single device, are the position estimate and uncertainty fed into the track fit.
5. Multiple 2D overlap areas from different detector planes → **line-of-best-fit** through the area centres (weighted by 1/σ²) replaces the current front/back endpoint method.
6. The TOF wall is single-layer; each bar gives an X constraint (bar width ± half-width) and an approximate Y position (bar centre ± 10 cm offset to account for light path timing difference along the bar).

---

## 2. Implementation Strategy: Where to Apply the Mapping

There are two places the module-centre mapping could be applied:

### Option A — In the ROOT Analysis Macros Only (RECOMMENDED)

Keep `JLabKSteppingAction.cc` and the ROOT ntuple format **unchanged**. In `KLong_save_vectors.C`, after reading `deviceID` and the raw hit, look up the module geometry using the device ID and **replace** the raw hit position with the module-centre-derived measurement.

**Advantages:**
- No recompilation or rerunning of Geant4 simulations.
- All existing ROOT files remain valid inputs.
- You can iteratively adjust the reconstruction without regenerating data.
- The true hit position is still available for debugging (it has just not been used in reconstruction).

**Disadvantages:**
- Slightly more code in the ROOT macro (geometry lookup tables).

### Option B — In JLabKSteppingAction.cc

Change what is stored in the ntuple from `postStepPoint->GetPosition()/cm` to the module-centre position (derived from `GetCopyNumber(0)`).

**Advantages:**
- Cleaner data: the ROOT file directly contains reconstructed positions.

**Disadvantages:**
- Requires rerunning all Geant4 simulations before any analysis changes can be tested.
- Loses the true hit position permanently.

**Decision: Implement Option A first.** If the reconstruction design is validated and performance is satisfactory, Option B can be considered for production runs.

---

## 3. Complete Geometry Reference

### 3.1 Tracker Straws (Device IDs 1–1952, globally unique)

There are 4 tracker stations, each at a different Z position and orientation. Each station's 488 straws are assigned **globally unique** device IDs encoding both the station number and the straw within the station.

| Station | Globally unique IDs | Z position (cm) | φ rotation | Measurement axis | Long axis |
|---------|--------------------|-----------------|-----------|-----------------|-----------|
| 0 (T1)  | 1–488              | T1 (default 240)| 0°        | **X**           | Y (±45 cm)|
| 1 (T2)  | 489–976            | T2 (default 250)| −90°      | **Y**           | X (±45 cm)|
| 2 (T3)  | 977–1464           | T3 (default 570)| +45°      | **U** (diagonal)| V-axis    |
| 3 (T4)  | 1465–1952          | T4 (default 580)| −45°      | **V** (diagonal)| U-axis    |

> **Note**: T1/T2/T3/T4 positions are overridden by the `detector_simulation_master_FIXED.sh` arguments (the `T1-`, `T2-`, etc. values in the job directory names). Use the actual Z value for reconstruction. If `T1=0` and `T2=0`, those stations are disabled.

**Straw X-centre positions** (same array regardless of φ, taken in station's local frame before rotation):

```
x[i], i=0..121:
-48.4, -47.6, -46.8, ..., +47.6 cm  (step 0.8 cm)
```

Explicit formula: `x[i] = -48.4 + i * 0.8`  (cm), where i = (copyNumber - 1) % 122

**Straw geometry:**
- Outer radius (transverse to tube axis): **0.40 cm** → position uncertainty σ ≈ 0.40/√3 ≈ 0.23 cm (uniform)
- Half-length along tube axis: **45 cm** → position unconstrained over ±45 cm along long axis

**Smearing rule for a straw hit:**
- Measurement axis: use straw centre X (or rotated equivalent), no additional smear needed (or use uniform ±0.4 cm to represent wire resolution)
- Long axis: sample uniformly over ±45 cm (the hit could be anywhere along the tube)

---

**Tracker module structure — 4 straw sub-layers per station:**

Each tracker station (one entry in the T1/T2/T3/T4 list, representing a single measurement orientation) is physically constructed from **4 straw sub-layers** at slightly different Z offsets within the station volume. The sub-layers are staggered in the measurement direction by fractions of the straw pitch (0.8 cm) to ensure a particle always passes through the sensitive region of at least one straw per sub-layer. This is a standard hermeticity design.

The Geant4 simulation assigns a single nominal Z position per station (e.g. T1 = 240 cm) that represents the centre of the 4-sub-layer stack. The 122 copy numbers per station are distributed across the 4 sub-layers combined (the exact assignment must be confirmed from `JLabKDetectorConstruction.cc`).

**Effect on position measurement — intersection approach:**

Because sub-layers are offset, a particle near the boundary between two straws in adjacent sub-layers will fire both straws. The particle's true position is then constrained to the **geometric intersection** of both fired straw constraint regions:

```
Straw A at x_A ± r:   constraint = [x_A − r,  x_A + r]
Straw B at x_B ± r:   constraint = [x_B − r,  x_B + r]
Intersection:                       [max(x_A−r, x_B−r),  min(x_A+r, x_B+r)]
Intersection centre:  0.5 × (lower + upper)
Intersection half-width:  0.5 × (upper − lower)
```

With straw radius r = 0.4 cm and pitch 0.8 cm, two straws in the same sub-layer are exactly touching (no overlap). However, straws from **different sub-layers offset by a half-pitch** (0.4 cm) will produce overlapping constraint regions for particles between them, yielding a non-zero intersection interval and a tighter position estimate.

If the intersection turns out to be empty (the regions are disjoint — two widely separated straws fired, suggesting different particles or delta rays), do not intersect them; treat each straw independently.

> **Confirmed from `JLabKDetectorConstruction.cc` (March 2026):**
> 1. **4 sub-layers per station box.** The straw placement loop runs `for i=0..121` and places four physical straws per i: at sub-layer Z offsets `strawPosZ[4] = {−1.5, −0.5, +0.5, +1.5} cm` relative to the station centre.
> 2. **Z offsets confirmed:** ±0.5 cm and ±1.5 cm from station nominal Z — small compared to the 10 cm gap between neighbouring station boxes. The current implementation treats all straws in a station as being at the single nominal Z, which remains the correct practical choice.
> 3. **Stagger confirmed:** Sub-layers 0 and 2 are at `strawPosX[i]` (no offset). Sub-layers 1 and 3 are at `strawPosX[i] + 0.4 cm` (half-pitch offset). This has been corrected in `tracker_local_x()` in both macros.
> 4. **Straw copy number mapping within one station box (straw_cn = 1–488):**
>    - Sub-layer 0 → straw_cn `i + 1`   → straw_cn 1–122
>    - Sub-layer 1 → straw_cn `i + 123` → straw_cn 123–244 (X positions +0.4 cm)
>    - Sub-layer 2 → straw_cn `i + 245` → straw_cn 245–366 (same X as sub-layer 0)
>    - Sub-layer 3 → straw_cn `i + 367` → straw_cn 367–488 (same X as sub-layer 1)
>
> **deviceID ambiguity — RESOLVED (March 2026):** The original ambiguity (IDs 1–488 appearing at all four Z positions) has been fixed. `JLabKSteppingAction.cc` now encodes both the parent-box copy number (`GetCopyNumber(1)`, values 620–623, uniquely identifying the station) and the straw-within-box copy number (`GetCopyNumber(0)`, 1–488) into a single globally unique ID: `deviceID = (GetCopyNumber(1) − 620) × 488 + GetCopyNumber(0)`. This gives IDs 1–1952, uniquely identifying every straw in the detector. Non-tracker detectors (PIZZA, FRI, TOF) are shifted by +1464 to avoid overlap. All simulations must be **recompiled and re-run** before these new IDs appear in ROOT files.

---

**Converting globally unique tracker deviceID to local X centre:**
```cpp
// Extract straw copy number within station box (1-488), then compute sub-layer:
int straw_cn  = ((deviceID - 1) % 488) + 1;  // 1-488 within station box
int sub_layer = (straw_cn - 1) / 122;        // 0, 1, 2, or 3
int straw_i   = (straw_cn - 1) % 122;        // 0..121
double local_x = -48.4 + straw_i * 0.8;     // cm, base position
// Sub-layers 1 and 3 are staggered +0.4 cm (half-pitch)
if (sub_layer == 1 || sub_layer == 3) local_x += 0.4;

// Station index (determines which Z position and φ rotation to use):
int station = (deviceID - 1) / 488;  // 0=T1/X, 1=T2/Y, 2=T3/U, 3=T4/V
```

**Converting local X to lab-frame position** (for stations with φ ≠ 0):
- Station 1 (φ=0°):  lab_x = local_x,  lab_y = free (smear ±45 cm)
- Station 2 (φ=−90°): lab_y = local_x,  lab_x = free (smear ±45 cm)  
  *(The box is rotated −90° in Z, so local X becomes lab Y)*
- Stations 3 & 4 (φ=±45°): project along rotated axes (u/v stereo):
  - φ=+45°: measures `u = (x+y)/√2`, so the constraint is `x + y = √2 · local_x`
  - φ=−45°: measures `v = (x−y)/√2`, so the constraint is `x − y = √2 · local_x`

> **Simplification for initial implementation**: Use only stations 1 and 2 (X-measuring and Y-measuring) for the two-detector coincidence. The u/v stereo stations can be added later as additional constraints in the line fit.

---

### 3.2 FRI Wall Strips (Device IDs 2001–2052)

Two FRI walls, 26 strips each. Wall 1 (IDs 2001–2026) measures **X** (strips run along Y). Wall 2 (IDs 2027–2052) is rotated 90° and measures **Y** (strips run along X).

**Strip X-centre positions** (in wall's local frame):

```
Wall 1 (i=0..25, IDs 2001-2026): stripPosX[0..25]
Wall 2 (i=26..51, IDs 2027-2052): stripPosX[26..51]
```

The position arrays are identical (wall 2 is a copy placed with 90° rotation):

| Index i (within wall) | Local X centre (cm) | Strip width (cm) | Half-length Y (cm) |
|----------------------|--------------------|-----------------|--------------------|
| 0, 25                | ±66                | 6 (wide)        | 37.25 (short)      |
| 1, 24                | ±60                | 6 (wide)        | 42.75              |
| 2, 23                | ±54                | 6 (wide)        | 53.75              |
| 3, 4, 21, 22         | ±48, ±42           | 6 (wide)        | 59.25              |
| 5, 20                | ±36                | 6 (wide)        | 64.75              |
| 6–9, 16–19           | ±30, ±24, ±18, ±12 | 6 (wide)        | 70.25 (full)       |
| 10, 15               | ±7.5               | 3 (inner)       | 70.25 (no cutout)  |
| 11, 14               | ±4.5               | 3 (inner)       | 70.25 (with cutout)|
| 12, 13               | ±1.5               | 3 (inner)       | 70.25 (with cutout)|

**Converting copy number to local strip centre X:**
```cpp
int fri_index = deviceID - 2001;       // 0..51
int wall      = fri_index / 26;         // 0 = Wall1, 1 = Wall2
int strip_i   = fri_index % 26;         // 0..25 within wall

double stripPosX[26] = {
  -66., -60., -54., -48., -42., -36., -30., -24., -18., -12., -7.5,
  -4.5, -1.5,  1.5,  4.5,  7.5,  12.,  18.,  24.,  30.,  36.,  42.,
   48.,  54.,  60.,  66.
};
double local_x = stripPosX[strip_i];
```

**Lab-frame interpretation:**
- Wall 1 (IDs 2001–2026): rotated 0°. Strips run along Y. Measures **lab X**.  
  `lab_x = local_x`,  `lab_y = free` (±half-length, see table above)
- Wall 2 (IDs 2027–2052): rotated +90° about Z. What was local X becomes lab Y.  
  `lab_y = local_x`,  `lab_x = free` (±half-length)

**FRI wall Z positions** (overridden by script, default):
- Wall 1: F1 (default 260 cm)
- Wall 2: F2 (default 270 cm)

**Smearing rule:**
- Measurement axis: strip centre ± half-width (uniform, ~1.5–3 cm)
- Long axis: free over ± strip half-length (see table)

---

### 3.3 PIZZA Detector (Device IDs 1953–2000)

The PIZZA detector provides a rough angular measure only; it is not used for the primary pion track fit. Its main role is a coincidence trigger to confirm a particle passed through.

- IDs 1953–1976: 24 flat pizza slices (15° each), at P1 position (default 215 cm)
- IDs 1977–2000: 24 conical pizza slices, at P2 position (default 230 cm)
- Outer radius: 37.2 cm, half-thickness: 0.15 cm

> PIZZA is not part of the tracking overhaul. Keep existing treatment.

---

### 3.4 TOF Wall Bars (Device IDs 2053–2070)

Single layer, 18 bars, at E1 position (default 600 cm).

| Bar index i | Device ID | Bar X centre (cm) | Bar half-width X (cm) | Bar type   | Y centre (cm) | Half-height Y (cm) |
|-------------|-----------|------------------|-----------------------|------------|---------------|--------------------|
| 0           | 2053      | −96              | 6                     | Full       | 0             | 130                |
| 1           | 2054      | −84              | 6                     | Full       | 0             | 130                |
| 2           | 2055      | −72              | 6                     | Full       | 0             | 130                |
| 3           | 2056      | −60              | 6                     | Full       | 0             | 130                |
| 4           | 2057      | −48              | 6                     | Full       | 0             | 130                |
| 5           | 2058      | −36              | 6                     | Full       | 0             | 130                |
| 6           | 2059      | −24              | 6                     | Full       | 0             | 130                |
| 7           | 2060      | −12              | 6                     | Full       | 0             | 130                |
| 8           | 2061      | 0                | 6                     | Short lower| −72           | 58.25              |
| 9           | 2062      | 0                | 6                     | Short upper| +72           | 57.25              |
| 10          | 2063      | +12              | 6                     | Full       | 0             | 130                |
| 11          | 2064      | +24              | 6                     | Full       | 0             | 130                |
| 12          | 2065      | +36              | 6                     | Full       | 0             | 130                |
| 13          | 2066      | +48              | 6                     | Full       | 0             | 130                |
| 14          | 2067      | +60              | 6                     | Full       | 0             | 130                |
| 15          | 2068      | +72              | 6                     | Full       | 0             | 130                |
| 16          | 2069      | +84              | 6                     | Full       | 0             | 130                |
| 17          | 2070      | +96              | 6                     | Full       | 0             | 130                |

**TOF bar X positions formula:** `bar_x = -96.0 + i * 12.0` cm, for `i = (deviceID - 2053)`

**Positional uncertainty:**
- X: uniform ±6 cm (bar half-width, no sub-bar X resolution assumed)
- Y: ±10 cm (assumed from light-path timing uncertainty along bar), applied symmetrically about bar Y centre

---

## 4. Step-by-Step Implementation Plan

All changes below are to `KLong_save_vectors.C`. After verifying each step, mirror the same change to `KLong_save_momentum_acceptance.C`.

---

### Step 1: Add Device ID Classification Functions

**Location**: Near the existing `[TO CHANGE]` marker at the top of the macro (around line 188, the `is_tracker` / `is_fri` etc. lambda definitions).

**What to do**: The existing lambdas identify device types. Extend them to know about the two FRI wall layers:

```cpp
// Updated ID ranges (globally unique, March 2026):
auto is_tracker = [](int id){ return id >=    1 && id <=  1952; };
auto is_pizza   = [](int id){ return id >= 1953 && id <=  2000; };
auto is_fri     = [](int id){ return id >= 2001 && id <=  2052; };
auto is_tof     = [](int id){ return id >= 2053 && id <=  2070; };

// Sub-classify FRI into wall 1 (X-measuring) and wall 2 (Y-measuring)
auto is_fri_wall1 = [](int id){ return id >= 2001 && id <= 2026; };  // measures X
auto is_fri_wall2 = [](int id){ return id >= 2027 && id <= 2052; };  // measures Y

// Sub-classify tracker stations
auto tracker_station = [](int id) -> int {
    if(id < 1 || id > 1952) return -1;
    return (id - 1) / 488;  // 0=T1/X, 1=T2/Y, 2=T3/U, 3=T4/V
};
```

**Why**: Downstream logic needs to know whether a hit constrains X or Y.

---

### Step 2: Build the Module Geometry Lookup Tables

**Location**: Add this block at the top of the event-loop analysis section, after opening the ROOT file and before the main loop.

**What to do**: Create lookup arrays/functions mapping device ID → module centre and measurement info.

```cpp
// ======================================================
// MODULE GEOMETRY LOOKUP TABLES
// ======================================================

// --- Tracker ---
// 122 straws per station, straw x in station-local frame
auto tracker_local_x = [](int copyNumber) -> double {
    // straw index within station: 0..121
    int i = (copyNumber - 1) % 122;
    return -48.4 + i * 0.8;  // cm
};
// straw half-length along tube (the unconstrained axis)
const double STRAW_HALF_LENGTH = 45.0;  // cm
// straw half-width (the measurement precision)
const double STRAW_HALF_WIDTH  = 0.40;  // cm (outer radius)

// --- FRI ---
const double FRI_STRIP_X[26] = {
    -66., -60., -54., -48., -42., -36., -30., -24., -18., -12., -7.5,
     -4.5, -1.5,  1.5,  4.5,  7.5,  12.,  18.,  24.,  30.,  36.,  42.,
     48.,  54.,  60.,  66.
};

// FRI strip half-heights (long axis, unconstrained), indexed by strip position 0..25
const double FRI_HALF_LENGTH[26] = {
    37.25, 42.75, 53.75, 59.25, 59.25, 64.75,   // i=0..5
    70.25, 70.25, 70.25, 70.25,                  // i=6..9
    70.25, 70.25, 70.25, 70.25, 70.25, 70.25,    // i=10..15 (inner strips)
    70.25, 70.25, 70.25, 70.25,                  // i=16..19
    64.75, 59.25, 59.25, 53.75, 42.75, 37.25     // i=20..25
};

// FRI strip half-widths (measurement axis precision)
auto fri_half_width = [](int strip_i) -> double {
    // inner strips (i=10..15): 3 cm wide → half-width 1.5 cm
    // outer strips: 6 cm wide → half-width 3.0 cm
    return (strip_i >= 10 && strip_i <= 15) ? 1.5 : 3.0;
};

auto fri_local_x = [&](int deviceID) -> double {
    int fri_index = deviceID - 2001;     // 0..51
    int strip_i = fri_index % 26;
    return FRI_STRIP_X[strip_i];
};

// --- TOF ---
auto tof_bar_x = [](int deviceID) -> double {
    int i = deviceID - 2053;
    return -96.0 + i * 12.0;  // cm
};
const double TOF_BAR_HALF_WIDTH = 6.0;   // cm (X measurement uncertainty)
const double TOF_Y_UNCERTAINTY  = 10.0;  // cm (Y from timing along bar)

// TOF bar Y centres (special for bars 8 and 9)
auto tof_bar_y = [](int deviceID) -> double {
    int i = deviceID - 2053;
    if(i == 8) return -72.0;  // short lower bar
    if(i == 9) return +72.0;  // short upper bar
    return 0.0;               // all others centred at y=0
};
```

---

### Step 3: Module-Centre Reconstruction Function

**Location**: Add a `struct ModuleMeasurement` and a function `get_module_measurement()` before the event loop.

**What to do**: For any hit, return the module-centre position and the range of uncertainty in each transverse coordinate.

```cpp
struct ModuleMeasurement {
    double x_centre;     // best-estimate X (cm)
    double y_centre;     // best-estimate Y (cm)
    double z_pos;        // Z position of the detector plane (cm)
    double x_half_range; // ±half-width in X (cm); large = unconstrained
    double y_half_range; // ±half-width in Y (cm); large = unconstrained
    bool   valid;        // false if deviceID not recognised
};

// Z positions of detector planes (match the job directory T/F/E values)
// These should be read from the filename or set as constants per run.
// For a run with T1=240, T2=250, T3=570, T4=580, F1=260, F2=270, E1=600:
const double Z_TRACKER[4] = { T1_z, T2_z, T3_z, T4_z };  // fill from run params
const double Z_FRI[2]     = { F1_z, F2_z };
const double Z_TOF        = E1_z;

ModuleMeasurement get_module_measurement(int deviceID,
                                         const double* z_tracker,
                                         const double* z_fri,
                                         double z_tof)
{
    ModuleMeasurement m; m.valid = false;

    if(deviceID >= 1 && deviceID <= 1952)
    {
        // --- TRACKER STRAW ---
        int station   = (deviceID - 1) / 488;           // 0=T1/X, 1=T2/Y, 2=T3/U, 3=T4/V
        double lx     = tracker_local_x(deviceID);      // local X within station frame
        m.z_pos       = z_tracker[station];
        m.valid       = (m.z_pos > 0.);                 // disabled if z=0

        if(station == 0) {  // φ=0°, measures X
            m.x_centre     = lx;
            m.y_centre     = 0.;
            m.x_half_range = STRAW_HALF_WIDTH;
            m.y_half_range = STRAW_HALF_LENGTH;
        }
        else if(station == 1) {  // φ=−90°, measures Y
            m.x_centre     = 0.;
            m.y_centre     = lx;
            m.x_half_range = STRAW_HALF_LENGTH;
            m.y_half_range = STRAW_HALF_WIDTH;
        }
        else if(station == 2) {  // φ=+45°, stereo U
            // Constraint: x + y = sqrt(2)*lx (treat as diagonal strip)
            // For initial implementation, assign to both with large cross-term uncertainty
            m.x_centre     = lx / TMath::Sqrt2();
            m.y_centre     = lx / TMath::Sqrt2();
            m.x_half_range = STRAW_HALF_LENGTH / TMath::Sqrt2();
            m.y_half_range = STRAW_HALF_LENGTH / TMath::Sqrt2();
        }
        else {  // station==3, φ=−45°, stereo V
            m.x_centre     =  lx / TMath::Sqrt2();
            m.y_centre     = -lx / TMath::Sqrt2();
            m.x_half_range = STRAW_HALF_LENGTH / TMath::Sqrt2();
            m.y_half_range = STRAW_HALF_LENGTH / TMath::Sqrt2();
        }
    }
    else if(deviceID >= 2001 && deviceID <= 2052)
    {
        // --- FRI STRIP ---
        int fri_index = deviceID - 2001;
        int wall      = fri_index / 26;    // 0=Wall1 (measures X), 1=Wall2 (measures Y)
        int strip_i   = fri_index % 26;
        double lx     = FRI_STRIP_X[strip_i];
        double hw     = fri_half_width(strip_i);
        double hl     = FRI_HALF_LENGTH[strip_i];
        m.z_pos       = z_fri[wall];
        m.valid       = (m.z_pos > 0.);

        if(wall == 0) {  // strips run along Y, measures X
            m.x_centre     = lx;
            m.y_centre     = 0.;
            m.x_half_range = hw;
            m.y_half_range = hl;
        }
        else {           // rotated 90°, strips run along X, measures Y
            m.x_centre     = 0.;
            m.y_centre     = lx;
            m.x_half_range = hl;
            m.y_half_range = hw;
        }
    }
    else if(deviceID >= 2053 && deviceID <= 2070)
    {
        // --- TOF BAR ---
        m.z_pos       = z_tof;
        m.valid       = (m.z_pos > 0.);
        m.x_centre    = tof_bar_x(deviceID);
        m.y_centre    = tof_bar_y(deviceID);
        m.x_half_range = TOF_BAR_HALF_WIDTH;
        m.y_half_range = TOF_Y_UNCERTAINTY;
    }

    return m;
}
```

**Note on Z values**: Extract the T1/T2/T3/T4/F1/F2/E1 values from the input filename string (they appear in the JOB_LOG directory name and ROOT file name). Add a small helper to parse these, or pass them as macro arguments.

---

### Step 4: Hit Grouping, Intra-Component Intersection, and Coincidence Requirement

**Location**: Inside the event loop, after all tracker and FRI hits have been collected for a pion.

**New concept**: For each detector plane (tracker station pair, FRI wall pair), separately intersect all X-measuring fired elements on the X axis, and all Y-measuring fired elements on the Y axis. The intersection of multiple constraint regions on the same axis gives a **tighter position estimate** than any single element alone. If only one element fires that axis, its constraint is used directly.

The coincidence requirement is then: the event is reconstructable only if, for at least one detector plane, both the X-intersection AND the Y-intersection are valid (non-empty) — meaning the 2D overlap area exists.

**Helper functions to add (before the event loop):**

```cpp
// 1D interval defined by a centre and half-width.
struct Interval1D {
    double centre;     // midpoint of the interval
    double half_width; // half-length of the interval
    bool   valid;      // false if interval is empty
};

// Compute the intersection of a set of 1D constraint regions.
// Each region is supplied as {centre, half_width}.
// Returns the tightest interval consistent with ALL supplied regions.
// Returns valid=false if any pair is disjoint (inconsistent hits).
Interval1D intersect_constraints(
    const std::vector<std::pair<double,double>>& regions)
{
    if(regions.empty()) return {0., 9999., false};
    double lo = regions[0].first - regions[0].second;
    double hi = regions[0].first + regions[0].second;
    for(size_t i = 1; i < regions.size(); ++i) {
        double new_lo = regions[i].first - regions[i].second;
        double new_hi = regions[i].first + regions[i].second;
        lo = std::max(lo, new_lo);
        hi = std::min(hi, new_hi);
        if(lo >= hi) return {0., 0., false};  // disjoint — reject
    }
    return {0.5*(lo + hi), 0.5*(hi - lo), true};
}
```

**What to do in the event loop — for each pion, for each detector plane:**

```cpp
// Group all X-measuring hits at the same nominal Z (same station/wall) into one list,
// and all Y-measuring hits likewise.

// Example: for tracker module 1 (stations T1=X, T2=Y)
std::vector<std::pair<double,double>> t1_x_constraints;  // {centre, half_width}
std::vector<std::pair<double,double>> t2_y_constraints;

for(auto& hit : all_hits_pip) {
    ModuleMeasurement m = get_module_measurement(hit.deviceID, Z_TRACKER, Z_FRI, Z_TOF);
    if(!m.valid) continue;
    int station = (hit.deviceID >= 1 && hit.deviceID <= 1952)
                  ? (hit.deviceID - 1) / 488 : -1;
    if(station == 0)  // T1: X-measuring
        t1_x_constraints.push_back({m.x_centre, m.x_half_range});
    if(station == 1)  // T2: Y-measuring
        t2_y_constraints.push_back({m.y_centre, m.y_half_range});
    // Repeat for FRI Wall 1 (x), FRI Wall 2 (y), T3, T4 ...
}

// Compute intersections for this module
Interval1D x_from_T1 = intersect_constraints(t1_x_constraints);
Interval1D y_from_T2 = intersect_constraints(t2_y_constraints);

// Coincidence: both axes must be constrained at this detector plane
bool module1_valid = x_from_T1.valid && y_from_T2.valid;
```

Repeat this grouping and intersection for every detector plane (tracker module 2 at T3/T4, FRI at F1/F2).

**Coincidence requirement**: proceed with track fitting only if at least **two distinct detector modules** each contribute a valid 2D overlap area, giving enough Z separation to define a track direction.

---

### Step 5: 2D Overlap Area Construction from Intersection Results

**Location**: Replace the existing `hits_pip.push_back(smeared_pos)` logic and the subsequent `pip_dir` calculation.

**Concept**: Each detector module (a T1/T2 pair, a T3/T4 pair, or a F1/F2 pair) now contributes at most **one 2D overlap area** per pion track. This area is the Cartesian product of the intersected X constraint (from the module's X-measuring elements) and the intersected Y constraint (from the module's Y-measuring elements):

```
2D overlap area = [x_centre ± x_half_width] × [y_centre ± y_half_width]
```

where `x_centre`, `x_half_width` come from `intersect_constraints(all_X_constraints_in_module)` and similarly for Y. This is the region in which the particle is guaranteed to have passed, given all hits in that module.

The key difference from the previous approach (which used the centre of a single device):
- **Before**: `tp.x = single_straw_x_centre`, `tp.sx = straw_radius`
- **After**: `tp.x = intersection_of_all_fired_X_straws.centre`, `tp.sx = intersection_half_width` (smaller when multiple elements fire)

**Building `TrackPoint`s from the `Interval1D` intersection results (Step 4 output):**

```cpp
struct TrackPoint {
    double x, y, z;   // centre of the 2D overlap area
    double sx, sy;    // half-widths of the overlap area (uncertainties)
};

std::vector<TrackPoint> track_points_pip;

// For each detector module that produced a valid 2D overlap area:
struct DetectorModule {
    Interval1D x_constraint;  // from intersect_constraints on X-measuring hits
    Interval1D y_constraint;  // from intersect_constraints on Y-measuring hits
    double z_nominal;         // average Z of the X and Y planes in this module
};

// Example: a module built from T1 (X, z=T1_z) and T2 (Y, z=T2_z):
DetectorModule mod1;
mod1.x_constraint = intersect_constraints(t1_x_constraints);
mod1.y_constraint = intersect_constraints(t2_y_constraints);
mod1.z_nominal    = 0.5 * (Z_TRACKER[0] + Z_TRACKER[1]);

if(mod1.x_constraint.valid && mod1.y_constraint.valid) {
    TrackPoint tp;
    tp.x  = mod1.x_constraint.centre;
    tp.y  = mod1.y_constraint.centre;
    tp.z  = mod1.z_nominal;
    tp.sx = mod1.x_constraint.half_width;
    tp.sy = mod1.y_constraint.half_width;
    track_points_pip.push_back(tp);
}

// Repeat for the FRI module (F1=X at z_fri[0], F2=Y at z_fri[1]):
DetectorModule fri_mod;
fri_mod.x_constraint = intersect_constraints(f1_x_constraints);
fri_mod.y_constraint = intersect_constraints(f2_y_constraints);
fri_mod.z_nominal    = 0.5 * (Z_FRI[0] + Z_FRI[1]);

if(fri_mod.x_constraint.valid && fri_mod.y_constraint.valid) {
    TrackPoint tp;
    tp.x  = fri_mod.x_constraint.centre;
    tp.y  = fri_mod.y_constraint.centre;
    tp.z  = fri_mod.z_nominal;
    tp.sx = fri_mod.x_constraint.half_width;
    tp.sy = fri_mod.y_constraint.half_width;
    track_points_pip.push_back(tp);
}
```

**Why this is strictly better than the previous single-device centre approach:**

| Scenario | Old `tp.sx` | New `tp.sx` |
|----------|------------|------------|
| One FRI outer strip fires | 3.0 cm | 3.0 cm (unchanged) |
| Two adjacent outer FRI strips fire | 3.0 cm (only first used) | ≤ 3.0 cm (intersection, narrower at boundaries) |
| Two non-adjacent FRI strips fire | 3.0 cm (first used) | intersection likely empty → module excluded |
| Tracker: one straw fires from sub-layer 1 | 0.4 cm | 0.4 cm (unchanged) |
| Tracker: two straws from offset sub-layers | 0.4 cm (only first used) | < 0.4 cm (intersection of overlapping regions) |

For the TOF wall there is no intersection to compute (single layer, one bar per pion per event); the bar geometry lookup from Step 3 is used directly as the X and Y uncertainty.

> **Note on degenerate cases**: If the intersection of all X constraints returns `valid=false` (inconsistent fired elements, e.g. two widely separated straws from separate particles), that module cannot contribute a TrackPoint. The event may still be reconstructable from the remaining modules.

---

### Step 6: Line-of-Best-Fit Track Reconstruction

**Location**: Replace the current `pip_dir = (hits_pip.back() - hits_pip.front()).Unit()` calculation.

**Concept**: With 2 or more `TrackPoint`s along Z, fit a straight line in (x vs z) and (y vs z) independently using weighted linear regression. The weights are `1/σ²` for each point.

```cpp
TVector3 fit_track_direction(const std::vector<TrackPoint>& pts)
{
    // Need at least 2 points to define a line
    if(pts.size() < 2) return TVector3(0,0,1);

    // Weighted linear regression: fit x = ax*z + bx, y = ay*z + by
    double sw=0, swz=0, swz2=0, swx=0, swzx=0, swy=0, swzy=0;

    for(const auto& p : pts)
    {
        double wx = 1.0 / (p.sx * p.sx + 1e-6);
        double wy = 1.0 / (p.sy * p.sy + 1e-6);
        double w  = 0.5 * (wx + wy);   // combined weight

        sw   += w;
        swz  += w * p.z;
        swz2 += w * p.z * p.z;
        swx  += wx * p.x;
        swzx += wx * p.z * p.x;
        swy  += wy * p.y;
        swzy += wy * p.z * p.y;
    }

    double denom = sw * swz2 - swz * swz;
    if(std::abs(denom) < 1e-12) return TVector3(0,0,1);

    double ax = (sw * swzx - swz * swx) / denom;   // dx/dz slope
    double ay = (sw * swzy - swz * swy) / denom;   // dy/dz slope

    return TVector3(ax, ay, 1.0).Unit();
}

// Track origin: evaluate at the mean Z of all points
TVector3 fit_track_origin(const std::vector<TrackPoint>& pts, TVector3 dir)
{
    // Use centroid of all track points as a point on the line
    double mx=0, my=0, mz=0;
    for(const auto& p : pts){ mx+=p.x; my+=p.y; mz+=p.z; }
    int n = pts.size();
    return TVector3(mx/n, my/n, mz/n);
}
```

**Replacing the existing `pip_dir` and `pip_start` assignments:**

```cpp
// OLD:
// TVector3 pip_dir  = (hits_pip.back() - hits_pip.front()).Unit();
// TVector3 pip_start = hits_pip.front();

// NEW:
TVector3 pip_dir   = fit_track_direction(track_points_pip);
TVector3 pip_start = fit_track_origin(track_points_pip, pip_dir);
```

The PoCA call (`closest_point_between_lines`) then proceeds unchanged.

---

### Step 7: TOF Wall Treatment

**Location**: The section that reads TOF hits and computes `hit_tof_pip_t` (timing) and `reco_x`, `reco_y` for the TOF hit used in kaon momentum calculation.

**What to do**: For TOF bars, the X position is constrained to the bar centre ±6 cm, and the Y position is the bar Y centre ±10 cm. The timing `hitT` is still used as before for velocity calculation. Only the positional part changes:

```cpp
// When processing a TOF hit (is_tof(deviceID)):
ModuleMeasurement tof_m = get_module_measurement(deviceID, Z_TRACKER, Z_FRI, Z_TOF);

double tof_x = tof_m.x_centre;   // bar centre X, no smear needed
double tof_y = tof_m.y_centre    // bar Y centre
               + randGen.Uniform(-TOF_Y_UNCERTAINTY, +TOF_Y_UNCERTAINTY);
               // OR: use Gaussian smear with σ=TOF_Y_UNCERTAINTY
double tof_z = Z_TOF;
// Use tof_x, tof_y, tof_z for the TOF hit position
// Use hitT as before for timing
```

> The ±10 cm Y uncertainty can be modelled as a Gaussian `Gaus(0, 10)` or as a uniform `Uniform(-10, 10)`. Gaussian is more physically motivated (timing resolution along bar → central hits more precise).

---

### Step 8: Mirror Changes to Acceptance Macro

After validating all steps above in `KLong_save_vectors.C`, apply the **identical set of changes** to `KLong_save_momentum_acceptance.C`. The reconstruction logic in both files is parallel; all the same lookup tables, coincidence requirement, overlap construction, and line-of-best-fit changes apply.

The only difference is the outer loop structure (`reco_flag` assignment vs event writing), which is unchanged.

---

## 5. Pipeline Integration Notes

### Extracting Z Positions from the Filename

Each ROOT file is named after the job configuration, e.g.:  
`T1-240_T2-250_T3-0_T4-0_P1-215_P2-230_F1-260_F2-270_E1-700`

The macro already receives the filename as input. Parse it like this:

```cpp
// Example: extract T1 value from filename string
auto extract_param = [](const TString& fname, const TString& key) -> double {
    TRegexp re(key + "-[0-9]+");
    TString token = fname(re);
    if(token.Length() == 0) return 0.;
    token.Remove(0, key.Length() + 1);
    return token.Atof();
};

double T1 = extract_param(fname, "T1");
double T2 = extract_param(fname, "T2");
double T3 = extract_param(fname, "T3");
double T4 = extract_param(fname, "T4");
double F1 = extract_param(fname, "F1");
double F2 = extract_param(fname, "F2");
double E1 = extract_param(fname, "E1");

const double Z_TRACKER[4] = { T1, T2, T3, T4 };
const double Z_FRI[2]     = { F1, F2 };
const double Z_TOF        = E1;
```

This ensures the reconstruction automatically uses the correct Z positions for every run configuration without hardcoding.

### `detector_simulation_master_FIXED.sh`

No changes required to the shell script if you are only modifying `KLong_save_vectors.C` and `KLong_save_momentum_acceptance.C`. The script already copies these files to each job's parsing directory.

If you rename the files or create new versions (e.g. `KLong_save_vectors_v2.C`), update the `cp` commands in the script that copy the parsing macros.

---

## 6. Validation Checklist

Work through these checks in order, ideally on a single small ROOT file before running the full pipeline.

### Check 1: Geometry lookup sanity
- Print `get_module_measurement(1, ...)` → should give `x_centre = -48.4`, `y_half_range = 45.0`
- Print `get_module_measurement(489, ...)` → should give `y_centre = -48.4`, `x_half_range = 45.0` (first straw of T2/Y station)
- Print `get_module_measurement(2001, ...)` → should give `x_centre = -66.0` (first FRI Wall 1 strip)
- Print `get_module_measurement(2027, ...)` → should give `y_centre = -66.0` (first FRI Wall 2 strip)
- Print `get_module_measurement(2061, ...)` → should give `x_centre = 0.0`, `y_centre = -72.0` (TOF bar i=8)

### Check 2: Coincidence requirement and intersection logic
- Run on a sample file, count events passing the coincidence criterion (at least two detector modules each contributing a valid 2D overlap area).
- Compare to rates in old macro. A modest reduction is expected (stricter requirement), but the rate should not drop to near zero.
- **Intersection sanity check**: add a debug printout for events where two FRI Wall-1 strips fire simultaneously. The reported `x_centre` and `x_half_width` should be narrower than either individual strip's half-width. Example: outer strip at x=−12 cm (hw=3 cm) and adjacent outer strip at x=−18 cm (hw=3 cm) give intersection [−15, −15] → `centre=−15, half_width=0` (boundary hit). A strip at −12 (hw=3) and −6 (hw=3) give [−9,−9] similarly. Any two non-adjacent strips give an empty intersection → that module is excluded.
- **Tracker intersection sanity check**: once sub-layer offsets are confirmed from `JLabKDetectorConstruction.cc`, verify that two straws separated by ~0.4 cm (one pitch offset) produce a non-empty intersection; two straws separated by ≥0.8 cm (no sub-layer overlap) produce an empty intersection.

### Check 3: Track direction comparison
- For events reconstructed by both old and new methods, compare the angle between `pip_dir_old` and `pip_dir_new`.
- Without detector misalignment, the directions should cluster near 0° difference with RMS matching the smearing scale.

### Check 4: PoCA vertex distribution
- Compare the `reco_vtx` distribution (from PoCA) between old and new methods.
- The distribution should be centred at the same place but with different widths (new method may be narrower if geometry constraints are better than random Gaussian smearing, or wider if the module long-axis uncertainty dominates).

### Check 5: Kaon momentum distribution
- Compare `kaon_p` distributions. Shape should be similar; any large shift indicates a systematic error in the new geometry lookup.

### Check 6: Multi-configuration consistency
- Run on two configurations with different detector Z positions (e.g. `T1-240...` and `T1-280...`) and verify that the reconstructed vertex distributions shift appropriately with Z position.

### Check 7: Mirror to acceptance macro
- After all above pass, run `KLong_save_momentum_acceptance.C` on the same sample file.
- Verify that the `reco_flag` efficiency matches expectations and that `true_mom_vec` is unchanged (it should be — it's from simulation, not reconstruction).

---

## 7. DSSSD Simplification and Arrival-Angle Correction (March 2026)

**Status**: Implemented in `KLong_save_vectors.C` and `KLong_save_momentum_acceptance.C`  
**Motivation**: Supervisor guidance (email, March 2026) — detailed below.

### 7.1 Background and Rationale

The original module-based reconstruction treated each tracker station independently:
- T1 (X-measuring, z=240) → one TrackPoint with tight X, unconstrained Y (±45 cm)
- T2 (Y-measuring, z=250) → one TrackPoint with tight Y, unconstrained X (±45 cm)

This gives two weakly-constrained points rather than one tight 2D pixel.  The
practical concerns raised by the JLab group are:

1. **Electronics**: No one at JLab expects all 488 straws per station to have
   individual readout. Realistically at most one layer of X-tubes and one layer of
   Y-tubes would be read out at any given station.
2. **Complexity**: The 4-sub-layer intersection logic is correct in principle but adds
   unnecessary code complexity for the current design study.
3. **FRI/ToF priority**: The FRI walls and TOF wall are far more likely to be
   implemented in the long-term experiment than the tracker itself; these should be
   the primary focus for reconstruction quality.

### 7.2 What Changed: Sub-Layer 0 Only

Tracker hits are now filtered to **sub-layer 0 only** before reconstruction.

Sub-layer 0 straws (straw_cn 1–122 within each 488-ID station block):

| Station | Sub-layer-0 ID range | Measures |
|---------|---------------------|---------|
| T1      | 1–122               | X       |
| T2      | 489–610             | Y       |
| T3      | 977–1098            | U (+45°)|
| T4      | 1465–1586           | V (−45°)|

Hits with straw_cn > 122 within a station block (sub-layers 1–3) are discarded
before the intersection logic runs.

### 7.3 What Changed: DSSSD Pixel Formation

Instead of producing one single-axis TrackPoint per station, stations are now
**paired as DSSSD units** — one tight 2D pixel per detector pair:

| DSSSD pair | X source | Y source | Mid-Z |
|------------|----------|----------|-------|
| Tracker    | T1 (key=0, X)  | T2 (key=1, Y)  | (T1+T2)/2 |
| FRI        | F1 (key=10, X) | F2 (key=11, Y) | (F1+F2)/2 |
| Stereo     | T3 (u, via unfolding) | T4 (v, via unfolding) | (T3+T4)/2 |

For the **stereo T3+T4 pair**, proper u/v unfolding is applied:

```
u = (x+y)/√2  [measured by T3]   v = (x-y)/√2  [measured by T4]
→  x = (u+v)/√2                   y = (u-v)/√2
   σ_x = σ_y = hypot(σ_u, σ_v) / √2
```

This replaces the previous approximation where T3 and T4 were treated with
equal x/y half-ranges and never formed a proper pixel.

When only one axis of a pair fires, a single-axis point is still recorded
(tight axis, free axis = STRAW\_HALF\_LENGTH or 70.25 cm) — this matches the
old behaviour for unpaired hits.

**Minimum requirement**: >= 2 TrackPoints are still required, but now each
point is typically a tight 2D pixel.  A typical reconstructable event produces:
- 1 tight pixel from the tracker DSSSD (T1+T2), OR only 1 axis-constrained
  point if T2 or T1 alone fires
- 1 tight pixel from the FRI DSSSD (F1+F2)
- Optionally 1 tight stereo pixel from T3+T4

### 7.4 What Changed: Arrival-Angle Correction

When T1 and T2 are 10 cm apart in z, a particle with transverse momentum
will have its X position at z=240 and Y position at z=250 at slightly
different transverse coordinates from the common pixel.  Correcting both to
the mid-z (z=245) requires knowing the slope dx/dz and dy/dz.

**Implementation**: The MC-truth kaon momentum direction is used:

```
truth_dxdz = true_px / true_pz
truth_dydz = true_py / true_pz

x_pixel = x_T1 + truth_dxdz × (z_mid − z_T1)
y_pixel = y_T2 + truth_dydz × (z_mid − z_T2)
```

This is an approved approximation for the simplified tracker model:
> *"In terms of adding an arrival angle, how do you currently implement the
>  arrival of a particle? Do you have direction information in the truth
>  4-vector? If so, can you just use that?"*
> — Supervisor email, March 2026

The same correction applies to the FRI DSSSD pair (F1 at z=260, F2 at z=270).
For the T3+T4 stereo pair the correction is omitted because the unfolding
geometry makes the angle terms couple across both axes non-trivially, and the
effect is small.

### 7.5 Low Acceptance: Interpretation from Parse Logs

Examination of the debug parse logs (20 debug events per job file) reveals that
**most events with 0 detector hits are caused by genuine geometric acceptance**,
not by reconstruction bugs:

- A K_L decaying at z >> 0 (before the first tracker at z=240) produces pions
  that must travel forward to reach the detectors. K_L with modest transverse
  momentum can produce pions that miss the ±48.4 cm tracker aperture entirely.
- K_L can also decay far downstream (past the FRI wall), producing pions that
  never reach any tracker plane.
- Events with nonzero tracker hits BUT PIZZA=NO typically have the pion
  scattered or absorbed before reaching the PIZZA detector.

The DSSSD simplification does not change geometric acceptance; it improves
reconstruction quality for events that *do* pass through the detectors.

---

## 8. TGraph-Based Track Fitting — Raw-Position Approach (March 2026)

This section documents the second major overhaul of the tracking code.  In
response to supervisor feedback the geometry-aware DSSSD pixelation approach
(Section 7) was replaced as the *default* by a simpler, more robust method
that fits straight lines directly to the raw Geant4 3D hit positions using
ROOT's `TGraphErrors`.  The DSSSD pixelation code is preserved in the
`_pixelation.C` companion files.

---

### 8.1 Motivation and Differences from the Pixelation Approach

The pixelation approach required:
- A large module-centre geometry lookup table (straw positions, FRI strip
  centres, TOF bar centres).
- An intersection/pairing algorithm to produce "pixels" from pairs of
  perpendicularly-oriented DSSSD strips (T1+T2, F1+F2, T3+T4 stereo).
- An arrival-angle correction to the projected strip measurement.
- Restriction to sub-layer 0 only (to avoid producing false cross-pixels
  from staggered sub-layers).
- Special T3+T4 stereo unfolding via rotation matrix.

**Why switch to TGraphErrors:**

| Concern | Pixelation | TGraphErrors |
|---------|-----------|-------------|
| Geometry tables needed? | Yes — full straw/strip/bar look-up | Only TOF bar centre |
| Sub-layer restriction? | Yes — sub-layer 0 only | No — all sub-layers used |
| T3/T4 stereo handling | u/v rotation, paired pixel | Raw (x,y) + inflated σ |
| Arrival-angle correction | Required | Not needed (raw Geant4 pos.) |
| Code complexity | ~400 lines of helper functions | ~80 lines for core function |
| Fit quality | Line-of-best-fit via least-squares over pixels | Weighted χ² fit via `pol1` |

The key insight is that the raw Geant4 positions already contain the correct
spatial information; the geometry lookup was only needed to *recover* the
measurement from an ID-based encoding.  Fitting a weighted line through all
hits simultaneously is statistically superior to using only a first/last pair.

---

### 8.2 The `TrackFit` Struct

Each reconstructed track is described by:

```cpp
struct TrackFit {
    TVector3 dir;    // unit direction vector  (dx/dz, dy/dz, 1).Unit()
    TVector3 orig;   // reference point at z = 0:  orig = (bx, by, 0)
    bool     valid;  // false when fit failed or z-span was insufficient
};
```

The parameterisation is:

$$\mathbf{r}(z) = \mathbf{orig} + z\,\frac{\mathbf{dir}}{dir_z}
               = (b_x + a_x z,\; b_y + a_y z,\; z)$$

where $a_x = dx/dz$, $a_y = dy/dz$ are the fitted slopes and $(b_x, b_y)$
are the fitted intercepts at $z=0$.  The unit direction is then
$(a_x, a_y, 1)\,/\,\|(a_x, a_y, 1)\|$.

The downstream PoCA (`closest_point_between_lines`) is unchanged and consumes
`TrackFit::dir` and `TrackFit::orig` identically to how it consumed the old
pixelation-derived direction/origin.

---

### 8.3 `fit_3dgraph_track()`: Implementation Details

**Signature:**
```cpp
TrackFit fit_3dgraph_track(const std::vector<HitInfo>& hits,
                           int    tof_devID,
                           double z_tof);
```

**Algorithm (step by step):**

1. **Accumulate data vectors.**  Loop over every hit in `hits`.  For each hit,
   extract the raw Geant4 position $(h.x, h.y, h.z)$ and decide which
   projection(s) it constrains based on its `deviceID` range (see §8.4).
   Push `(z, value, sigma)` into one or both of:
   - `(xz, xv, xe)` — the x(z) data set
   - `(yz, yv, ye)` — the y(z) data set

2. **Add the TOF bar anchor.**  If `tof_devID` is valid (2053–2070) and
   `z_tof > 0`, append a point to *both* data sets using the geometry-derived
   bar centre `(tof_bar_x_centre, tof_bar_y_centre)` with uncertainties
   `TOF_BAR_HALF_WIDTH = 6 cm` (X) and `TOF_Y_UNCERTAINTY = 10 cm` (Y).
   This anchors the far end of the fit to the known bar position.

3. **Validity pre-check.**  Return `FAIL` immediately if either data set has
   fewer than 2 entries.

4. **Z-span check.**  Return `FAIL` if all x-projection points lie at the same
   z-plane (within 1 cm), or all y-projection points likewise.  This guards
   against a degenerate fit with undefined slope.

5. **Construct `TGraphErrors` and fit.**
   ```cpp
   TGraphErrors gx(N, xz.data(), xv.data(), nullptr, xe.data());
   gx.Fit("pol1", "Q");        // quiet mode
   TF1* fx = gx.GetFunction("pol1");
   double bx = fx->GetParameter(0);   // intercept
   double ax = fx->GetParameter(1);   // slope dx/dz
   ```
   The same is done for `gy`/`gy` to obtain $(b_y, a_y)$.

6. **Sanity check.**  If any of the four parameters is non-finite (fit
   diverged), return `FAIL`.

7. **Assemble and return `TrackFit`.**
   ```cpp
   TVector3 dir = TVector3(ax, ay, 1.0).Unit();
   TVector3 orig(bx, by, 0.0);
   return {dir, orig, true};
   ```

---

### 8.4 Per-Detector Uncertainty Table

The error assigned to each hit point determines the weight it receives in the
χ² minimisation.  The values are chosen to reflect the physical measurement
resolution of each sub-detector.

| Detector | ID Range | Projection(s) | σ (cm) | Rationale |
|----------|----------|---------------|--------|-----------|
| T1 (station 0) | 1–488     | x(z) only | 0.40 | Straw half-width; straw runs ‖ Y so only X is measured |
| T2 (station 1) | 489–976   | y(z) only | 0.40 | Straw runs ‖ X so only Y is measured |
| T3 (station 2) | 977–1464  | x(z) **and** y(z) | 0.40 × √2 ≈ 0.566 | Stereo +45°; see §8.5 |
| T4 (station 3) | 1465–1952 | x(z) **and** y(z) | 0.40 × √2 ≈ 0.566 | Stereo −45°; see §8.5 |
| FRI Wall 1     | 2001–2026 | x(z) only | 1.5 or 3.0 | Strip half-width; inner 6 strips = 3 cm, middle 6 = 1.5 cm |
| FRI Wall 2     | 2027–2052 | y(z) only | 1.5 or 3.0 | Same strip scheme, perpendicular orientation |
| TOF bar        | 2053–2070 | x(z) **and** y(z) | 6.0 (X), 10.0 (Y) | Bar half-width / timing resolution |

The constant `STRAW_HALF_WIDTH = 0.40 cm` is defined at file scope.

For FRI strips, `fri_half_width(strip_i)` returns 1.5 cm for strips 10–15
(the central, narrower strips) and 3.0 cm for strips 0–9 and 16–25.

---

### 8.5 T3/T4 Stereo Uncertainty Propagation

T3 straws run at +45° and T4 straws run at −45° relative to the Z-axis in
the XY-plane.  Each straw measures the perpendicular displacement in its
rotated coordinate frame with precision $\sigma_0 = 0.40$ cm.

Calling the rotated measurement axis $\hat{u}$ (where $|\hat{u}|=1$ and
$\hat{u}$ is perpendicular to the straw), the lab-frame components are:

$$\sigma_x = \sigma_0|\cos\phi|, \quad \sigma_y = \sigma_0|\sin\phi|$$

For $\phi = \pm 45°$: $|\cos\phi| = |\sin\phi| = 1/\sqrt{2}$, so projecting
the constraint onto lab X forces $\sigma_x = \sigma_0/\sqrt{2}$.  But because
the TGraph approach uses the raw Geant4 hit position (not a projected
measurement), the hit is entered into **both** the x(z) and y(z) fits.  The
error to assign is then the uncertainty on the raw position arising from the
±rotated-measurement precision, which — propagated onto either lab axis by
the rotation — gives:

$$\sigma_{\rm lab} = \sigma_0 \times \sqrt{2}$$

This is the value stored in `xe` and `ye` for T3/T4 hits, implemented as
`STRAW_HALF_WIDTH * SQRT2` where `SQRT2 = std::sqrt(2.0)`.

**Contrast with the pixelation approach**, where T3+T4 was treated as a
stereo pair: the two straw measurements were geometrically inverted
(via a 2×2 rotation matrix) to obtain a pixel centre in lab (x, y), and only
that pixel — not both raw hits — was passed to the fitter.  The TGraph
approach foregoes this inversion; instead both raw hits contribute, each with
their inflated uncertainty, letting the χ² minimiser absorb the redundancy.

---

### 8.6 TOF Bar: Geometry Lookup vs Raw Position

Unlike all other detector hits (which use the raw Geant4 `h.x, h.y`
positions), the TOF bar uses the **geometry-derived bar centre**:

```cpp
tof_bar_x_centre(tof_devID)  // = -96 + (id - 2053) * 12  cm
tof_bar_y_centre(tof_devID)  // = 0 cm for most bars; ±72 cm for bars 8 and 9
```

This is used for two reasons:
1. The TOF scintillator bar is large ($\approx 12$ cm wide, $\approx 100$ cm
   tall) and the Geant4 hit position within it has little track-fitting value.
2. The bar position narrows down which of 18 bars was hit; the bar *centre* is
   the best available position estimate given the bar geometry.

The assigned uncertainties ($\sigma_x = 6$ cm half-bar-width,
$\sigma_y = 10$ cm) ensure the TOF point contributes a soft anchor at the far
end of the fit rather than dominating it.

---

### 8.7 Fit Validity Conditions

A `TrackFit` is marked `valid = false` (returned as `FAIL`) under any of the
following conditions:

| Condition | Rationale |
|-----------|-----------|
| Fewer than 2 x-projection points | Cannot fit a line with a single point |
| Fewer than 2 y-projection points | Same |
| All x-projection points at the same z (within 1 cm) | Slope is undefined (division by zero in χ² sense) |
| All y-projection points at the same z (within 1 cm) | Same |
| ROOT `pol1` fit returns null function pointer | Fit did not converge |
| Any of $a_x, a_y, b_x, b_y$ is non-finite | Fit diverged numerically |

Events whose pion produces a `valid = false` track for either the positive or
negative pion are dropped from the output tree (same behaviour as the
pixelation approach).

---

### 8.8 All Sub-Layers Enabled

The pixelation approach required filtering to **sub-layer 0 only** (122 straws
per station) to avoid generating spurious cross-pixels from the staggered
sub-layers.  With the TGraph approach this restriction is removed.

All four sub-layers of each tracker station (488 straws in T1, 488 in T2,
488 in T3, 488 in T4) contribute hit points.  More points improve the
weighted least-squares fit; the staggered sub-layers simply add more
constraints at slightly different z-positions (the sub-layer z-offsets are
small compared to the station spacing but are handled naturally by the fit).

In code, the ID ranges used cover the full straw complement:
- T1: IDs 1–488 (all 4 sub-layers)
- T2: IDs 489–976
- T3: IDs 977–1464
- T4: IDs 1465–1952

No masking or modular-arithmetic sub-layer extraction is performed.

---

### 8.9 File Structure After Overhaul

After all three overhauls (original → DSSSD simplification → TGraph), the
repository contains the following data-parsing macros:

| File | Reconstruction method | Status |
|------|-----------------------|--------|
| `KLong_save_vectors.C` | TGraph (`fit_3dgraph_track`), all sub-layers | **Current default** |
| `KLong_save_momentum_acceptance.C` | TGraph (mirrors vectors), all sub-layers | **Current default** |
| `KLong_save_vectors_pixelation.C` | DSSSD pixelation, sub-layer 0, angle correction | Preserved backup |
| `KLong_save_momentum_acceptance_pixelation.C` | DSSSD pixelation (mirrors vectors) | Preserved backup |

Both defaults compile cleanly under ROOT 6.32.20.  The output branch
structure is identical across both approaches, so all downstream histogram
and combination macros work without modification.

---

## Appendix: Files to Modify

| File | Changes | Status |
|------|---------|--------|
| `KLong_save_vectors.C` | **Current default** — TGraph approach (§8): `fit_3dgraph_track`, all 4 sub-layers, `TrackFit` struct. Also contains all Steps 1–7 and Section 7 DSSSD changes underneath. | ✅ Implemented |
| `KLong_save_momentum_acceptance.C` | **Current default** — mirrors `KLong_save_vectors.C` with TGraph reconstruction; records all pi+pi-pi0 events with `reco_flag`. | ✅ Implemented |
| `KLong_save_vectors_pixelation.C` | **Backup** — DSSSD pixelation approach: module-centre geometry lookup, DSSSD pairing, sub-layer 0 filter, arrival-angle correction, stereo T3+T4 unfolding. Preserved from Section 7 overhaul. | ✅ Preserved |
| `KLong_save_momentum_acceptance_pixelation.C` | **Backup** — DSSSD acceptance version of the pixelation approach. | ✅ Preserved |
| `JLabKSteppingAction.cc` | **Updated** (March 2026) — encodes `(GetCopyNumber(1)−620)×488 + GetCopyNumber(0)` for straws; shifts PIZZA/FRI/TOF by +1464 | ✅ Done |
| `JLabKDetectorConstruction.cc` | **No change required** | — |
| `detector_simulation_master_FIXED.sh` | No change unless files are renamed | — |

## Appendix: Device ID Quick Reference

| ID Range   | Detector           | Measures          | TGraph fit | Notes |
|------------|--------------------|-------------------|------------|-------|
| 1–488      | Tracker station 0 (T1) | X (φ=0°)     | x(z) only | All 4 sub-layers; σ = 0.4 cm. Pixelation backup uses sub-layer 0 only (IDs 1–122) |
| 489–976    | Tracker station 1 (T2) | Y (φ=−90°)   | y(z) only | All 4 sub-layers; σ = 0.4 cm. Pixelation backup uses IDs 489–610 |
| 977–1464   | Tracker station 2 (T3) | U (φ=+45°)   | x(z) **and** y(z) | All 4 sub-layers; σ = 0.4√2 cm. Pixelation pairs T3+T4 for stereo pixel |
| 1465–1952  | Tracker station 3 (T4) | V (φ=−45°)   | x(z) **and** y(z) | All 4 sub-layers; σ = 0.4√2 cm. Pixelation pairs with T3 |
| 1953–1976  | PIZZA flat         | angular           | — | Not used in track fit |
| 1977–2000  | PIZZA conical      | angular           | — | Not used in track fit |
| 2001–2026  | FRI Wall 1         | X (φ=0°)          | x(z) only | σ = 1.5 cm (strips 10–15) or 3.0 cm (all others) |
| 2027–2052  | FRI Wall 2         | Y (φ=+90°)        | y(z) only | Same strip σ scheme |
| 2053–2070  | TOF Wall           | X (bar), Y (±10 cm) | x(z) **and** y(z) | Bar centre from geometry; σ_x = 6 cm, σ_y = 10 cm |
