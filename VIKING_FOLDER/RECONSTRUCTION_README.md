# K_L Decay Reconstruction — Annotated Technical Reference

> **Purpose:** This document walks through every step of the event reconstruction
> in `KLong_save_vectors.C`, tracing each numerical quantity back to its exact
> source and spelling out every equation with full annotation.
> The focus is the **decay vertex determination**, which is the critical
> intermediate that feeds directly into the kaon momentum calculation.
> A systematic overestimation of the vertex z coordinate has been identified
> as the source of an underestimated kaon momentum.

---

## Table of Contents

1. [Detector Layout and Coordinate System](#1-detector-layout-and-coordinate-system)
2. [Data Sources — Where Every Number Comes From](#2-data-sources--where-every-number-comes-from)
3. [Stage 1 — Event Selection from Ntuple2](#3-stage-1--event-selection-from-ntuple2)
4. [Stage 2 — Track Fitting (TGraphErrors pol1)](#4-stage-2--track-fitting-tgrapherrors-pol1)
   - 4.1 [What goes into the fit and what doesn't](#41-what-goes-into-the-fit-and-what-doesnt)
   - 4.2 [x(z) and y(z) fits](#42-xz-and-yz-fits)
   - 4.3 [What the fit output means](#43-what-the-fit-output-means)
   - 4.4 [The TrackFit struct](#44-the-trackfit-struct)
5. [Stage 3 — Decay Vertex via Point of Closest Approach (PoCA)](#5-stage-3--decay-vertex-via-point-of-closest-approach-poca)
   - 5.1 [Two-line PoCA derivation](#51-two-line-poca-derivation)
   - 5.2 [What is actually computed](#52-what-is-actually-computed)
   - 5.3 [The midpoint convention and its implications](#53-the-midpoint-convention-and-its-implications)
   - 5.4 [Why small slope errors produce large z errors](#54-why-small-slope-errors-produce-large-z-errors)
6. [Stage 4 — PIZZA Hit Position](#6-stage-4--pizza-hit-position)
7. [Stage 5 — TOF Hit Position](#7-stage-5--tof-hit-position)
8. [Stage 6 — Pion Velocity](#8-stage-6--pion-velocity)
9. [Stage 7 — Kaon Decay Time (Back-Propagation)](#9-stage-7--kaon-decay-time-back-propagation)
10. [Stage 8 — Kaon Momentum](#10-stage-8--kaon-momentum)
11. [Error Propagation Summary and Bias Pathways](#11-error-propagation-summary-and-bias-pathways)

---

## 1. Detector Layout and Coordinate System

The beam travels in the **+z direction**. All positions are in the lab frame in
**centimetres**, all times in **nanoseconds** unless otherwise stated.

The kaon is produced at the **origin** `(0, 0, 0)` (hard-coded assumption;
see [Stage 8](#10-stage-8--kaon-momentum)).

```
z = 0           kaon production point (FIXED: hard-coded as origin)
z ~ 50–300 cm   expected range of K_L decay vertices
z ~ 240 cm      T1 (tracker station 0, measures X)
z ~ 250 cm      T2 (tracker station 1, measures Y)
z ~ 260 cm      F1 (FRI Wall 1, measures X) — depends on config tag
z ~ 270 cm      F2 (FRI Wall 2, measures Y) — depends on config tag
z ~ 210–230 cm  PIZZA detectors (TOF start)  — depends on config tag (P1, P2)
z ~ 680 cm      T3 (tracker station 2, stereo +45°)
z ~ 690 cm      T4 (tracker station 3, stereo -45°)
z ~ 700 cm      TOF wall (TOF stop)           — from E1 tag in filename
```

**Device ID ranges** (from `JLabKDetectorConstruction.cc`):

| Range        | Detector          | Measures |
|--------------|-------------------|----------|
| 1 – 488      | T1                | X        |
| 489 – 976    | T2                | Y        |
| 977 – 1464   | T3 (stereo +45°)  | X and Y  |
| 1465 – 1952  | T4 (stereo -45°)  | X and Y  |
| 1953 – 1976  | PIZZA flat        | time     |
| 1977 – 2000  | PIZZA conical     | time     |
| 2001 – 2026  | FRI Wall 1        | X        |
| 2027 – 2052  | FRI Wall 2        | Y        |
| 2053 – 2070  | TOF bars          | time, X  |

---

## 2. Data Sources — Where Every Number Comes From

The ROOT file contains three TTrees.

### Ntuple1 — MC Truth

One row per event. Provides the **ground truth** for validation only; none of
these quantities are fed into the reconstruction.

| Branch             | Type   | What it is                                      |
|--------------------|--------|-------------------------------------------------|
| `evtNb`            | double | Event number (integer stored as double)         |
| `kaonDK_momX/Y/Z`  | double | True kaon momentum at decay point (GeV/c)       |
| `kaonDK_posX/Y/Z`  | double | True K_L decay vertex position (cm)             |

Stored as `TruthInfo { px, py, pz, vx, vy, vz }` in `truth_map[evtNb]`.

### Ntuple2 — Decay Products

One row per secondary particle per event. Used only for **event selection**.

| Branch                    | Type   | What it is                      |
|---------------------------|--------|---------------------------------|
| `evtNb`                   | double | Event number                    |
| `DKparticle_PDGEncoding`  | double | PDG code of each decay product  |

The macro selects only events where the set of decay products contains
**PDG 211 (π⁺), −211 (π⁻), and 111 (π⁰)** — i.e. `K_L → π⁺ π⁻ π⁰`.

### Ntuple3 — Raw Detector Hits

One row per detector hit. This is the **sole source of all reconstructed
quantities**. Every hit carries:

| Branch        | Type   | What it is                                            |
|---------------|--------|-------------------------------------------------------|
| `evtNb`       | double | Event number                                          |
| `hitX`        | double | Geant4 x position of the hit in lab frame (cm)        |
| `hitY`        | double | Geant4 y position of the hit in lab frame (cm)        |
| `hitZ`        | double | Geant4 z position of the hit in lab frame (cm)        |
| `hitT`        | double | Geant4 hit time in lab frame (ns)                     |
| `PDGEncoding` | double | PDG code of the particle that made the hit            |
| `deviceID`    | double | Integer device ID identifying which detector fired    |

> **Critical note:** `hitX`, `hitY`, `hitZ` are the **raw Geant4 truth positions**
> of the particle when it passed through the detector. There is **no digitisation
> model** — no smearing of x/y/z for tracker or FRI hits. The only smearing
> applied is to **timing**: `t_smeared = Gaus(hitT, 1.5 ps)` for PIZZA and TOF.
> TOF bar x/y are derived from the bar geometry, not from hitX/hitY directly.

---

## 3. Stage 1 — Event Selection from Ntuple2

```cpp
// For each event, collect all PDG codes from Ntuple2.
// Accept the event if and only if it has at least one of each: 211, -211, 111.
std::map<int, std::vector<int>> event_products;
// ... fill from Ntuple2 ...
for (const auto& kv : event_products) {
    bool has_pip = false, has_pim = false, has_pi0 = false;
    for (int code : kv.second) { ... }
    if (has_pip && has_pim && has_pi0)
        selected_events.push_back(kv.first);
}
```

Only hits from `selected_events` are loaded from Ntuple3. This is a
truth-level selection — it uses MC decay product information, not detector
hits.

---

## 4. Stage 2 — Track Fitting (TGraphErrors pol1)

### 4.1 What goes into the fit and what doesn't

For each pion (π⁺ and π⁻ separately), the hits are split by particle PDG code
from Ntuple3's `PDGEncoding` branch. The following hits are **included** in the
track fit:

| Detector | IDs        | Contributes to | σ used in fit                          |
|----------|------------|----------------|----------------------------------------|
| T1       | 1–488      | x(z) fit only  | `STRAW_HALF_WIDTH` = 0.40 cm           |
| T2       | 489–976    | y(z) fit only  | `STRAW_HALF_WIDTH` = 0.40 cm           |
| T3       | 977–1464   | x(z) and y(z)  | `STRAW_HALF_WIDTH × √2` = 0.566 cm    |
| T4       | 1465–1952  | x(z) and y(z)  | `STRAW_HALF_WIDTH × √2` = 0.566 cm    |
| FRI W1   | 2001–2026  | x(z) fit only  | `fri_half_width(strip_i)` = 1.5 or 3.0 cm |
| FRI W2   | 2027–2052  | y(z) fit only  | `fri_half_width(strip_i)` = 1.5 or 3.0 cm |

The following hits are **excluded** from the track fit:

| Detector | IDs        | Reason for exclusion                                           |
|----------|------------|----------------------------------------------------------------|
| PIZZA    | 1953–2000  | Position used for timing/velocity only; not hit in track fit   |
| TOF bars | 2053–2070  | Bar centre has ±6 cm X and ±10 cm Y uncertainty; at ~460 cm   |
|          |            | lever arm from front trackers this strongly biases slope,      |
|          |            | degrading backward vertex extrapolation                        |

The T3/T4 stereo hit values `hitX` and `hitY` are both the **raw Geant4
position**. Because the straw is rotated ±45°, the measurement precision 
projects onto both X and Y in the lab frame, each with σ = σ_straw × √2.
This is an approximation — the true constraint is along the perpendicular to
the straw axis in the rotated frame.

### 4.2 x(z) and y(z) fits

Two independent `TGraphErrors` objects are built:

```
gx:  data points (z_i, x_i) with y-error σ_{x,i}
gy:  data points (z_i, y_i) with y-error σ_{y,i}
```

where `z_i = h.z`, `x_i = h.x`, `y_i = h.y` are all **raw Geant4 values**
directly from Ntuple3.

Each is fitted with a first-order polynomial (ROOT `"pol1"`):

$$x(z) = b_x + a_x \cdot z$$
$$y(z) = b_y + a_y \cdot z$$

ROOT minimises the chi-squared:

$$\chi^2_x = \sum_i \frac{(x_i - b_x - a_x z_i)^2}{\sigma_{x,i}^2}$$

The fit is run with the `"Q"` option (quiet, no printout). The four parameters
$(a_x, b_x, a_y, b_y)$ encode the complete 3D track.

### 4.3 What the fit output means

The fit result means the track satisfies:

- At `z = 0`: the track passes through point `(b_x, b_y, 0)` — the **reference origin**
- The track direction in the x-plane has slope `a_x = dx/dz`
- The track direction in the y-plane has slope `a_y = dy/dz`

The **3D direction vector** is:

$$\hat{d} = \text{Unit}(a_x,\; a_y,\; 1)$$

Note the z-component of the direction raw vector is fixed to 1 before
normalisation. This means the parametrisation implicitly assumes the track
moves in the +z direction (which is valid for all pion tracks after kaon
decay).

A point on the track at a given z is:

$$\mathbf{r}(z) = (b_x + a_x z,\; b_y + a_y z,\; z)$$

Or equivalently, using the unit direction vector, a point at parameter $t$ is:

$$\mathbf{r}(t) = \mathbf{orig} + t\,\hat{d} \qquad \text{where } \mathbf{orig} = (b_x,\; b_y,\; 0)$$

The parameter $t$ needed to reach a given $z$ is:

$$t = \frac{z - \mathbf{orig}.z}{\hat{d}.z} = \frac{z}{d_z}$$

since `orig.Z() = 0` by construction.

### 4.4 The TrackFit struct

```cpp
struct TrackFit {
    TVector3 dir;    // unit direction: (ax, ay, 1).Unit()
    TVector3 orig;   // reference point: (bx, by, 0)  — always z=0
    bool     valid;
};
```

The fit is rejected (`valid = false`) if:
- Fewer than 2 points are available in either gx or gy
- All points in either gx or gy lie at the same z (degenerate — slope undefined)
- ROOT fit returns non-finite parameters

---

## 5. Stage 3 — Decay Vertex via Point of Closest Approach (PoCA)

### 5.1 Two-line PoCA derivation

The decay vertex is estimated as the **midpoint of the shortest segment**
connecting the two reconstructed pion track lines. This is the standard
PoCA algorithm for two skew lines in 3D.

Given:
- Line 1:  $\mathbf{P}_1 + s\,\hat{v}_1$ where $\mathbf{P}_1 = \mathbf{orig}_{\pi^+}$, $\hat{v}_1 = \hat{d}_{\pi^+}$
- Line 2:  $\mathbf{P}_2 + t\,\hat{v}_2$ where $\mathbf{P}_2 = \mathbf{orig}_{\pi^-}$, $\hat{v}_2 = \hat{d}_{\pi^-}$

Define `w0 = P1 − P2`. The shortest connecting segment has parameters:

$$a = \hat{v}_1 \cdot \hat{v}_1, \quad b = \hat{v}_1 \cdot \hat{v}_2, \quad c = \hat{v}_2 \cdot \hat{v}_2$$
$$d = \hat{v}_1 \cdot \mathbf{w}_0, \quad e = \hat{v}_2 \cdot \mathbf{w}_0$$
$$\Delta = ac - b^2$$

$$s_c = \frac{b e - c d}{\Delta}, \qquad t_c = \frac{a e - b d}{\Delta}$$

The two closest points are:

$$\mathbf{Q}_1 = \mathbf{P}_1 + s_c\,\hat{v}_1 \qquad \mathbf{Q}_2 = \mathbf{P}_2 + t_c\,\hat{v}_2$$

The reconstructed decay vertex is their midpoint:

$$\mathbf{v}_{\text{reco}} = \frac{1}{2}(\mathbf{Q}_1 + \mathbf{Q}_2)$$

### 5.2 What is actually computed

Expanding the z-component explicitly:

$$z_{\text{vertex}} = \frac{1}{2}\left( [\mathbf{orig}_{\pi^+}]_z + s_c\,[\hat{d}_{\pi^+}]_z \;+\; [\mathbf{orig}_{\pi^-}]_z + t_c\,[\hat{d}_{\pi^-}]_z \right)$$

Since **both** `orig.Z() = 0` by construction:

$$\boxed{z_{\text{vertex}} = \frac{1}{2}\left( s_c\,d_z^{\pi^+} + t_c\,d_z^{\pi^-} \right)}$$

Where $d_z = \frac{1}{\sqrt{a_x^2 + a_y^2 + 1}}$ is the z-component of the unit
direction vector. The PoCA parameters $s_c$ and $t_c$ are the **arc-length along
the unit vector** from the `z=0` reference plane to each track's closest point.

### 5.3 The midpoint convention and its implications

The midpoint is exact only if the two tracks are perfectly consistent —
i.e., they intersect at a single point. For real reconstructed tracks the
two closest points $\mathbf{Q}_1$ and $\mathbf{Q}_2$ are not the same, and the midpoint
is a compromise.

If one track's closest point has $z_1$ and the other has $z_2$, then
$z_{\text{vertex}} = (z_1 + z_2)/2$. If both tracks are biased in the same
direction (e.g., both fitted slopes `a_x` are too shallow), both $z_1$ and
$z_2$ shift in the same direction, and the midpoint shifts accordingly.

### 5.4 Why small slope errors produce large z errors

The fitted track must be extrapolated **backward** from the detector region
to the decay vertex region. The detector planes span approximately:

| Region         | z range           |
|----------------|-------------------|
| Front trackers | 240–270 cm        |
| Rear trackers  | 680–700 cm        |
| Decay vertices | typically 50–300 cm |

The track is fitted to hits at z ≈ 240–700 cm and extrapolated to z ≈ 50–300 cm,
a backward extrapolation of up to ~200 cm from the nearest hit.

A slope error of $\delta a_x$ in the x-fit results in a lateral displacement at
the vertex:

$$\delta x_{\text{vertex}} = \delta a_x \cdot (z_{\text{vertex}} - z_{\text{ref}})$$

For z-fits (the PoCA parameter), the error in the vertex z comes from both tracks'
slope errors combining in the PoCA:

$$\delta z_{\text{vertex}} \sim \frac{\delta a_x \cdot \Delta x_{\text{sep}}}{\Delta a_x}$$

where $\Delta x_{\text{sep}}$ is the transverse separation between the two tracks
at the vertex and $\Delta a_x$ is the opening angle in x. At small opening angles
(boosted kaon, high momentum), even a small absolute slope error causes a large
fractional shift in the PoCA z parameter.

> **Implication for the known bias:** If `a_x` and `a_y` are systematically
> underestimated (tracks fitted as more forward/parallel than they truly are),
> the PoCA will place the vertex further downstream (+z) than the true vertex.
> The extrapolated track pair appears to converge later — pushed toward larger z.
> This is precisely the kind of bias that would be produced if the high-z tracker
> hits (T3, T4 at z ≈ 680–690 cm) carry too much weight relative to the
> front trackers.

---

## 6. Stage 4 — PIZZA Hit Position

The PIZZA detector provides the **TOF start** — the time at which the pion
passes through the pizza-shaped tracking detector near the beam axis.

**What is stored from Ntuple3:**

```cpp
// For each PIZZA hit (deviceID 1953–2000) for the pion:
double st = randGen.Gaus(hitT, 0.0015);  // smear time by 1.5 ps
if (!ev.has_pip_pizza || st < ev.pip_pizza_time) {
    ev.pip_pizza_time = st;       // earliest smeared hit time (ns)
    ev.pip_pizza_x    = hitX;     // raw Geant4 position (cm)
    ev.pip_pizza_y    = hitY;     // raw Geant4 position (cm)
    ev.pip_pizza_z    = hitZ;     // raw Geant4 z position (cm)
}
```

Only the **earliest** PIZZA hit (smallest smeared time) is kept per pion.

**What is then USED in reconstruction:**

The raw PIZZA x/y position is **discarded**. Instead, the PIZZA position is
recomputed by evaluating the fitted track at `z = pip_pizza_z`:

```cpp
TVector3 pip_pizza_pos = eval_track_at_z(fit_pip, pip_pizza_z);
```

Where `eval_track_at_z` evaluates:

$$\mathbf{r}_{\text{PIZZA}} = \mathbf{orig} + t_{\text{PIZZA}}\,\hat{d}$$

$$t_{\text{PIZZA}} = \frac{z_{\text{PIZZA}} - \mathbf{orig}.z}{\hat{d}.z} = \frac{z_{\text{PIZZA}}}{d_z}$$

$$\mathbf{r}_{\text{PIZZA}} = (b_x + a_x \cdot z_{\text{PIZZA}},\; b_y + a_y \cdot z_{\text{PIZZA}},\; z_{\text{PIZZA}})$$

The only PIZZA quantity that survives into the momentum calculation is:
- `pip_pizza_z` — the raw Geant4 z position of the PIZZA hit (defines which z-plane to evaluate the track at)
- `pip_pizza_time` — the smeared PIZZA hit time (ns)

> **Note:** `pip_pizza_z` is the raw Geant4 hit z — it reflects where along the z axis the
> pion actually hit the PIZZA. For angled PIZZA geometry (conical section) this z is
> not necessarily flat. However it is only used as an argument to `eval_track_at_z`,
> so any z variation in the PIZZA itself is captured.

---

## 7. Stage 5 — TOF Hit Position

The TOF wall provides the **TOF stop** — the time at which the pion reaches the
far detector wall.

**What is stored from Ntuple3:**

```cpp
// For each TOF hit (deviceID 2053–2070) for the pion:
double st = randGen.Gaus(hitT, 0.0015);  // smear time by 1.5 ps
if (!ev.has_pip_tof || st < ev.pip_tof_time) {
    ev.pip_tof_time     = st;                          // earliest smeared time (ns)
    ev.pip_tof_deviceID = id;                          // which bar fired
    ev.pip_tof_x = randGen.Gaus(tof_bar_x_centre(id), TOF_BAR_HALF_WIDTH);
    //  ^^^--- x is NOT taken from hitX; it is sampled from bar geometry
    ev.pip_tof_y = randGen.Gaus(tof_bar_y_centre(id), TOF_Y_UNCERTAINTY);
    //  ^^^--- y is NOT taken from hitY; sampled with large 10 cm sigma
    ev.pip_tof_z = z_tof;
    //  ^^^--- z is NOT taken from hitZ; it comes from the E1 filename tag
}
```

The TOF hit position is **not** from Geant4 coordinates. It is:

| Component | Source                                                         | Uncertainty |
|-----------|----------------------------------------------------------------|-------------|
| x         | `Gaus(tof_bar_x_centre(deviceID), 6.0 cm)`                    | ±6.0 cm     |
| y         | `Gaus(tof_bar_y_centre(deviceID), 10.0 cm)`                   | ±10.0 cm    |
| z         | `z_tof` — parsed from `E1-<value>` in the **filename**        | 0 cm (exact)|

`tof_bar_x_centre(id)` is computed as:
$$x_{\text{bar}} = -96.0 + (id - 2053) \times 12.0 \qquad \text{(cm)}$$

`tof_bar_y_centre(id)` is:
- Bar index 8 (id=2061): $y = -72$ cm
- Bar index 9 (id=2062): $y = +72$ cm
- All others: $y = 0$ cm

`z_tof` is parsed from the filename with:
```cpp
size_t pos = fname.find("E1-");
z_tof = stod(fname.substr(pos + 3, ...));   // e.g. "E1-700" → z_tof = 700.0 cm
```

> **Critical note:** The TOF z position is an exact value from the filename.
> It is the **nominal** detector face position, not the actual Geant4 hit z.
> If the TOF bars are not perfectly at `z_tof` in the simulation, this introduces
> a fixed offset in the TOF z used for velocity calculation.

---

## 8. Stage 6 — Pion Velocity

For each pion, the velocity is computed from the path length between PIZZA and TOF
divided by the measured transit time:

$$v_\pi = \frac{L_{\pi \to \text{TOF}}}{t_{\text{TOF}} - t_{\text{PIZZA}}}$$

In code:
```cpp
double pip_track_cm = (pip_tof_pos - pip_pizza_pos).Mag();   // |TOF_pos - PIZZA_pos| in cm
double pip_dt_ns    = pip_tof_time - pip_pizza_time;         // Δt in ns

double pip_v = (pip_track_cm * 1e-2) / (pip_dt_ns * 1e-9);  // converted to m/s
```

Where:
- `pip_tof_pos` = `(pip_tof_x, pip_tof_y, pip_tof_z)` — the TOF position from Stage 5
- `pip_pizza_pos` = result of `eval_track_at_z(fit_pip, pip_pizza_z)` — from Stage 4
- `pip_tof_time` = smeared Geant4 TOF hit time (ns)
- `pip_pizza_time` = smeared Geant4 PIZZA hit time (ns)

The path length is the **3D straight-line distance** from the track-extrapolated
PIZZA position to the smeared TOF bar position. It is **not** the arc length along
the track itself, but since pion tracks are nearly straight (negligible magnetic
field curvature in this geometry), this is a valid approximation.

> **Source of each quantity in the velocity:**
>
> | Quantity | Comes from |
> |---|---|
> | `pip_tof_pos.x` | `Gaus(bar_x_centre, 6 cm)` sampled at hit-collection time |
> | `pip_tof_pos.y` | `Gaus(bar_y_centre, 10 cm)` sampled at hit-collection time |
> | `pip_tof_pos.z` | `z_tof` from filename E1 tag |
> | `pip_pizza_pos.x` | `bx + ax * z_pizza` — from track fit |
> | `pip_pizza_pos.y` | `by + ay * z_pizza` — from track fit |
> | `pip_pizza_pos.z` | `z_pizza` — raw Geant4 PIZZA hit z |
> | `pip_tof_time` | `Gaus(Geant4 hitT, 1.5 ps)` — smeared Geant4 time |
> | `pip_pizza_time` | `Gaus(Geant4 hitT, 1.5 ps)` — smeared Geant4 time |

---

## 9. Stage 7 — Kaon Decay Time (Back-Propagation)

The kaon decay time is estimated by back-propagating each pion from its
PIZZA hit position to the reconstructed decay vertex, using the measured
pion velocity from Stage 6.

For π⁺:

$$t_{\text{decay}}^{\pi^+} = t_{\text{PIZZA}}^{\pi^+} - \frac{d(\mathbf{r}_{\text{PIZZA}}^{\pi^+},\; \mathbf{v}_{\text{reco}})}{v_{\pi^+}}$$

In code:
```cpp
double pip_path_cm    = (pip_pizza_pos - decay_vertex).Mag();
double pip_decay_time = pip_pizza_time * 1e-9
                        - (pip_path_cm * 1e-2 / pip_v);     // in seconds
```

And for π⁻ identically. The kaon decay time is the average:

$$t_{\text{decay}} = \frac{1}{2}\left(t_{\text{decay}}^{\pi^+} + t_{\text{decay}}^{\pi^-}\right)$$

> **Every quantity in this equation and its source:**
>
> | Symbol | Code variable | Source |
> |---|---|---|
> | $t_{\text{PIZZA}}^{\pi^+}$ | `pip_pizza_time * 1e-9` | Smeared Geant4 hit time, converted to seconds |
> | $\mathbf{r}_{\text{PIZZA}}^{\pi^+}$ | `pip_pizza_pos` | Track fit evaluated at raw Geant4 PIZZA z |
> | $\mathbf{v}_{\text{reco}}$ | `decay_vertex` | PoCA midpoint of two fitted tracks |
> | $\|\mathbf{r}_{\text{PIZZA}} - \mathbf{v}_{\text{reco}}\|$ | `pip_path_cm` | 3D Euclidean distance (cm) |
> | $v_{\pi^+}$ | `pip_v` | Velocity from Stage 6 (m/s) |

### The decay time depends critically on z_vertex

The path $d(\mathbf{r}_{\text{PIZZA}}, \mathbf{v}_{\text{reco}})$ is dominated by its z-component
because the PIZZA is at z ≈ 210–230 cm and the vertex is at z ≈ 50–250 cm.

If $z_{\text{vertex}}$ is overestimated (as identified), then:
- $\mathbf{r}_{\text{PIZZA}} - \mathbf{v}_{\text{reco}}$ has a **smaller z-component** than it should
- The 3D path from PIZZA to vertex is **underestimated**
- The back-propagation subtracts **less time** than it should
- $t_{\text{decay}}$ is therefore **overestimated** (the kaon appears to decay later)
- With a later (larger) $t_{\text{decay}}$ but the same flight path, the kaon appears to travel more slowly
- This directly produces **lower reco momentum**

This is the primary bias pathway from vertex z error to momentum underestimation.

---

## 10. Stage 8 — Kaon Momentum

The kaon is assumed to have been produced at rest at the origin `(0, 0, 0)`:

```cpp
TVector3 kaon_prod(0., 0., 0.);   // HARD-CODED: kaon production point
```

The kaon's flight path to the decay vertex:

$$L_K = |\mathbf{v}_{\text{reco}} - (0,0,0)| = |\mathbf{v}_{\text{reco}}|$$

The kaon's average velocity:

$$v_K = \frac{L_K}{t_{\text{decay}}}$$

The kaon's relativistic β:

$$\beta_K = \frac{v_K}{c}$$

(clamped to 0.9999 if $\beta_K \geq 1$)

The Lorentz factor:

$$\gamma_K = \frac{1}{\sqrt{1 - \beta_K^2}}$$

The reconstructed kaon momentum:

$$\boxed{p_K = \gamma_K \cdot m_K \cdot \beta_K \cdot c}$$

where $m_K = 0.497611\ \text{GeV}/c^2$ (K_L PDG mass).

> **Every quantity and its source:**
>
> | Symbol | Code variable | Source |
> |---|---|---|
> | $L_K$ | `flight_cm` | $|\mathbf{v}_{\text{reco}}|$ — the PoCA vertex from Stage 3 |
> | $t_{\text{decay}}$ | `kaon_decay_time` | Average back-propagated decay time from Stage 7, in seconds |
> | $v_K$ | `kaon_v` | `flight_cm * 1e-2 / kaon_decay_time` (m/s) |
> | $c$ | `2.99792458e8` | Hard-coded speed of light (m/s) |
> | $m_K$ | `0.497611` | Hard-coded K_L PDG mass (GeV/c²) |

---

## 11. Error Propagation Summary and Bias Pathways

### Primary bias chain (vertex z overestimation → momentum underestimation)

```
Track slope a_x or a_y is biased
    → PoCA parameter sc or tc shifts
        → z_vertex = (sc * dz_pip + tc * dz_pim) / 2  is  OVERESTIMATED
            → flight distance PIZZA→vertex  is  UNDERESTIMATED
                → back-propagated Δt removed from t_PIZZA  is  too small
                    → t_decay  is  OVERESTIMATED
                        → v_K = flight_K / t_decay  is  UNDERESTIMATED
                            → p_K = γ·m·β  is  UNDERESTIMATED  ✓
```

### Secondary effects on momentum magnitude

| Error source | Effect on $p_K$ |
|---|---|
| z_vertex overestimated | $p_K$ underestimated (dominant, as above) |
| $v_\pi$ overestimated (path too long or $\Delta t$ too short) | More time subtracted in back-prop → $t_{\text{decay}}$ decreases → $p_K$ increases |
| TOF z wrong (from filename E1 tag) | Affects both $L_{\pi}$ (velocity) and cannot compensate itself |
| PIZZA z wrong (raw Geant4 hitZ) | Shifts PIZZA position used in both velocity and back-prop |
| Kaon production assumed at origin | If kaon produced off-axis, $L_K$ is biased |

### Key places to instrument for debugging

1. **PoCA closest-point separation** — `(Q1 - Q2).Mag()`: large separations indicate poor vertex resolution or track mismatch
2. **Per-pion decay time consistency** — `pip_decay_time` vs `pim_decay_time` before averaging: large discrepancies indicate inconsistent velocity or vertex z
3. **True vs reco z_vertex** — compare `decay_vertex.Z()` vs `truth_map[ev].vz`: if reco z > true z systematically, this confirms the slope bias
4. **Track slope values** — print `ax, ay` for each pion and compare to `(true_vertex - PIZZA_hit).X / (true_vertex - PIZZA_hit).Z` to check for slope bias
5. **Track evaluation at true vertex z** — evaluate `fit_pip` at `truth_map[ev].vz` and compare (x,y) position with the truth vertex (x,y) to isolate whether the error is in z or in the transverse position

---

*Generated from `KLong_save_vectors.C` — reconstruction version with TOF excluded from track fit.*
