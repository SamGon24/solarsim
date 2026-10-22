# SolarSim: physics and algorithms

This document explains the equations behind every module, in the order data flows through the
program: units → data → initial conditions → forces → integration → analysis → rendering.

---

## 1. Units and time scales

| Quantity | Unit |
|---|---|
| length | astronomical unit, 1 AU = 149 597 870.700 km (IAU 2012, exact) |
| time | day = 86 400 SI seconds, on the **TDB** time scale (≈ TT, differences < 2 ms) |
| mass | gravitational parameter **GM** in AU³/day² |

GM is used everywhere instead of G and M separately. GM_☉ is known to about 10 significant digits from
planetary radar ranging, while G alone is only known to about 5. Conversion:
`GM[AU³/d²] = GM[km³/s²] · 86400² / AU_km³`, so GM_☉ = 2.959122×10⁻⁴ AU³/d² = k², with k the Gaussian
gravitational constant 0.01720209895.

**Civil time.** User dates are UTC. Internally `TT = UTC + ΔT`, where ΔT is the Espenak–Meeus
polynomial (≈ 69 s in 2026), blended with the observed IERS trend after 2015 (`astro/JulianDate.cpp`).
Julian dates follow Meeus ch. 7, using the Julian calendar before 1582-10-15 and the Gregorian calendar
after.

---

## 2. Data sources (`astro/BodyCatalog.cpp`)

| Data | Source |
|---|---|
| GM of Sun, planets and Moon | JPL DE440 (Park et al. 2021) |
| radii, rotation periods, pole RA/Dec | IAU WGCCRE 2015 (Archinal et al. 2018) |
| planetary orbits | E. M. Standish, *Keplerian Elements for Approximate Positions of the Major Planets*, Table 1 (1800–2050): a, e, I, L, ϖ, Ω and their rates per century |
| Moon | ELP-2000/82 truncated as in Meeus ch. 47 (Tables 47.A and 47.B, about 120 terms, plus additive Venus, Jupiter and flattening terms) |
| Galilean moons | Meeus ch. 44 mean longitudes, which preserve the Laplace resonance λ_Io − 3λ_Eu + 2λ_Ga = 180° |
| bright stars | J2000 RA/Dec, V magnitude, B−V colour index (Yale Bright Star Catalogue, rounded) |

Planet GM values exclude their large moons. A moon that is not simulated has its GM *folded back* into
the parent, so the total mass of every planetary system is always conserved.

---

## 3. Orbital mechanics building blocks (`astro/Kepler.cpp`, `astro/OrbitalElements.cpp`)

**Kepler's equation.**

- Elliptic: `M = E − e sin E`, solved by Halley's method (cubic convergence) with Danby's starting
  value `E₀ = M + 0.85 e sign(sin M)`. It converges for 0 ≤ e < 1, including e = 0.999.
- Hyperbolic: `M = e sinh H − H`, solved by Halley's method with damping.
- Parabolic: Barker's equation, solved in closed form, `tan(ν/2) = Y − 1/Y` with
  `Y = ∛(W/2 + √(W²/4+1))`.

**One element set for all conics.** Elements are stored as (q, e, i, Ω, ω, t_p): perihelion distance and
time of perihelion. This works for ellipses, parabolas and hyperbolas alike. The semi-major axis
`a = q/(1−e)` and mean anomaly `M = n(t − t_p)` are derived from them.

**State ⇄ elements.** These use the eccentricity vector `e = ((v² − μ/r) r − (r·v) v)/μ` and the node
vector `n = ẑ × h`. The singular cases are handled explicitly:

- circular orbits: argument of latitude instead of ν
- equatorial orbits: longitude of periapsis instead of ω
- retrograde equatorial orbits: angle signs flipped

**Universal-variable propagator** (Battin and Vallado). The Stumpff functions

```
c₂(ψ) = (1 − cos √ψ)/ψ,   c₃(ψ) = (√ψ − sin √ψ)/ψ^{3/2}    (series near 0, cosh/sinh for ψ < 0)
```

let one Newton iteration on the universal anomaly χ cover every conic. The new state then follows from
the Lagrange coefficients f, g, ḟ, ġ. This propagator is used by the Wisdom–Holman integrator and by
the tests.

**Lambert's problem** (`astro/Lambert.cpp`). Given r₁, r₂ and a time of flight Δt, the solver finds
v₁ and v₂. It uses the universal-variable formulation with bisection on ψ (Bate–Mueller–White, Vallado
algorithm 58), which is slow but very robust for single-revolution transfers.

---

## 4. Reference frames (`astro/Coordinates.cpp`)

The inertial frame is the **ecliptic and mean equinox of J2000**: x points to the vernal equinox and z
to the north ecliptic pole.

- **Equatorial ⇄ ecliptic:** rotation about x by ε₀ = 84381.406″ (IAU 2006).
- **Precession** (IAU 1976): `P = R_z(z) R_y(−θ) R_z(ζ)`, with
  ζ = 2306.2181″T + 0.30188T² + 0.017998T³ (and similar for z and θ). This reproduces Meeus example
  21.b to 10⁻⁵°, which is checked in the tests.
- **Nutation:** abridged IAU 1980 (the four largest terms, about 0.5″).
- **Apparent solar longitude** = geometric longitude + general precession p_A + Δψ − 20.4898″/R
  (annual aberration).
- **Horizontal coordinates:** GMST (IAU 1982) → local sidereal time → hour angle → alt/az, with Bennett
  refraction and topocentric parallax (the observer's position vector on the Earth's surface).

---

## 5. Initial conditions (`physics/SystemBuilder.cpp`, `physics/Calibration.cpp`)

1. Evaluate each orbit model at the epoch, relative to the parent body.
   - Standish: L, ϖ, Ω, I, e, a at time T → osculating conic with μ = GM_☉ + GM_system.
2. **Barycentric subsystems.** Standish's "EM Bary" row is the Earth–Moon barycentre, and planetary
   mean orbits likewise describe system barycentres. A planet with simulated moons is therefore displaced
   by `−Σ GMₖ rₖ / GM_system` and the moons are placed around it.
3. **Solar-System barycentre:** subtract the centre-of-mass position and velocity, so that P_total = 0.
4. **Differential correction.** Mean elements and truncated series yield velocities that are
   inconsistent with the dynamics at the 10⁻⁵ level. For a near-circular orbit,
   `dn/n = −3/2 · da/a = −3 dv/v`, so such errors make a body drift linearly in longitude. The
   calibration step
   - integrates a reduced system (Sun, planets, Moon) for 1 year (Moon) or 30 years (planets),
   - fits the longitude difference to the reference theory with least squares, giving the slope dn,
   - multiplies the body's velocity by `1 + dn/(3n)` (a body that runs ahead needs *more* energy),
   - repeats this twice and applies the factors to the full system.

   For the Moon this changes the error in lunar-phase times from about 110 min to about 1 min.

**Synthetic populations** are massless test particles drawn from realistic distributions:

- **Main belt:** a ∈ [2.1, 3.35] AU, with the Kirkwood gaps (4:1, 3:1, 5:2, 7:3, 2:1) left empty;
  e and i are Rayleigh-distributed (σ_e = 0.09, σ_i = 7.5°).
- **Trojans:** a ≈ a_J, λ = λ_J ± 60° + N(0, 9°).
- **Kuiper belt:** classical objects (a ∈ [42, 47.8] AU, cold and hot inclination components) plus
  plutinos in the 3:2 resonance, with resonant angle `φ = 3λ − 2λ_N − ϖ` librating about 180°.

**Intruder star:** a hyperbolic orbit with `e = 1 + q v∞²/μ`.

---

## 6. Forces (`physics/Forces.cpp`)

**Newtonian gravity:**

```
a_i = Σ_j GM_j (r_j − r_i) / (|r_j − r_i|² + ε²)^{3/2}
```

Massive pairs are summed once (Newton's third law). Test particles only feel the massive bodies, which
costs O(N_test · N_massive), and that loop runs in parallel on a persistent thread pool. Pair exclusions
(source → target) implement patched-conic flights, and a "skip source" option implements the
Wisdom–Holman interaction Hamiltonian.

**General relativity (1PN).** This is the Schwarzschild test-particle limit of the
Einstein–Infeld–Hoffmann equations in harmonic coordinates, relative to the Sun:

```
a_GR = GM / (c² r³) · [ (4 GM/r − v²) r + 4 (r·v) v ]
```

For Mercury this gives the classic `Δϖ = 6πGM / (c² a (1−e²))` per orbit, which is 42.98″ per century.
Because the force depends on velocity, the Newtonian energy used in diagnostics is no longer exactly
conserved when `--gr` is on, which is expected.

---

## 7. Integrators (`physics/Integrators.cpp`)

| Method | Scheme (step h) |
|---|---|
| Euler | x += h v; v += h a(x) |
| Symplectic Euler | v += h a(x); x += h v |
| Leapfrog (DKD) | x += h/2 v; v += h a(x); x += h/2 v |
| RK4 | classical 4-stage on y = (x, v) |
| Yoshida 4 | drift c₁ kick d₁ drift c₂ kick d₂ drift c₃ kick d₃ drift c₄, with w₁ = 1/(2 − 2^{1/3}), w₀ = −2^{1/3}/(2 − 2^{1/3}), c₁ = c₄ = w₁/2, c₂ = c₃ = (w₀+w₁)/2, d₁ = d₃ = w₁, d₂ = w₀ |
| Dormand–Prince 5(4) | 7 stages, embedded 4th-order error estimate, max-norm per body with atol/rtol, step factor `0.9·err^{−1/5}` clamped to [0.2, 5] |
| Wisdom–Holman (DH) | `H = H_Kepler + H_interaction + H_Sun` in democratic heliocentric coordinates (Duncan, Levison & Lee 1998). Each step: Sun-drift(h/2) · kick(h/2) · Kepler(h) · kick(h/2) · Sun-drift(h/2) |

**Why symplectic matters.** A symplectic map conserves a nearby "shadow Hamiltonian", so the energy error
oscillates but stays bounded forever, while RK4 and DOPRI5 drift secularly. Run `solarsim compare` to see
this: leapfrog stays at about 3×10⁻⁸ and Yoshida at about 2×10⁻¹¹, both bounded, while RK4 grows
steadily. Symplectic methods also have their own small precession error, which is why the precession
experiment uses dt = 0.04 d. With dt = 0.25 d, Yoshida adds a spurious −74″/cy to Mercury.

**Automatic step.** dt = (shortest circular period at any body's periapsis) / 150. That gives 0.16 d with
the Moon, 0.012 d with the Galilean moons, and about 0.4 d without the Moon (set by Mercury's perihelion).

---

## 8. Analysis

### Event finding (`analysis/Events.cpp`)
The search samples the trajectory every 0.5 day and keeps snapshots. A root of a geometric function (for
example `λ_planet − λ_☉ − 180°` for an opposition) is bracketed between samples and refined with **Brent's
method**, where every function evaluation re-integrates a copy of the snapshot. Extrema (greatest
elongation, perihelion) are refined with **golden-section search** the same way. Angle functions are
wrapped to (−180°, 180°], and jumps larger than 90° are rejected as wrap-arounds.

### Observations (`analysis/Observation.cpp`)
- light-time iteration `τ = |x(t−τ) − o|/c` with linear back-propagation,
- elongation, phase angle, illuminated fraction `(1 + cos i)/2`,
- magnitudes from Meeus ch. 41 (Saturn includes the ring-tilt term −2.60 sin|B| + 1.25 sin²B, with B from
  the IAU pole),
- angular diameter `2 asin(R/Δ)`.

### Lagrange points (`analysis/Lagrange.cpp`)
In the rotating frame with the barycentre at the origin, the collinear points solve

```
x − (1−μ)(x+μ)/|x+μ|³ − μ(x−1+μ)/|x−1+μ|³ = 0
```

which is solved with Brent's method on three brackets. L4 and L5 sit at `(½ − μ, ±√3/2)`.

- Routh's criterion: L4/L5 are stable if `μ < ½(1 − √(23/27)) = 0.0385209`.
- Hill sphere: `a(1−e)∛(m/3M)`. Laplace sphere of influence: `a (m/M)^{2/5}`.
- The stability experiment releases massless probes with co-rotating velocities `v_B + ω × r`.

### Mission design (`analysis/Transfer.cpp`)

1. **Hohmann:** `a_t = (r₁+r₂)/2`, Δv₁ = √(μ/r₁)(√(2r₂/(r₁+r₂)) − 1), TOF = π√(a_t³/μ), required phase
   angle `φ = π − n₂·TOF`.
2. **Porkchop:** Lambert arcs on a grid of launch date × TOF, using planet states from a daily N-body table
   with cubic Hermite interpolation. The cost is
   `Δv_dep + Δv_cap = [√(v∞² + 2μ/r_p) − √(μ/r_p)]_dep + [same]_arr` (Oberth effect from circular orbits).
3. **Single shooting:** solve `F(v₀) = r_sc(T) − r_target(T) − b = 0` with Newton iterations, using a 3×3
   finite-difference Jacobian (δv = 1.7 m/s) of the full N-body propagation. Planetary gravity on the
   probe is switched off here, which is the patched-conic assumption.
4. **B-plane:** the patched-conic impact parameter is `b = r_p √(1 + 2μ/(r_p v∞²))`. The target's
   gravity bends the approach long before its sphere of influence, so a **secant iteration on b**, run
   with full-gravity flights, adjusts it until the requested periapsis is reached. Real missions do the
   same job with trajectory-correction manoeuvres.
5. **DLA:** the declination of the departure asymptote v∞ in the equatorial frame.

### Precession (`analysis/Precession.cpp`)
The code samples the osculating heliocentric longitude of perihelion ϖ = Ω + ω every 2 days, unwraps it,
and fits a least-squares line, giving the slope in ″/century and its standard error.

---

## 9. Rendering

### Display mapping (`render/DisplayMapping.cpp`)
Heliocentric distances are compressed radially, with `r' = √r` (default), `ln(1+r/0.5)/ln 3`, or linear.
Every mode keeps 1 AU at 1. Satellite offsets are expanded logarithmically,
`D = R_p (2.6 + 2.2 log₂(d/d_min))`, so moons sit outside the planet's (exaggerated) disc. Radii are
exaggerated as `0.022 (R/R⊕)^{0.45}`, and `--true-scale --scale linear` restores physical proportions.

### Ray casting (`render/SceneRenderer.cpp`)
- **Spheres:** per pixel, intersect the ray `o + t d` with the sphere (numerically stable perpendicular-
  distance form), with coverage-based anti-aliasing at the limb.
- **Surface coordinates:** the normal is rotated into the body frame (spin axis from the IAU pole,
  rotation angle W = W₀ + 2π t/P), then fed to procedural 3-D value noise or fBm (Earth's continents and
  clouds, gas-giant zonal bands with a Great Red Spot at −22° latitude, ice ridges, Io's volcanoes).
- **Shading:**
  - Lambert lighting from the Sun, with a softened terminator for atmospheric bodies
  - Blinn–Phong ocean glint
  - Fresnel-like atmospheric rim
  - solar limb darkening `I(μ) = 1 − 0.6(1 − μ)` (Eddington approximation)
- **Rings:** ray–plane intersection in the equatorial plane, with an opacity profile using the real
  radii of Saturn's C ring, B ring, Cassini division (117 580–122 170 km), A ring and Encke gap. The
  planet's shadow on the rings and the rings' shadow on the planet are both found with secondary rays.
- **Comets:** an ion tail along the anti-solar direction and a curved dust tail (Bézier between the
  anti-solar and the lagging direction). Activity scales as 1/r², and the effect switches off beyond 5 AU.
- **Post-processing:** bright-pass, ¼-resolution separable Gaussian blur (two radii), additive bloom,
  ACES filmic tone mapping (Narkowicz fit), sRGB encoding.

### Terminal output (`render/TerminalRenderer.cpp`)
Each character cell is two pixels: '▀' with the top pixel as foreground and the bottom pixel as
background colour. Colour codes are only emitted when they change, and the screen is redrawn with a
single cursor-home write to avoid flicker.

### File formats (`io/`)
- **PNG:** zlib stream with LZ77 (32 KiB window, hash chains) and fixed-Huffman DEFLATE blocks, adaptive
  per-row filters (None/Sub/Up/Average/Paeth), CRC-32 and Adler-32.
- **GIF89a:** median-cut palette per frame (15-bit histogram, split the box with the largest
  range·√population at the weighted median), optional Floyd–Steinberg dithering (serpentine scan),
  variable-width LZW (9–12 bit codes, clear on a full dictionary), NETSCAPE2.0 looping.
- **SVG:** vector diagrams. **CSV:** RFC 4180 quoting.

---

## 10. References

- Standish, E. M. *Keplerian Elements for Approximate Positions of the Major Planets*, JPL SSD.
- Park, R. S. et al. (2021) *The JPL Planetary and Lunar Ephemerides DE440 and DE441*, AJ 161, 105.
- Archinal, B. A. et al. (2018) *Report of the IAU WG on Cartographic Coordinates and Rotational Elements: 2015*.
- Meeus, J. (1998) *Astronomical Algorithms*, 2nd ed., Willmann-Bell.
- Vallado, D. A. *Fundamentals of Astrodynamics and Applications*.
- Battin, R. H. *An Introduction to the Mathematics and Methods of Astrodynamics*.
- Yoshida, H. (1990) *Construction of higher order symplectic integrators*, Phys. Lett. A 150, 262.
- Dormand, J. R. & Prince, P. J. (1980) *A family of embedded Runge–Kutta formulae*, J. Comp. Appl. Math. 6, 19.
- Wisdom, J. & Holman, M. (1991) *Symplectic maps for the N-body problem*, AJ 102, 1528.
- Duncan, M., Levison, H. & Lee, M. H. (1998) *A multiple time step symplectic algorithm…*, AJ 116, 2067.
- Einstein, A. (1915) *Erklärung der Perihelbewegung des Merkur aus der allgemeinen Relativitätstheorie*.
- Heckbert, P. (1982) *Color image quantization for frame buffer display*, SIGGRAPH.
- RFC 1950 / 1951 (zlib / DEFLATE), PNG specification (ISO/IEC 15948), GIF89a specification.
