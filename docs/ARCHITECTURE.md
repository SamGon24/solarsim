# SolarSim: architecture

## Layering

```
                      ┌──────────────────────────────┐
                      │ main.cpp  (dispatch)          │
                      └──────────────┬───────────────┘
                      ┌──────────────▼───────────────┐
                      │ app/   one file per command   │  AppContext: options → objects
                      └─┬─────────┬─────────┬────────┘
          ┌─────────────▼──┐ ┌────▼─────┐ ┌─▼──────────────┐
          │ analysis/      │ │ render/  │ │ io/            │
          │ events, mission│ │ renderer │ │ PNG GIF SVG CSV│
          │ Lagrange, ...  │ │ terminal │ │ body files     │
          └───────┬────────┘ └────┬─────┘ └───────┬────────┘
                  │          ┌────▼─────┐         │
                  └─────────►│ physics/ │◄────────┘
                             │ N-body   │
                             └────┬─────┘
                             ┌────▼─────┐   ┌───────────┐
                             │ astro/   │   │ platform/ │ terminal & OS
                             └────┬─────┘   └───────────┘
                             ┌────▼─────┐
                             │ core/    │ math, parsing, threads
                             └──────────┘
```

Every layer depends only on the layers below it. `platform/` is the only place with OS-specific code
(`#ifdef _WIN32`).

All code lives in the `solsim` namespace, with sub-namespaces `math`, `str`, `log`, `constants`,
`astro`, `physics`, `analysis`, `render`, `io`, `platform` and `app`.

## Key types

| Type | File | Role |
|---|---|---|
| `Vector3<T>`, `Mat3` | core/Vec3.hpp, core/Mat3.hpp | templates for geometry |
| `NBodySystem` | physics/NBodySystem.hpp | structure-of-arrays storage: massive bodies, then massless named bodies, then particles |
| `ForceModel` (abstract) → `NewtonianGravity`, `PostNewtonianGR` | physics/Forces.hpp | add accelerations for arbitrary (x, v) arrays |
| `ForceSet` | physics/Forces.hpp | ordered list of force models (deep-copyable via `clone()`) |
| `Integrator` (abstract) → 7 implementations | physics/Integrators.hpp | `step(system, forces, dt)` returns the time actually advanced |
| `Simulation` | physics/Simulation.hpp | owns system + forces + integrator, the clock, invariants; copyable (snapshots) |
| `BodyDefinition`, `BodyCatalog` | astro/BodyCatalog.hpp | built-in data and orbit models |
| `SceneRenderer` | render/SceneRenderer.hpp | draws a `Simulation` through a `Camera` into a `Framebuffer` |
| `TerminalRenderer` | render/TerminalRenderer.hpp | `Image8` + labels → ANSI string |
| `GifEncoder`, `MedianCutQuantizer` | io/GifEncoder.hpp | animated GIF writer |
| `AppContext` | app/AppContext.hpp | translates CLI options into builder, simulation and render settings |

### Polymorphism and run-time type information

- `ForceModel` and `Integrator` are classic abstract base classes with virtual `clone()`.
  `Simulation`'s copy constructor clones them, which is how the event finder takes snapshots.
- `ForceSet::gravity()` uses `dynamic_cast` to find the Newtonian model (the Wisdom–Holman integrator and
  the mission planner need it).
- Errors are reported with exceptions (`UsageError` for bad input, `std::runtime_error` otherwise), which
  `main()` turns into exit codes 2 and 1.

## Data flow of a frame (view / render)

```
Simulation::advance(days)
   └─ Integrator::step(...)  ×n
        └─ ForceSet::compute → NewtonianGravity::accumulate (ThreadPool for particles)
                             → PostNewtonianGR::accumulate
SceneRenderer::record()         trail ring buffers (positions relative to the parent)
SceneRenderer::render()
   stars → bodies (ray cast, depth) → grid → orbits → trails → particles → rings → comet tails → glow
finishFrame()                   bloom + ACES tone map → Image8
drawOverlay() / TerminalRenderer::compose()
GifEncoder::addFrame() / writeImage() / writeStdout()
```

## Extending

**A new integrator.**
1. Derive from `physics::Integrator` and implement `name`, `displayName`, `order`, `isSymplectic`,
   `forceEvaluationsPerStep`, `step` and `clone`.
2. Register it in `createIntegrator()` and `integratorNames()` in `physics/Integrators.cpp`.
3. `compare`, `view` (`i` key) and `--integrator` pick it up automatically.

**A new force** (for example solar J2 or radiation pressure). Derive from `physics::ForceModel`, *add* to
`acc[i]` in `accumulate`, and insert it in `Simulation::configureForces()` behind an option.

**A new command.**
1. Write `int cmdXxx(AppContext&)` in `src/app/CmdXxx.cpp`.
2. Declare it in `AppContext.hpp`.
3. Register the name in `registerOptions()` (with usage text) and in the dispatch map in `main.cpp`.

**A new body.** Prefer a data file (`--bodies`). Built-ins go in `buildCatalog()` in
`astro/BodyCatalog.cpp`.

## Testing

`tests/` contains a tiny self-registering test framework (`TestFramework.hpp`) and 45 tests:

- **core:** vectors and matrices, root finders, least squares, RNG statistics, CLI and INI parsing,
  thread-pool coverage, Julian dates (Meeus examples 7.a, 7.b, 7.e), time-scale parsing
- **astro:** Kepler equations for every e, element round-trips for all conic types, the universal
  propagator, Lambert, frames, precession (Meeus 21.b), Standish Earth position, ELP Moon (Meeus 47.a)
- **physics:** **measured convergence order** of every integrator, bounded vs drifting energy, time
  reversibility, barycentric construction, mass folding, populations bound, GR precession (42.98″/cy),
  Lagrange points
- **io:** CRC-32 and Adler-32 check values; DEFLATE verified by an **independent inflater** written in the
  test; GIF LZW verified by an **independent decoder**; PNG chunk CRCs and full filter reversal;
  median-cut quality; body-file loader
- **render and analysis:** camera projection and ray consistency, mapping monotonicity, font, terminal
  output, Hohmann numbers, the June 2024 solstice and the 2024-04-08 eclipse new moon against published
  times, the Sun's apparent RA/Dec at J2000

Run them with `./build/solarsim_tests [filter]`, `make tests` or `(cd build && ctest)`.
