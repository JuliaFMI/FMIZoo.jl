This FMU is a modified version of BouncingBall from the Modelica Reference-FMUs repository.

Original License (2-Clause BSD):
Copyright (c) 2024, Modelica Association Project "FMI".
All rights reserved.

--------------------------------------------------------------------------------
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
--------------------------------------------------------------------------------

Modifications:
Modified to an LS-SA compatible version of BouncingBall. These modifications are provided under the MIT License (see FMIZoo.jl root).

# NBouncingBalls (LS-SA reference FMU)

`NUM_BALLS` independent bouncing balls in one FMI 3.0 FMU (ME + CS), implementing the proposed
layered standard for sensitivity analysis (LS-SA). One FMU per ball count (compile time constant).

## Model

Per ball `i` (value reference `100 + (i-1)*12 + offset`):

| offset | variable                     | equation / role                                   |
|-------:|------------------------------|---------------------------------------------------|
| 0      | `h<i>`                       | continuous state, output                          |
| 1      | `der(h<i>)`                  | `= v`                                             |
| 2      | `v<i>`                       | continuous state, output                          |
| 3      | `der(v<i>)`                  | `= g + k*time - drag*(1-e)*v*abs(v)`              |
| 4,5,6  | `e<i>`, `k<i>`, `drag_coefficient<i>` | tunable parameters                       |
| 7      | `bounce_count<i>` (Int32)    | discrete state (LS-SA: settable in Event Mode)    |
| 8,9    | `pre(h<i>)`, `pre(v<i>)`     | LS-SA `previous` variables                        |
| 10     | `error_h_start<i>`           | LS-SA error indicator `= -h`                      |
| 11     | `event_indicator_ground<i>`  | event indicator `= h + ground_flexibility*(1-e)`  |

Globals: `time` (0), `g` (1), `v_min` (2), `ground_flexibility` (3), `base_event_period` (4), `nextEventTime` (5).

Events, within one event iteration: `pre(x) := x` for all balls; at a **time event** all balls are
kicked (`v += 2`) and the next one is scheduled at
`time + base_event_period * (1 + 0.1 * mean(pre(h)))`, i.e. the next event time depends on the
state; a **state event** (`h <= ground` and `v < 0`) resets `h` to the ground and sets `v := -e*v`
(`v := 0` if `abs(v) < v_min`).

The `ModelStructure` declares exact `dependencies` (FMI 3.0 knowns: states, parameters, time) for all
outputs, state derivatives and event indicators: `der(h<i>)` on `v<i>`; `der(v<i>)` on `time`, `g`, `v<i>`,
`e<i>`, `k<i>`, `drag_coefficient<i>`; `event_indicator_ground<i>` on `ground_flexibility`, `h<i>`, `e<i>`;
the outputs on themselves. Initial unknowns have no `dependencies` attribute (may depend on all knowns).
Importers can use them for sparse (colored) finite differences.

## Coupled variant (NCoupledBouncingBalls)

`NCoupledBouncingBalls_<N>_K<K>`: the same balls, each coupled to its `K` nearest neighbors (by index) on
each side, `0 < |j - i| <= K` (built with `NUM_NEIGHBORS = K`; `K = 0` is exactly the model above):

* `der(v<i>)` additionally contains the force of hardening springs with dampers,
  `1/K * sum_j [ (c_i + c_j)/2 * phi(h_j - h_i) + 0.1 * (v_j - v_i) ]` with `phi(x) = x + x^3`.
  `c_i` is the tunable parameter `spring_stiffness<i>` (offset 12, i.e. value reference
  `100 + (i-1)*13 + offset` in this variant, start 1).
* Impact transfer: after all bounces of an event iteration, every bounced (or clamped) ball `i` passes
  `0.2/K` of its velocity change `dv_i = v+ - v-` to each of its neighbors (`v_j += 0.2/K * dv_i`).
  `dv_i > 0`, so no further bounce is caused.

The springs are nonlinear, so the coupling costs O(N*K) per evaluation for every importer, the same as the
number of nonzeros of the Jacobian. The `dependencies` of `der(v<i>)` list `h`, `v`, `c` of the neighbors
(dense for `K = N - 1`). LS-SA: the fast paths evaluate the coupling terms of the directional / adjoint
derivatives spring by spring, each spring once for both balls, exactly like `fmi3GetContinuousStateDerivatives`
evaluates the forces (O(N*K) for a whole Jacobian product, O(K) for a single ball; the generic path evaluates
them per unknown); derivatives the coupling does not enter (e.g. w.r.t. time) cost nothing extra. The event map
includes the impact transfer, evaluated from the list of balls that impacted (O(N + impacts * K)).

Available for N = 3, 10, 30, 100, 300, 1000 and K = 1, 8, N - 1 (dense); K is capped at N - 1
(e.g. `NCoupledBouncingBalls_3_K2` is both "K = 8" and "dense" for N = 3).

## LS-SA features

* `fmi3GetDirectionalDerivative` / `fmi3GetAdjointDerivative` provide analytic partial derivatives for
  **all** Float64 variables (every variable depends on itself): w.r.t. states, `time` and tunable
  parameters in every mode, and of event indicators.
* Event map (paper Sec. 4.2): `d x+ / d pre(x)` of the last event iteration (identity for balls not
  affected by the event). `d x+ / d p` (`e`, `ground_flexibility`) is returned in **Event Mode only**;
  in Continuous-Time Mode the outputs `h`, `v` do not depend on parameters.
* Next event time (paper Sec. 4.4): `d nextEventTime / d {pre(h<i>), time, base_event_period}` is
  cached when the time event is handled (zero if the last event iteration did not schedule a new
  time event).
* Spec strictness: tunable parameters can only be set in Event Mode (ME) / Step Mode (CS);
  derivatives w.r.t. non-Float64 variables return `fmi3Error`.

## Performance: optimized on both sides

This FMU is used to benchmark FMI 3.0 *with* LS-SA against the best workarounds in *plain* FMI 3.0.
Such a comparison is only meaningful if **both** sides are implemented efficiently, so the functions a
plain FMI 3.0 importer needs received the same attention as the LS-SA functions. All FMUs of a
platform are compiled with the same compiler and flags (x86_64-windows: MSVC 2022, `/O2`, static
runtime; x86_64-linux: GCC, `-O3`). The cost tables below were measured with the Windows binaries.

**Standard FMI 3.0 functions** (what finite-difference workarounds call):

| function | implementation | cost |
|---|---|---|
| `fmi3SetContinuousStates`, `fmi3GetContinuousStates` | one loop over the balls, no copies, no allocation | O(N) |
| `fmi3GetContinuousStateDerivatives`, `fmi3GetEventIndicators` | evaluated on demand in one loop, nothing cached or flagged dirty | O(N) |
| `fmi3GetFloat64`, `fmi3SetFloat64`, `fmi3GetInt32`, `fmi3SetInt32` | value references decoded arithmetically (no lookup table, no search), derived variables computed on demand | O(1) per value |
| `fmi3SetTime`, `fmi3EnterEventMode`, `fmi3EnterContinuousTimeMode` | flags only | O(1) |
| `fmi3UpdateDiscreteStates` | one pass over the balls | O(N) |
| `fmi3GetFMUState` | one `memcpy` into the FMU state passed by the importer (reused, no reallocation) | O(N) |
| `fmi3SetFMUState` | copies the model data and event flags only (no self-copies of solver buffers) | O(N) |

The FMU state is lean: **64 bytes per ball** (states, parameters, `pre` values, discrete state and one
byte of LS-SA event information). Derived quantities (state derivatives, event and error indicators)
are not stored. All LS-SA scratch data (buffers, plan cache) lives *outside* the FMU state, so the
snapshots the workarounds depend on are not inflated by LS-SA, and LS-SA bookkeeping during event
handling is limited to one byte per ball (see lazy event map below).

**LS-SA functions** (`fmi3GetDirectionalDerivative`, `fmi3GetAdjointDerivative`):

* *Analytic and sparse*: closed-form partial derivatives, only structural non-zeros are evaluated;
  cost O(nUnknowns + nKnowns) (the rows / columns of `nextEventTime`, `time`, `g`,
  `ground_flexibility` are inherently dense, O(N)).
* *Cached derivative plans*: every (unknowns, knowns) pair is analysed once (validation and pattern
  detection) and cached per instance (8 entries, keyed by the list contents via `memcmp`). Importers
  repeat the same pattern in every step; only the seed changes. Build option `LSSA_PLAN_CACHE=OFF`
  analyses every call instead (used to benchmark "cold" calls).
* *Streaming fast paths*: lists that are block structured over the balls (e.g. `h1, v1, h2, v2, ...`,
  all event indicators, all `e<i>`, or a single global variable) are processed type-major: one tight
  loop over the balls per variable type, seeds read in place through precomputed column pointers and
  strides; no scatter, no per-element value-reference decoding, no branches on list membership.
* *Struct-of-arrays generic path*: arbitrary lists are scattered into one contiguous scratch array per
  variable type, evaluated, and cleared again. Both paths share the formulas and are checked against
  each other.
* *Lazy event map*: event handling only records which branch each ball took (1 byte); the partial
  derivatives of the event map are evaluated when they are requested.

**What remains inherent**: FMI passes value reference lists with every derivative call, so reading
these lists once is a lower bound for a call. And every FMI call from Julia goes through FMIImport.jl
(status checks, mode tracking), which costs the same per call for both approaches.

Measured per-call costs (N = 1000, called from Julia through FMIImport.jl):

Standard FMI 3.0 functions (identical code in both builds; ranges over two runs):

| call | time |
|---|---:|
| `fmi3SetContinuousStates` (2000 states) | 1.2 - 1.6 µs |
| `fmi3GetContinuousStateDerivatives` | 1.2 - 1.5 µs |
| `fmi3GetEventIndicators` | 0.9 - 1.2 µs |
| `fmi3GetFMUState` (reused state) | 1.4 - 2.3 µs |
| `fmi3SetFMUState` | 1.4 - 2.2 µs |
| event iteration (`fmi3EnterEventMode`, `fmi3UpdateDiscreteStates`, `fmi3EnterContinuousTimeMode`) | 1.5 - 1.9 µs |

LS-SA functions:

| call | plan cached | no plan cache (cold) |
|---|---:|---:|
| JVP `dz/dx` (1000 x 2000) | 1.5 µs | 4.3 µs |
| VJP `dz/dx` (1000 x 2000) | 2.2 µs | 5.3 µs |
| JVP `dxdot/de` (2000 x 1000) | 4.2 µs | 7.2 µs |
| VJP `dxdot/de` (2000 x 1000) | 1.8 µs | 4.9 µs |
| single entry `dz1/dh1` | 0.25 µs | 0.26 µs |

Measured on an Intel Core Ultra 7 165U (hybrid CPU: the core type the OS schedules on causes most of the
variation), Julia 1.12.4. For comparison: a single finite-difference evaluation of the event indicators needs at
least `fmi3SetContinuousStates` + `fmi3GetEventIndicators` (about 2.5 µs); a complete LS-SA JVP costs less.

## Versions

* `0.1.0`: initial prototype (MinGW GCC 8.1, N = 1, 10, 100, 1000). Dense derivative evaluation,
  partials of `nextEventTime` w.r.t. states and of `x+` w.r.t. `pre(x)` missing, no
  `providesAdjointDerivatives` flag.
* `0.2.0`: everything above (N = 1, 3, 10, 30, 100, 300, 1000), binaries for x86_64-windows
  (MSVC 2022 `/O2`, static runtime) and x86_64-linux (GCC 15 `-O3`, needs glibc >= 2.14), declared
  `dependencies`; coupled variant `NCoupledBouncingBalls_<N>_K<K>` (N = 3 ... 1000, K = 1, 8, N - 1).
  Note the changed time event scheduling (first kick no longer at `0.7 * base_event_period`).

## Building

`build_fmus.ps1` (Windows, MSVC) and `build_fmus.sh` (Linux, GCC) build all ball counts with CMake and
add the binaries of their platform to `models/bin/LSSAReferenceFMUs/<version>/3.0/NBouncingBalls_<N>.fmu`;
binaries of the other platform in an existing FMU are kept. Both scripts stop if the existing FMU has a
different `modelDescription.xml` (its binaries were built from different sources). Variants go to a
separate build and output directory, e.g. without the plan cache:
`build_fmus.ps1 -BuildDir build-nocache -OutDir <dir> -CMakeArgs "-DLSSA_PLAN_CACHE=OFF"` or
`./build_fmus.sh -b build-linux-nocache -o <dir> -- -DLSSA_PLAN_CACHE=OFF`. The coupled variant:
`build_fmus.ps1 -Counts 3,10,30,100,300,1000 -Neighbors 1,8,-1` / `./build_fmus.sh -n "3 10 30 100 300 1000" -k "1 8 -1"`
(`-1`: dense).
