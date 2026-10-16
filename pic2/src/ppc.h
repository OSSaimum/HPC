#pragma once
// Hybrid PIC loop roadmap: 3-step iterated Crank-Nicolson (Screenshot 2).
//
// One step reads EXACTLY like the slide:
//
//   auto p1 = predictor1(state, params);     // E, B, n, u predicted
//   auto p2 = predictor2(state, p1, params); // refined prediction
//   return corrector(state, p2, params);     // corrected — done
//
// Abstract rule: states in, states out, <=2 args each. Stages never
// chain (every push restarts FROM time n). E is diagnosed by Ohm's law,
// never time-advanced; Faraday (dB/dt = -curl E) is the only field
// update. Smoothing/binomial, YeeGrid, and Ohm arrive in later phases —
// their call sites are already fixed below so later work only FILLS IN
// bodies, never refactors signatures.
//
// Phase map:
//   phase 1 (now):  Population + Moments + GridGeom (deposit only)
//   phase 2:        + YeeGrid (gather fields->particles, Faraday, average)
//   phase 3:        + Ohm (E from moments + B)
//   phase 4:        + bodies of predictor1/predictor2/corrector/advance
#include "moments.hpp"
#include "population.hpp"

namespace boris {

// Everything a step needs, in one object (never dt/params scattered).
struct StepParams {
    double dt{0.0};
    // Phase 3+: OhmParams ohm;  Phase 2+: int smooth{1};
};

// Everything the plasma IS at one time level.
// Phase 2+: add YeeGrid fields (E, B staggered).
struct HybridState {
    Population ions;  // r, v, q, m, w
    Moments moments;  // n, u nodal (deposited from ions)
};

// One predicted time level: same content, distinct type so prediction
// vs corrected state stay distinct in signatures.
struct Stage {
    Population particles;
    Moments moments;
    // Phase 2+: YeeGrid fields;
};

// Consistent initialization: moments = deposit(ions).
// Phase 3+: also E = Ohm(moments, B).
//   (1 argstruct + 1 population = 2 args)
HybridState init_state(const Population& ions, const GridGeom& geom,
                       const StepParams& p);

// predictor1: state at n -> predicted stage.
// Body (phase 4): B* = B^n - dt*curl(E^n); push in <E,B>^{n+1/2};
// moments of pushed ions; E* = Ohm(moments, B*).
Stage predictor1(const HybridState& s, const StepParams& p);

// predictor2: refined prediction from (n, p1).
// Body (phase 4): <E,B> = average(n, p1); B** = B^n - dt*curl(<E>);
// push in <E,B>; moments; E** = Ohm(moments, B**).
Stage predictor2(const HybridState& s, const Stage& p1, const StepParams& p);

// corrector: final state from (n, p2).
// Body (phase 4): B^{n+1} = B^n - dt*curl(<E>) [trapezoidal -> 2nd order];
// push in <E,B>; moments; E^{n+1} = Ohm(final).
HybridState corrector(const HybridState& s, const Stage& p2,
                      const StepParams& p);

// One full hybrid step: predictor1 -> predictor2 -> corrector.
HybridState advance(const HybridState& s, const StepParams& p);

}  // namespace boris
