// Phase-1 bodies of the PPC roadmap: deposit-only passes.
//
// These compile and run NOW (no fields/Ohm yet). Each stage re-deposits
// the time-n ions so the predictor1/predictor2/corrector vocabulary is
// already callable and testable. Phase 4 replaces each body with
// faraday + push in half-step fields + Ohm — signatures UNCHANGED.
#include "ppc.h"

namespace boris {

HybridState init_state(const Population& ions, const GridGeom& geom,
                       const StepParams& p) {
    (void)p;  // dt unused until pushes land (phase 4)
    HybridState s;
    s.ions = ions;
    s.moments.alloc(geom);
    s.moments.deposit(ions);
    return s;
}

Stage predictor1(const HybridState& s, const StepParams& p) {
    (void)p;
    Stage p1;
    p1.particles = s.ions;  // phase 4: push in <E,B>^{n+1/2} first
    p1.moments.alloc(s.moments.grid);
    p1.moments.deposit(p1.particles);
    return p1;
}

Stage predictor2(const HybridState& s, const Stage& p1, const StepParams& p) {
    (void)p1;  // phase 4: average(n, p1) fields refine the push
    (void)p;
    Stage p2;
    p2.particles = s.ions;  // phase 4: push in refined <E,B>
    p2.moments.alloc(s.moments.grid);
    p2.moments.deposit(p2.particles);
    return p2;
}

HybridState corrector(const HybridState& s, const Stage& p2,
                      const StepParams& p) {
    (void)p2;  // phase 4: trapezoidal Faraday + final push + Ohm
    (void)p;
    HybridState next;
    next.ions = s.ions;
    next.moments.alloc(s.moments.grid);
    next.moments.deposit(next.ions);
    return next;
}

HybridState advance(const HybridState& s, const StepParams& p) {
    Stage p1 = predictor1(s, p);
    Stage p2 = predictor2(s, p1, p);
    return corrector(s, p2, p);
}

}  // namespace boris
