#include "boris.h"

namespace boris {

// Non-relativistic Boris: half-E kick, B rotation, half-E kick, centered drift.
void boris_push(State& s, const Fields& f, const Params& p) {
    const double q_m = p.q / p.m;
    const double h = 0.5 * p.dt;

    Vec3 v_minus = s.v + f.E * (q_m * h);

    Vec3 t = f.B * (q_m * h);
    Vec3 v_rot = t * (2.0 / (1.0 + t.norm2()));

    Vec3 v_prime = v_minus + cross(v_minus, t);
    Vec3 v_plus = v_minus + cross(v_prime, v_rot);
    Vec3 v_new = v_plus + f.E * (q_m * h);

    s.x += (s.v + v_new) * (0.5 * p.dt);
    s.v = v_new;
}

} // namespace boris
