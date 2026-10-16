#pragma once
#include "ensemble.h"

namespace boris {

// B(x) = B0*(1 + alpha*x) z-hat, E = 0. Divergence-free.
// Grad-B drift: v_dy ~= m*v_perp^2*alpha / (2*q*B0); needs alpha*r_L << 1.
struct LinGradB {
    double B0{1.0};
    double alpha{0.1};

    Fields eval(const Vec3& x) const {
        return {{0, 0, 0}, {0, 0, B0 * (1.0 + alpha * x.x)}};
    }
};

// Per-particle field from model.eval(); midpoint keeps 2nd order
inline void push_ensemble_field(Ensemble& e, const LinGradB& model, double dt) {
    for (std::size_t i = 0; i < e.size(); ++i) {
        State s{e.x[i], e.v[i]};
        Vec3 x_mid = s.x + s.v * (0.5 * dt);
        Fields f = model.eval(x_mid);
        Params p{e.q[i], e.m[i], dt};
        boris_push(s, f, p);
        e.x[i] = s.x;
        e.v[i] = s.v;
    }
}

} // namespace boris
