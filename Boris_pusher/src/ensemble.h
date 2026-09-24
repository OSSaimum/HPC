#pragma once
// N-particle storage: struct-of-arrays. Particle i is (x[i], v[i], q[i], m[i]).
#include "boris.h"
#include <cstddef>
#include <vector>

namespace boris {

struct Ensemble {
    std::vector<Vec3> x;   // x[i] = position of particle i
    std::vector<Vec3> v;   // v[i] = velocity of particle i
    std::vector<double> q; // q[i] = charge of particle i
    std::vector<double> m; // m[i] = mass of particle i

    std::size_t size() const { return x.size(); }
};

// One step dt for all particles in the same uniform field via boris_push.
inline void push_ensemble(Ensemble& e, const Fields& f, double dt) {
    for (std::size_t i = 0; i < e.size(); ++i) {
        State s{e.x[i], e.v[i]};
        Params p{e.q[i], e.m[i], dt};
        boris_push(s, f, p);
        e.x[i] = s.x;
        e.v[i] = s.v;
    }
}

} // namespace boris
