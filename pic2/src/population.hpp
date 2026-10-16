#pragma once
// Ion population: particles + their macroparticle weights, in one object.
//
// Abstract rule: weights live INSIDE the population (the slide's
// n_ij = sum_p w_p * S). Pushes ignore w, deposits read it, and each
// function takes only what it needs — deposit(pop) is one argument.
#include "ensemble.h"
#include <vector>

namespace boris {

struct Population {
    Ensemble e;             // positions, velocities, charges, masses
    std::vector<double> w;  // w[p] = physical ions represented by macro-p p

    std::size_t size() const { return e.size(); }

    // Uniform weights for mean density n0 over the domain volume:
    // w = n0 * V / Np, so a uniform plasma deposits n == n0.
    void uniform_weights(double n0, double volume) {
        w.assign(e.size(), n0 * volume / (double)e.size());
    }

    double total_weight() const {
        double s = 0.0;
        for (double v : w) s += v;
        return s;
    }
};

}  // namespace boris
