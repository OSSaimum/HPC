#pragma once
// Ion moments: particle -> node deposition (CIC / B-spline order 1).
//
// Screenshot "Accumulate moments":
//   N_j = sum_p w_p * W(r_p - r_j) / dV
//   u_j = sum_p w_p v_p * W(r_p - r_j) / N_j
//
// W is the CIC tent (each particle touches 8 surrounding nodes).
// Abstract rule: Moments OWNS its GridGeom; deposit(pop) and
// alloc(geom) are METHODS (<=1 arg), never free functions with
// origin/d/nx/ny/nz lists. Serial loops by design.
#include "grid_geom.h"
#include "population.hpp"
#include <cmath>
#include <vector>

namespace boris {

struct Moments {
    GridGeom grid;          // nodes this instance lives on
    std::vector<double> n;  // n[node] number DENSITY (weight per volume)
    std::vector<Vec3> u;    // u[node] bulk velocity

    void alloc(const GridGeom& g) {
        grid = g;
        n.assign(g.ncells(), 0.0);
        u.assign(g.ncells(), {0, 0, 0});
    }

    // Deposit one population onto OUR nodes (zeroed first).
    // Near-vacuum nodes (n < 1e-10 of mean) get u = 0: the quotient
    // momentum/n is meaningless there. Uniform plasmas never hit it.
    void deposit(const Population& pop) {
        std::size_t nn = grid.ncells();
        n.assign(nn, 0.0);
        u.assign(nn, {0, 0, 0});
        const double dV = grid.cell_volume();
        for (std::size_t p = 0; p < pop.size(); ++p) {
            double gx = (pop.e.x[p].x - grid.origin.x) / grid.d.x;
            double gy = (pop.e.x[p].y - grid.origin.y) / grid.d.y;
            double gz = (pop.e.x[p].z - grid.origin.z) / grid.d.z;
            int i0 = (int)std::floor(gx), j0 = (int)std::floor(gy),
                k0 = (int)std::floor(gz);
            double fx = gx - i0, fy = gy - j0, fz = gz - k0;
            for (int a = 0; a < 2; ++a)
                for (int b = 0; b < 2; ++b)
                    for (int c = 0; c < 2; ++c) {
                        double s = (a ? fx : 1 - fx) * (b ? fy : 1 - fy) *
                                   (c ? fz : 1 - fz);
                        std::size_t id = grid.idx(i0 + a, j0 + b, k0 + c);
                        n[id] += pop.w[p] * s / dV;
                        u[id] += pop.e.v[p] * (pop.w[p] * s / dV);
                    }
        }
        for (std::size_t id = 0; id < nn; ++id)
            if (n[id] > 0.0) u[id] = u[id] / n[id];
        double floor = 1e-10 * pop.total_weight() / grid.volume();
        for (std::size_t id = 0; id < nn; ++id)
            if (n[id] < floor) {
                n[id] = 0.0;
                u[id] = {0, 0, 0};
            }
    }

    // Charge conservation: sum(n)*dV == sum(w) up to roundoff.
    double total() const {
        double s = 0.0;
        for (double v : n) s += v;
        return s * grid.cell_volume();
    }
};

// CIC tent: W(d) = max(0, 1 - |d|), d in cell units. Support [-1, 1].
inline double cic_shape(double d) {
    double a = std::abs(d);
    return a < 1.0 ? 1.0 - a : 0.0;
}

}  // namespace boris
