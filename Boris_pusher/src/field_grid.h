#pragma once
// E/B on a uniform 3D node mesh over a periodic domain; trilinear gather.
//   f = c0 + c1*fx + c2*fy + c3*fz + c4*fx*fy + c5*fy*fz + c6*fz*fx + c7*fx*fy*fz
// with (fx,fy,fz) in [0,1) and c0..c7 from the 8 corners (prototype form).
// One sample() reads 8 nodes x 6 components: memory-bound, SoA-friendly.
#include "ensemble.h"
#include <cmath>
#include <vector>

namespace boris {

// Trilinear interpolation of 8 corners, order subgrid[ix, iy, iz].
inline double tri8(double c000, double c100, double c010, double c001,
                   double c110, double c011, double c101, double c111,
                   double fx, double fy, double fz) {
    double c0 = c000;
    double c1 = c100 - c0;
    double c2 = c010 - c0;
    double c3 = c001 - c0;
    double c4 = c110 - c010 - c1;
    double c5 = c011 - c001 - c2;
    double c6 = c101 - c001 - c1;
    double c7 = c111 - c011 - c101 - c110 + c100 + c010 + c001 - c000;
    return c0 + c1 * fx + c2 * fy + c3 * fz + c4 * fx * fy +
           c5 * fy * fz + c6 * fz * fx + c7 * fx * fy * fz;
}

struct FieldGrid {
    int nx{0}, ny{0}, nz{0}; // node counts per axis
    Vec3 origin{};           // domain start, e.g. (0, 0, 0)
    Vec3 d{1, 1, 1};         // node spacing per axis
    std::vector<Vec3> E;     // E[node], collocated with B
    std::vector<Vec3> B;     // B[node], flat index (i*ny + j)*nz + k

    static int mod(int a, int n) {
        int r = a % n;
        return r < 0 ? r + n : r;
    }

    std::size_t idx(int i, int j, int k) const {
        return (std::size_t)((mod(i, nx) * ny + mod(j, ny)) * nz + mod(k, nz));
    }

    // Fold a position into [origin, origin + n*d).
    Vec3 wrap(const Vec3& x) const {
        auto fold = [](double v, double o, double step, int n) {
            double t = (v - o) / step;
            double f = t - std::floor(t / n) * n;
            return o + f * step;
        };
        return {fold(x.x, origin.x, d.x, nx),
                fold(x.y, origin.y, d.y, ny),
                fold(x.z, origin.z, d.z, nz)};
    }

    // Trilinear gather of all 6 components (periodic).
    Fields sample(const Vec3& x) const {
        Vec3 xw = wrap(x);
        double gx = (xw.x - origin.x) / d.x;
        double gy = (xw.y - origin.y) / d.y;
        double gz = (xw.z - origin.z) / d.z;
        int i0 = (int)std::floor(gx), j0 = (int)std::floor(gy), k0 = (int)std::floor(gz);
        double fx = gx - i0, fy = gy - j0, fz = gz - k0;

        const Vec3* e[2][2][2];
        const Vec3* b[2][2][2];
        for (int a = 0; a < 2; ++a)
            for (int c = 0; c < 2; ++c)
                for (int d_ = 0; d_ < 2; ++d_) {
                    std::size_t id = idx(i0 + a, j0 + c, k0 + d_);
                    e[a][c][d_] = &E[id];
                    b[a][c][d_] = &B[id];
                }
        Fields f;
        f.E.x = tri8(e[0][0][0]->x, e[1][0][0]->x, e[0][1][0]->x, e[0][0][1]->x,
                     e[1][1][0]->x, e[0][1][1]->x, e[1][0][1]->x, e[1][1][1]->x,
                     fx, fy, fz);
        f.E.y = tri8(e[0][0][0]->y, e[1][0][0]->y, e[0][1][0]->y, e[0][0][1]->y,
                     e[1][1][0]->y, e[0][1][1]->y, e[1][0][1]->y, e[1][1][1]->y,
                     fx, fy, fz);
        f.E.z = tri8(e[0][0][0]->z, e[1][0][0]->z, e[0][1][0]->z, e[0][0][1]->z,
                     e[1][1][0]->z, e[0][1][1]->z, e[1][0][1]->z, e[1][1][1]->z,
                     fx, fy, fz);
        f.B.x = tri8(b[0][0][0]->x, b[1][0][0]->x, b[0][1][0]->x, b[0][0][1]->x,
                     b[1][1][0]->x, b[0][1][1]->x, b[1][0][1]->x, b[1][1][1]->x,
                     fx, fy, fz);
        f.B.y = tri8(b[0][0][0]->y, b[1][0][0]->y, b[0][1][0]->y, b[0][0][1]->y,
                     b[1][1][0]->y, b[0][1][1]->y, b[1][0][1]->y, b[1][1][1]->y,
                     fx, fy, fz);
        f.B.z = tri8(b[0][0][0]->z, b[1][0][0]->z, b[0][1][0]->z, b[0][0][1]->z,
                     b[1][1][0]->z, b[0][1][1]->z, b[1][0][1]->z, b[1][1][1]->z,
                     fx, fy, fz);
        return f;
    }

    // Uniform B0 z-hat: periodic extension equals analytic field everywhere,
    // so a boundary crosser must match the folded analytic trajectory.
    void fill_uniform(double B0) {
        E.assign((std::size_t)nx * ny * nz, {0, 0, 0});
        B.assign((std::size_t)nx * ny * nz, {0, 0, B0});
    }

    // Sample B(x) = B0*(1 + alpha*x) z-hat, E = 0; node (i,j,k) at origin + (i*dx, j*dy, k*dz).
    void fill_lingrad(double B0, double alpha) {
        E.assign((std::size_t)nx * ny * nz, {0, 0, 0});
        B.assign((std::size_t)nx * ny * nz, {0, 0, 0});
        for (int i = 0; i < nx; ++i)
            for (int j = 0; j < ny; ++j)
                for (int k = 0; k < nz; ++k) {
                    double x = origin.x + i * d.x;
                    std::size_t id = (std::size_t)((i * ny + j) * nz + k);
                    B[id] = {0, 0, B0 * (1.0 + alpha * x)};
                }
    }
};

// One step dt per particle from the grid-sampled midpoint field; wrap positions.
inline void push_ensemble_grid(Ensemble& e, const FieldGrid& g, double dt) {
    for (std::size_t i = 0; i < e.size(); ++i) {
        State s{e.x[i], e.v[i]};
        Vec3 x_mid = s.x + s.v * (0.5 * dt);
        Fields f = g.sample(x_mid);
        Params p{e.q[i], e.m[i], dt};
        boris_push(s, f, p);
        e.x[i] = g.wrap(s.x);
        e.v[i] = s.v;
    }
}

} // namespace boris
