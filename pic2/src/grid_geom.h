#pragma once
// Grid geometry: the 5 numbers every grid object needs, in one place.
//
// Abstract rule: geometry travels as ONE object, never as
// (origin, d, nx, ny, nz) argument lists. Moments OWNS one;
// deposits and gathers take it via the owner.
#include "vec3.h"
#include <cstddef>

namespace boris {

struct GridGeom {
    Vec3 origin{};    // domain start, e.g. (0, 0, 0)
    Vec3 d{1, 1, 1};  // cell size per axis
    int nx{0}, ny{0}, nz{0};  // cell counts per axis (= node counts, periodic)

    std::size_t ncells() const {
        return (std::size_t)nx * (std::size_t)ny * (std::size_t)nz;
    }

    double volume() const { return d.x * nx * d.y * ny * d.z * nz; }

    double cell_volume() const { return d.x * d.y * d.z; }

    static int mod(int a, int n) {
        int r = a % n;
        return r < 0 ? r + n : r;
    }

    std::size_t idx(int i, int j, int k) const {
        return (std::size_t)((mod(i, nx) * ny + mod(j, ny)) * nz + mod(k, nz));
    }
};

}  // namespace boris
