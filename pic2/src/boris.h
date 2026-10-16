#pragma once
#include "vec3.h"

namespace boris {

struct State {
    Vec3 x{};
    Vec3 v{};
};

struct Fields {
    Vec3 E{};
    Vec3 B{};
};

struct Params {
    double q{1.0};
    double m{1.0};
    double dt{0.0};
};

void boris_push(State& s, const Fields& f, const Params& p);

} // namespace boris
