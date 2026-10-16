// pybind11 module for pic2 phase 1: moments + PPC vocabulary.
//
//   import pic2_ext as px
//   pop = px.Population.from_arrays(x, v, q, m, w)  # (Np,3),(Np,3),(Np,),(Np,),(Np,)
//   g = px.Geom(); g.nx = g.ny = g.nz = 8; g.origin = ...; g.d = ...
//   mo = px.Moments(); mo.setup(g); mo.deposit(pop)
//
// Abstract rule: geometry and populations travel as OBJECTS, so no
// function takes more than 2 args (except from_arrays, where the 5
// numpy arrays must enter together — the one necessary exception).
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "moments.hpp"
#include "ppc.h"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace boris;

namespace {

Population arrays_to_population(py::array_t<double> x0, py::array_t<double> v0,
                                py::array_t<double> q, py::array_t<double> m,
                                py::array_t<double> w) {
    if (x0.ndim() != 2 || x0.shape(1) != 3 || v0.ndim() != 2 ||
        v0.shape(1) != 3)
        throw std::runtime_error("x0/v0 must be (Np,3)");
    if (x0.shape(0) != v0.shape(0) || q.ndim() != 1 || m.ndim() != 1 ||
        w.ndim() != 1 || q.shape(0) != x0.shape(0) ||
        m.shape(0) != x0.shape(0) || w.shape(0) != x0.shape(0))
        throw std::runtime_error("inconsistent Np in x0/v0/q/m/w");
    py::ssize_t np = x0.shape(0);
    auto rx = x0.unchecked<2>(), rv = v0.unchecked<2>();
    auto rq = q.unchecked<1>(), rm = m.unchecked<1>(), rw = w.unchecked<1>();
    Population p;
    p.e.x.resize((std::size_t)np);
    p.e.v.resize((std::size_t)np);
    p.e.q.resize((std::size_t)np);
    p.e.m.resize((std::size_t)np);
    p.w.resize((std::size_t)np);
    for (py::ssize_t i = 0; i < np; ++i) {
        std::size_t k = (std::size_t)i;
        p.e.x[k] = {rx(i, 0), rx(i, 1), rx(i, 2)};
        p.e.v[k] = {rv(i, 0), rv(i, 1), rv(i, 2)};
        p.e.q[k] = rq(i);
        p.e.m[k] = rm(i);
        p.w[k] = rw(i);
    }
    return p;
}

py::array_t<double> vec_to_array(const std::vector<double>& v) {
    py::array_t<double> a(std::vector<py::ssize_t>{(py::ssize_t)v.size()});
    auto r = a.mutable_unchecked<1>();
    for (std::size_t i = 0; i < v.size(); ++i) r((py::ssize_t)i) = v[i];
    return a;
}

py::array_t<double> vec3_to_array(const std::vector<Vec3>& v) {
    py::array_t<double> a(
        std::vector<py::ssize_t>{(py::ssize_t)v.size(), 3});
    auto r = a.mutable_unchecked<2>();
    for (std::size_t i = 0; i < v.size(); ++i) {
        r((py::ssize_t)i, 0) = v[i].x;
        r((py::ssize_t)i, 1) = v[i].y;
        r((py::ssize_t)i, 2) = v[i].z;
    }
    return a;
}

}  // namespace

PYBIND11_MODULE(pic2_ext, mod) {
    mod.doc() = "pic2 phase 1: CIC moments + PPC vocabulary";

    py::class_<Vec3>(mod, "Vec3")
        .def(py::init<>())
        .def(py::init<double, double, double>())
        .def_readwrite("x", &Vec3::x)
        .def_readwrite("y", &Vec3::y)
        .def_readwrite("z", &Vec3::z);

    py::class_<GridGeom>(mod, "Geom")
        .def(py::init<>())
        .def_readwrite("origin", &GridGeom::origin)
        .def_readwrite("d", &GridGeom::d)
        .def_readwrite("nx", &GridGeom::nx)
        .def_readwrite("ny", &GridGeom::ny)
        .def_readwrite("nz", &GridGeom::nz)
        .def("ncells", &GridGeom::ncells)
        .def("volume", &GridGeom::volume)
        .def("cell_volume", &GridGeom::cell_volume);

    py::class_<Population>(mod, "Population")
        .def(py::init<>())
        .def_static("from_arrays", &arrays_to_population, py::arg("x"),
                    py::arg("v"), py::arg("q"), py::arg("m"), py::arg("w"))
        .def("uniform_weights", &Population::uniform_weights, py::arg("n0"),
             py::arg("volume"))
        .def("total_weight", &Population::total_weight)
        .def("size", &Population::size);

    py::class_<Moments>(mod, "Moments")
        .def(py::init<>())
        .def("setup", [](Moments& mo, const GridGeom& g) { mo.alloc(g); })
        .def("deposit", [](Moments& mo, const Population& p) { mo.deposit(p); })
        .def("total", &Moments::total)
        .def("n_array",
             [](const Moments& mo) { return vec_to_array(mo.n); })
        .def("u_array",
             [](const Moments& mo) { return vec3_to_array(mo.u); });

    // One-call deposit: deposit_moments(pop, geom) -> {n, u}. 2 args.
    mod.def(
        "deposit_moments",
        [](const Population& pop, const GridGeom& g) {
            Moments mo;
            mo.alloc(g);
            mo.deposit(pop);
            return py::dict("n"_a = vec_to_array(mo.n),
                            "u"_a = vec3_to_array(mo.u));
        },
        py::arg("pop"), py::arg("geom"));

    // PPC vocabulary (phase-1 bodies: deposit-only passes).
    py::class_<StepParams>(mod, "StepParams").def(py::init<>()).def_readwrite(
        "dt", &StepParams::dt);
    py::class_<HybridState>(mod, "HybridState")
        .def(py::init<>())
        .def_readwrite("ions", &HybridState::ions)
        .def_readwrite("moments", &HybridState::moments);
    py::class_<Stage>(mod, "Stage")
        .def_readonly("particles", &Stage::particles)
        .def_readonly("moments", &Stage::moments);
    mod.def("init_state", [](const Population& pop, const GridGeom& g,
                             const StepParams& p) { return init_state(pop, g, p); });
    mod.def("predictor1",
            [](const HybridState& s, const StepParams& p) { return predictor1(s, p); });
    mod.def("predictor2", [](const HybridState& s, const Stage& p1,
                             const StepParams& p) { return predictor2(s, p1, p); });
    mod.def("corrector", [](const HybridState& s, const Stage& p2,
                            const StepParams& p) { return corrector(s, p2, p); });
    mod.def("advance", [](const HybridState& s, const StepParams& p) {
        return advance(s, p);
    });
}
