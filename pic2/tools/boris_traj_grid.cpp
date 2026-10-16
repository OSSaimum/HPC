// Ensemble with grid-sampled field -> trajectory CSV t,pid,x,y,z,vx,vy,vz.
// boris_traj_grid out.csv init.csv MODEL B0 alpha n ox oy oz L dt nsteps
// MODEL = lingrad | uniform; grid: n^3 nodes over periodic [o, o+L)^3.
#include "field_grid.h"
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>

int main(int argc, char** argv) {
    if (argc != 13) {
        std::fprintf(stderr,
                     "usage: %s out.csv init.csv MODEL B0 alpha n ox oy oz L dt nsteps\n",
                     argv[0]);
        return 2;
    }
    bool is_lingrad = std::strcmp(argv[3], "lingrad") == 0;
    bool is_uniform = std::strcmp(argv[3], "uniform") == 0;
    if (!is_lingrad && !is_uniform) {
        std::fprintf(stderr, "unknown field model '%s' (lingrad | uniform)\n", argv[3]);
        return 2;
    }

    // Read the ensemble (skip the header line).
    boris::Ensemble e;
    {
        std::ifstream in(argv[2]);
        if (!in) {
            std::fprintf(stderr, "cannot open %s\n", argv[2]);
            return 1;
        }
        std::string line;
        std::getline(in, line);
        double px, py, pz, vx, vy, vz, q, m;
        while (in) {
            std::getline(in, line);
            if (line.empty()) continue;
            for (char& c : line)
                if (c == ',') c = ' ';
            if (std::sscanf(line.c_str(), "%lf %lf %lf %lf %lf %lf %lf %lf",
                             &px, &py, &pz, &vx, &vy, &vz, &q, &m) == 8) {
                e.x.push_back({px, py, pz});
                e.v.push_back({vx, vy, vz});
                e.q.push_back(q);
                e.m.push_back(m);
            }
        }
    }
    if (e.size() == 0) {
        std::fprintf(stderr, "no particles found in %s\n", argv[2]);
        return 1;
    }

    // Build the grid from the named analytic model.
    double B0 = std::atof(argv[4]), alpha = std::atof(argv[5]);
    int n = std::atoi(argv[6]);
    double ox = std::atof(argv[7]), oy = std::atof(argv[8]), oz = std::atof(argv[9]);
    double L = std::atof(argv[10]);
    double dt = std::atof(argv[11]);
    int nsteps = std::atoi(argv[12]);

    boris::FieldGrid g;
    g.nx = g.ny = g.nz = n;
    g.origin = {ox, oy, oz};
    g.d = {L / n, L / n, L / n};
    if (is_lingrad)
        g.fill_lingrad(B0, alpha);
    else
        g.fill_uniform(B0);

    std::ofstream o(argv[1]);
    o << "t,pid,x,y,z,vx,vy,vz\n";
    char buf[256];
    auto dump_block = [&](double t) {
        for (std::size_t i = 0; i < e.size(); ++i) {
            std::snprintf(buf, sizeof buf, "%g,%zu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                          t, i, e.x[i].x, e.x[i].y, e.x[i].z,
                          e.v[i].x, e.v[i].y, e.v[i].z);
            o << buf;
        }
    };
    dump_block(0.0);
    for (int step = 1; step <= nsteps; ++step) {
        push_ensemble_grid(e, g, dt);
        dump_block(step * dt);
    }
    std::printf("advanced %zu particles x %d steps (grid %d^3) -> %s\n",
                e.size(), nsteps, n, argv[1]);
    return 0;
}
