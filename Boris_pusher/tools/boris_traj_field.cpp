// Ensemble with per-particle analytic field -> trajectory CSV t,pid,x,y,z,vx,vy,vz.
// boris_traj_field out.csv init.csv lingrad B0 alpha dt nsteps
#include "fields_model.h"
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>

int main(int argc, char** argv) {
    if (argc != 8) {
        std::fprintf(stderr,
                     "usage: %s out.csv init.csv lingrad B0 alpha dt nsteps\n",
                     argv[0]);
        return 2;
    }
    if (std::strcmp(argv[3], "lingrad") != 0) {
        std::fprintf(stderr, "unknown field model '%s' (only 'lingrad' yet)\n", argv[3]);
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

    boris::LinGradB model{std::atof(argv[4]), std::atof(argv[5])};
    double dt = std::atof(argv[6]);
    int nsteps = std::atoi(argv[7]);

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
        push_ensemble_field(e, model, dt);
        dump_block(step * dt);
    }
    std::printf("advanced %zu particles x %d steps (lingrad) -> %s\n",
                e.size(), nsteps, argv[1]);
    return 0;
}
