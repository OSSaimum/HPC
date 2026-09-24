// Ensemble in a shared uniform field -> trajectory CSV t,pid,x,y,z,vx,vy,vz.
// boris_traj_n out.csv init.csv Ex Ey Ez Bx By Bz dt nsteps
// init.csv rows: x,y,z,vx,vy,vz,q,m (header line first).
#include "ensemble.h"
#include <cstdio>
#include <fstream>
#include <string>

int main(int argc, char** argv) {
    if (argc != 11) {
        std::fprintf(stderr,
                     "usage: %s out.csv init.csv Ex Ey Ez Bx By Bz dt nsteps\n",
                     argv[0]);
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

    int a = 3;
    boris::Fields f{{std::atof(argv[a]), std::atof(argv[a + 1]), std::atof(argv[a + 2])},
                    {std::atof(argv[a + 3]), std::atof(argv[a + 4]), std::atof(argv[a + 5])}};
    double dt = std::atof(argv[a + 6]);
    int nsteps = std::atoi(argv[a + 7]);

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
        push_ensemble(e, f, dt);
        dump_block(step * dt);
    }
    std::printf("advanced %zu particles x %d steps -> %s\n", e.size(), nsteps, argv[1]);
    return 0;
}
