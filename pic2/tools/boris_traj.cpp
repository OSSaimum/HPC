// One particle -> trajectory CSV t,x,y,z,vx,vy,vz.
// boris_traj out.csv x0 y0 z0 vx0 vy0 vz0 Ex Ey Ez Bx By Bz q m dt n
#include "boris.h"
#include <cstdio>
#include <fstream>

int main(int argc, char** argv) {
    if (argc != 18) {
        std::fprintf(stderr,
                     "usage: %s out.csv x0 y0 z0 vx0 vy0 vz0 "
                     "Ex Ey Ez Bx By Bz q m dt n\n",
                     argv[0]);
        return 2;
    }
    int a = 2;
    boris::State s{{std::atof(argv[a]), std::atof(argv[a + 1]), std::atof(argv[a + 2])},
                   {std::atof(argv[a + 3]), std::atof(argv[a + 4]), std::atof(argv[a + 5])}};
    boris::Fields f{{std::atof(argv[a + 6]), std::atof(argv[a + 7]), std::atof(argv[a + 8])},
                    {std::atof(argv[a + 9]), std::atof(argv[a + 10]), std::atof(argv[a + 11])}};
    boris::Params p{std::atof(argv[a + 12]), std::atof(argv[a + 13]), std::atof(argv[a + 14])};
    int n = std::atoi(argv[a + 15]);

    std::ofstream o(argv[1]);
    o << "t,x,y,z,vx,vy,vz\n";
    char buf[256];
    std::snprintf(buf, sizeof buf, "0,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                  s.x.x, s.x.y, s.x.z, s.v.x, s.v.y, s.v.z);
    o << buf;
    for (int i = 1; i <= n; ++i) {
        boris_push(s, f, p);
        std::snprintf(buf, sizeof buf, "%g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                      i * p.dt, s.x.x, s.x.y, s.x.z, s.v.x, s.v.y, s.v.z);
        o << buf;
    }
    return 0;
}
