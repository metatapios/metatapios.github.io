#include "sph_gravity.hpp"
#include <cstdlib>
#include <iomanip>
#include <iostream>
int main(int argc,char** argv) {
    if(argc!=6) return 2;
    try {
        SphGravity g;g.load(argv[1],std::atoi(argv[2]));
        Vec3 q{std::atof(argv[3]),std::atof(argv[4]),std::atof(argv[5])};
        Vec3 a=g.acceleration(q);
        std::cout<<std::setprecision(17)<<g.potential(q)<<" "<<a.x<<" "<<a.y<<" "<<a.z<<"\n";
        return 0;
    } catch(const std::exception& e) { std::cerr<<e.what()<<"\n";return 3; }
}
