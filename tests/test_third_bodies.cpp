#include "third_bodies.hpp"
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>
namespace fs=std::filesystem;
static void expect(bool okay, const std::string& why) {
    if(!okay) throw std::runtime_error(why);
}
static double relative(const Vec3& a,const Vec3& b) {
    return (a-b).norm()/std::max({1e-200,a.norm(),b.norm()});
}
int main() {
 try {
    using namespace thirdbody;
    expect(select("none")==std::array<bool,3>({false,false,false}),"none disables all");
    expect(select("earth,SUN,jupiter")==std::array<bool,3>({true,true,true}),"body selectors");
    bool bad=false;
    try {select("earth,mars");} catch(const std::invalid_argument&){bad=true;}
    expect(bad,"reject unsupported body");
    Perturbations none("none","none",""), earth("earth","analytic",""),
                  sun("sun","analytic",""),jup("jupiter","analytic",""),
                  all("all","analytic","");
    Vec3 q{1738000,335000,-275000};
    for(double t:{0.,2800.,1234567.,31536000.}) {
        expect(none.acceleration(q,t).norm()==0,"disabled perturbations have no force");
        const Vec3 R=earth.positions(t)[EARTH];
        const Vec3 d=R-q;
        const Vec3 expected=(d*(1/std::pow(d.norm(),3.))-
                             R*(1/std::pow(R.norm(),3.)))*GM_EARTH;
        expect(relative(earth.acceleration(q,t),expected)<5e-11,"missing indirect term or wrong sign");
        expect(relative(all.acceleration(q,t),earth.acceleration(q,t)+
                   sun.acceleration(q,t)+jup.acceleration(q,t))<1e-11,
               "sum of Earth/Sun/Jupiter forces");
        // Acceleration must vanish at the Moon center in lunar-centered coordinates.
        expect(all.acceleration(Vec3{},t).norm()<1e-18,"nonzero force at lunar center");
        const double cos_i=earth.reference_normal(t).dot(earth.reference_normal(t));
        expect(std::fabs(cos_i-1)<1e-12,"reference orbital plane has nonunit normal");
        // For r/R << 1 the exact tide tends to the usual quadrupole form.
        const double invR3=1/std::pow(R.norm(),3.);
        const Vec3 n=R*(1/R.norm());
        const Vec3 quadrupole=(n*(3*n.dot(q))-q)*(GM_EARTH*invR3);
        expect(relative(earth.acceleration(q,t),quadrupole)<0.02,
               "incorrect leading quadrupole tidal acceleration");
    }
    const fs::path file=fs::temp_directory_path()/"mascons_3body_polynomial.csv";
    {std::ofstream f(file);
      f<<"# synthetic cubic ephemeris for interpolation test\n"
          "time_s,earth_x_m,earth_y_m,earth_z_m,sun_x_m,sun_y_m,sun_z_m,jupiter_x_m,jupiter_y_m,jupiter_z_m\n";
      for(int j=-2;j<=7;++j) {
          double t=j*100.;
          // Each component is at most cubic, should interpolate to roundoff.
          f<<t;
          for(int b=0;b<3;++b)for(int k=0;k<3;++k)
              f<<","<<std::fixed<<std::setprecision(8)
               <<(1.0e10*(b+1)+1.e9*k+(b+1)*k*1.e5+31*t+0.02*t*t+1e-5*t*t*t);
          f<<'\n';
      }
    }
    Perturbations tab("earth,sun,jupiter","file",file.string());
    for(double t:{-195.,-10.,100.,250.,633.,699.}) {
        auto pos=tab.positions(t);
        for(int b=0;b<3;++b)for(int k=0;k<3;++k) {
            const double val=1e10*(b+1)+1e9*k+(b+1)*k*1e5+31*t+0.02*t*t+1e-5*t*t*t;
            const double got=(k==0?pos[b].x:(k==1?pos[b].y:pos[b].z));
            expect(std::fabs(got-val)<0.001,"cubic ephemeris interpolation error");
        }
    }
    bad=false;
    try{tab.positions(701.);}catch(const std::runtime_error&){bad=true;}
    expect(bad,"out-of-coverage ephemeris must fail");
    fs::remove(file);
    std::cout<<"PASS third-body differential gravity, quadrupole limit, selectors, "
               "reference plane, cubic table and bounds\n";
    return 0;
 } catch(const std::exception& e) {
    std::cerr<<"FAIL third bodies: "<<e.what()<<"\n";
    return 1;
 }
}
