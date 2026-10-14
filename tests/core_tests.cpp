#include "vec3.hpp"
#include "sph_gravity.hpp"
#include "topography.hpp"
#include "integrator.hpp"
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace fs=std::filesystem;
static void check(bool c,const char* m) { if(!c) throw std::runtime_error(m); }
static bool close(double a,double b,double eps) {
    return std::fabs(a-b) <= eps * std::max({1.0,std::fabs(a),std::fabs(b)});
}
static void test_harmonics() {
    auto path=fs::temp_directory_path()/"mascons_unit_grav.bin";
    constexpr double GM=4.9028e12,R=1738000.0;
    std::array<double,6> C={1.,1e-5,2e-5,-1e-4,1e-5,2e-5};
    std::array<double,6> S={0.,0.,-2e-5,0.,0.,-1e-5};
    {std::ofstream f(path,std::ios::binary);
      int32_t degree=2;f.write((char*)&degree,sizeof(degree));
      f.write((char*)&GM,sizeof(GM));f.write((char*)&R,sizeof(R));
      for(int i=0;i<6;++i){f.write((char*)&C[i],8);f.write((char*)&S[i],8);}}
    SphGravity g;g.load(path.string(),2);
    fs::remove(path);
    // Explicit independent unassociated formulas (4pi norm, csphase=+1)
    double lat=0.37,lon=0.91,r=R+50000,s=std::sin(lat),c=std::cos(lat),z=R/r;
    double P10=std::sqrt(3.)*s,P11=std::sqrt(3.)*c;
    double P20=0.5*std::sqrt(5.)*(3*s*s-1);
    double P21=std::sqrt(15.)*s*c,P22=0.5*std::sqrt(15.)*c*c;
    double sum=1+z*(P10*C[1]+P11*(C[2]*std::cos(lon)+S[2]*std::sin(lon)))+
        z*z*(P20*C[3]+P21*(C[4]*std::cos(lon)+S[4]*std::sin(lon))+
        P22*(C[5]*std::cos(2*lon)+S[5]*std::sin(2*lon)));
    Vec3 x{r*c*std::cos(lon),r*c*std::sin(lon),r*s};
    check(close(g.potential(x),GM/r*sum,1e-14),"unnormalized tesseral coefficients");
    for(double ph:{lat,1.2,M_PI/2, -M_PI/2, M_PI/2-1e-9}) {
        Vec3 p{r*std::cos(ph)*std::cos(lon),r*std::cos(ph)*std::sin(lon),r*std::sin(ph)};
        double rel=(g.acceleration(p)-g.acceleration_fd(p)).norm()/g.acceleration_fd(p).norm();
        check(rel<2e-7,"gravity gradient differs from finite-difference near poles");
    }
}
static void test_topography_normal() {
    auto path=fs::temp_directory_path()/"mascons_unit_topo.bin";
    {std::ofstream f(path,std::ios::binary);
      int32_t ny=16,nx=32;f.write((char*)&ny,4);f.write((char*)&nx,4);
      for(int j=0;j<ny;++j)for(int i=0;i<nx;++i){
          double r=1737400.+100.*j+200.*i;f.write((char*)&r,8);
      }}
    Topography t;t.load(path.string());fs::remove(path);
    double lat=0.2,lon=1.1,r=t.surface_radius(lat,lon);
    Vec3 normal=t.surface_normal(lat,lon);
    check(close(normal.norm(),1,1e-14),"normal not unit");
    // Normal must annihilate infinitesimal surface tangents.
    double cl=std::cos(lat),sl=std::sin(lat),co=std::cos(lon),so=std::sin(lon);
    Vec3 er{cl*co,cl*so,sl}, elat{-sl*co,-sl*so,cl}, elon{-so,co,0};
    auto grad=t.surface_sample(lat,lon);
    Vec3 tl=er*grad.dlat+elat*r;
    Vec3 tn=er*grad.dlon+elon*(r*cl);
    check(std::fabs(normal.dot(tl))/tl.norm()<1e-13,"incorrect latitude surface normal");
    check(std::fabs(normal.dot(tn))/tn.norm()<1e-13,"incorrect longitude surface normal");
}
static std::array<double,2> forced(double dt,bool yoshida) {
    // a(t)=sin(t) -> v(t)=1-cos(t), x(t)=t-sin(t)
    auto a=[](const Vec3&,double t){return Vec3{std::sin(t),0,0};};
    Vec3 q,v,aa=a(q,0);double t=0;
    int N=static_cast<int>(std::round(6.4/dt));
    for(int i=0;i<N;++i) {
        if(yoshida)yoshida4_step(q,v,t,dt,a,aa);
        else verlet_step(q,v,t,dt,a,aa);
    }
    return {std::fabs(q.x-(t-std::sin(t))),std::fabs(v.x-(1-std::cos(t)))};
}
static void test_integrators(){
    for(bool y:{false,true}){
        auto coarse=forced(0.2,y),fine=forced(0.1,y);
        double desired=y?11.0:3.0;
        check(coarse[0]/fine[0]>desired,"integrator does not converge at expected order (x)");
        check(coarse[1]/fine[1]>desired,"integrator does not converge at expected order (v)");
    }
    auto a=[](const Vec3& x,double){return Vec3{-x.x,0,0};};
    Vec3 q{1,0,0},v{0,0,0},aa=a(q,0);double t=0,h=0.1;
    bool caught=false;
    try {adaptive_step(q,v,t,h,a,aa,1e-20,0.1,0.1);}
    catch(const std::runtime_error&){caught=true;}
    check(caught,"adaptive solver must reject unattainable tolerance at dt_min");
}
static void test_dopri_dense_and_adaptivity(){
    auto force=[](const Vec3& q,double){return Vec3{-q.x,-q.y,-q.z};};
    const Vec3 q0{1,0.25,0},v0{0,1,0};
    const Vec3 a0=force(q0,0);
    auto trial=dopri_trial(q0,v0,a0,0,0.2,force,1e-10,1e-4,1e-6);
    Vec3 qend{},vend{};
    trial.dense.evaluate(1.,qend,vend);
    check((trial.q-qend).norm()<1e-13,"DOPRI dense endpoint position mismatch");
    check((trial.v-vend).norm()<1e-13,"DOPRI dense endpoint velocity mismatch");
    Vec3 qstart{},vstart{};trial.dense.evaluate(0.,qstart,vstart);
    check((q0-qstart).norm()<1e-14 && (v0-vstart).norm()<1e-14,
          "DOPRI dense initial state mismatch");
    // Verify the adaptive controller actually adjusts accepted steps, and
    // accuracy tightens when rtol/atol are decreased.
    auto propagate=[&](double rtol,double pos_abs,double vel_abs){
        Vec3 q=q0,v=v0,a=a0;double t=0,h=0.8;
        int nstep=0;
        while(t<4.0-1e-12){
            DopriDense dense;
            const double cap=std::min(1.0,4.0-t);
            h=std::min(h,cap);
            dopri_adaptive_step(q,v,t,h,force,a,rtol,pos_abs,vel_abs,
                                std::min(1e-5,cap),cap,dense);
            ++nstep;
        }
        Vec3 expected=q0*std::cos(4.)+v0*std::sin(4.);
        return std::pair<double,int>((q-expected).norm(),nstep);
    };
    const auto loose=propagate(1e-5,1e-3,1e-5);
    const auto tight=propagate(1e-11,1e-7,1e-9);
    check(tight.first<loose.first && tight.second>loose.second,
          "adaptive tolerance does not improve error/step count");
    check(tight.first<1e-8,"adaptive long trajectory fails accuracy test");
}
int main(){
    try {
        test_harmonics();std::cout<<"PASS independent 4pi harmonic normalization and polar gradient\n";
        test_topography_normal();std::cout<<"PASS interpolated terrain surface normal\n";
        test_integrators();std::cout<<"PASS second/fourth-order time-dependent integration and adaptive failure\n";
        test_dopri_dense_and_adaptivity();std::cout<<"PASS DOPRI dimensional error controller and continuous extension\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr<<"FAIL: "<<e.what()<<"\n";return 1;
    }
}
