// Time-aware velocity Verlet / Yoshida-4 / error-controlled Verlet.
// Inertial Hamiltonian H(q,p,t) = |p|^2/2 - U(Rz(-omega*t)q).
// accel(q,t) returns the INERTIAL acceleration at absolute time t.
// Fixed-step Verlet and Yoshida compositions are canonical symplectic maps
// for time-dependent separable Hamiltonians when evaluated at stage times.
// Adaptive step doubling is NOT symplectic; step-size variation is not either.
#pragma once
#include "vec3.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

// The time argument is in/out, and a is the acceleration at (q,t).
template <class Accel>
inline void verlet_step(Vec3& q, Vec3& v, double& t, double h, const Accel& accel,
                        Vec3& a) {
    v += a * (0.5 * h);
    q += v * h;
    t += h;
    a = accel(q, t);
    v += a * (0.5 * h);
}

template <class Accel>
inline void yoshida4_step(Vec3& q, Vec3& v, double& t, double h,
                          const Accel& accel, Vec3& a) {
    constexpr double c = 1.2599210498948731648;
    constexpr double w1 = 1.0 / (2.0 - c);
    constexpr double w0 = -c / (2.0 - c);
    verlet_step(q, v, t, w1 * h, accel, a);
    verlet_step(q, v, t, w0 * h, accel, a);
    verlet_step(q, v, t, w1 * h, accel, a);
}

// Dormand-Prince 5(4) for the first-order state (q,v). Returns the
// fifth-order estimate and embedded fourth-order error. Unlike Yoshida, this
// integration is adaptive and is NOT symplectic.
#include <array>
template <class Accel>
inline double rk45_trial(Vec3& q, Vec3& v, double t, double h, const Accel& accel) {
    constexpr double nodes[7]={0,1./5,3./10,4./5,8./9,1,1};
    constexpr double coeff[7][7]={
        {0},
        {1./5},
        {3./40,9./40},
        {44./45,-56./15,32./9},
        {19372./6561,-25360./2187,64448./6561,-212./729},
        {9017./3168,-355./33,46732./5247,49./176,-5103./18656},
        {35./384,0,500./1113,125./192,-2187./6784,11./84}
    };
    constexpr double w5[7]={35./384,0,500./1113,125./192,-2187./6784,11./84,0};
    constexpr double w4[7]={5179./57600,0,7571./16695,393./640,-92097./339200,187./2100,1./40};
    std::array<Vec3,7> kq{},kv{};
    for (int i=0;i<7;++i) {
        Vec3 qi=q,vi=v;
        for (int j=0;j<i;++j) {
            qi+=kq[j]*(h*coeff[i][j]);
            vi+=kv[j]*(h*coeff[i][j]);
        }
        kq[i]=vi;
        kv[i]=accel(qi,t+h*nodes[i]);
    }
    Vec3 q5=q,v5=v,eq{},ev{};
    for (int i=0;i<7;++i) {
        q5+=kq[i]*(h*w5[i]);
        v5+=kv[i]*(h*w5[i]);
        eq+=kq[i]*(h*(w5[i]-w4[i]));
        ev+=kv[i]*(h*(w5[i]-w4[i]));
    }
    q=q5;v=v5;
    return std::max(eq.norm()/std::max(1.0,q.norm()),ev.norm()/std::max(1.0,v.norm()));
}

// Adaptive DOPRI5(4). rtol applies independently to position and velocity,
// using state magnitudes as scales. Throws explicitly if the tolerance cannot
// be achieved at the minimum step size; never loops indefinitely.
template <class Accel>
inline double adaptive_step(Vec3& q,Vec3& v,double& t,double& h_next,
                            const Accel& accel,Vec3& a,double rtol,
                            double h_min,double h_max) {
    if(!(rtol>0) || !(h_min>0) || !(h_max>=h_min))
        throw std::invalid_argument("invalid adaptive tolerances/limits");
    double h=std::clamp(h_next,h_min,h_max);
    for(int rejected=0;rejected<64;++rejected) {
        Vec3 qp=q,vp=v;
        double err=rk45_trial(qp,vp,t,h,accel);
        if(!std::isfinite(err)) throw std::runtime_error("non-finite adaptive error");
        double factor=std::clamp(0.9*std::pow(rtol/std::max(err,1e-300),0.2),0.2,4.0);
        if(err<=rtol) {
            q=qp;v=vp;t+=h;a=accel(q,t);
            h_next=std::clamp(h*factor,h_min,h_max);
            return h;
        }
        if(h<=h_min*(1.0+1e-12))
            throw std::runtime_error("adaptive error tolerance impossible at minimum step size");
        h=std::clamp(h*std::min(0.9,factor),h_min,h_max);
    }
    throw std::runtime_error("adaptive solver exceeded 64 rejected steps");
}

// DOPRI5(4) with per-component dimensional error control and continuous
// (quartic) dense output. Uses the stage at the accepted endpoint as the new
// acceleration (FSAL) so no redundant gravitational evaluation is required.
// Units: position [m], velocity [m/s], time [s]. Not symplectic.
struct DopriDense {
    Vec3 q0{}, v0{};
    double t0=0, h=0;
    std::array<Vec3,7> kq{},kv{};
    void evaluate(double theta,Vec3& q,Vec3& v) const {
        if(theta<0 || theta>1)throw std::invalid_argument("DOPRI dense theta outside [0,1]");
        // Shampine fourth-order continuous extension, used by scipy RK45.
        constexpr double P[7][4]={
            {1.,-8048581381./2820520608.,8663915743./2820520608.,-12715105075./11282082432.},
            {0.,0.,0.,0.},
            {0.,131558114200./32700410799.,-68118460800./10900136933.,87487479700./32700410799.},
            {0.,-1754552775./470086768.,14199869525./1410260304.,-10690763975./1880347072.},
            {0.,127303824393./49829197408.,-318862633887./49829197408.,701980252875./199316789632.},
            {0.,-282668133./205662961.,2019193451./616988883.,-1453857185./822651844.},
            {0.,40617522./29380423.,-110615467./29380423.,69997945./29380423.}
        };
        q=q0;v=v0;
        const double x2=theta*theta,x3=x2*theta,x4=x3*theta;
        for(int i=0;i<7;++i){
            const double b=P[i][0]*theta+P[i][1]*x2+P[i][2]*x3+P[i][3]*x4;
            q+=kq[i]*(h*b);v+=kv[i]*(h*b);
        }
    }
};
struct DopriTry {
    Vec3 q{},v{},a{};
    DopriDense dense{};
    double error=0;
};

template <class Accel>
inline DopriTry dopri_trial(const Vec3& q,const Vec3& v,const Vec3& a,
                            double t,double h,const Accel& accel,
                            double rtol,double atol_pos_m,double atol_vel_mps) {
    constexpr double nodes[7]={0,1./5,3./10,4./5,8./9,1,1};
    constexpr double A[7][7]={
        {0},{1./5},{3./40,9./40},{44./45,-56./15,32./9},
        {19372./6561,-25360./2187,64448./6561,-212./729},
        {9017./3168,-355./33,46732./5247,49./176,-5103./18656},
        {35./384,0,500./1113,125./192,-2187./6784,11./84}
    };
    constexpr double W5[7]={35./384,0,500./1113,125./192,-2187./6784,11./84,0};
    constexpr double W4[7]={5179./57600,0,7571./16695,393./640,-92097./339200,187./2100,1./40};
    DopriTry out;out.dense.q0=q;out.dense.v0=v;out.dense.t0=t;out.dense.h=h;
    auto& kq=out.dense.kq;auto& kv=out.dense.kv;
    kq[0]=v;kv[0]=a;
    for(int i=1;i<7;++i){
        Vec3 qi=q,vi=v;
        for(int j=0;j<i;++j){qi+=kq[j]*(h*A[i][j]);vi+=kv[j]*(h*A[i][j]);}
        kq[i]=vi;kv[i]=accel(qi,t+h*nodes[i]);
    }
    Vec3 eq{},ev{};out.q=q;out.v=v;
    for(int i=0;i<7;++i){
        out.q+=kq[i]*(h*W5[i]);out.v+=kv[i]*(h*W5[i]);
        eq+=kq[i]*(h*(W5[i]-W4[i]));ev+=kv[i]*(h*(W5[i]-W4[i]));
    }
    out.a=kv[6];
    auto component_error=[&](double x0,double x1,double e,double atol){
        return std::fabs(e)/(atol+rtol*std::max(std::fabs(x0),std::fabs(x1)));
    };
    out.error=std::max({
        component_error(q.x,out.q.x,eq.x,atol_pos_m),
        component_error(q.y,out.q.y,eq.y,atol_pos_m),
        component_error(q.z,out.q.z,eq.z,atol_pos_m),
        component_error(v.x,out.v.x,ev.x,atol_vel_mps),
        component_error(v.y,out.v.y,ev.y,atol_vel_mps),
        component_error(v.z,out.v.z,ev.z,atol_vel_mps)
    });
    return out;
}

template <class Accel>
inline double dopri_adaptive_step(Vec3& q,Vec3& v,double& t,double& dt_next,
                                  const Accel& accel,Vec3& a,
                                  double rtol,double atol_pos_m,double atol_vel_mps,
                                  double min_dt,double max_dt,DopriDense& accepted_dense,
                                  double* accepted_error=nullptr) {
    if(!(rtol>0 && atol_pos_m>0 && atol_vel_mps>0 && min_dt>0 && max_dt>=min_dt))
        throw std::invalid_argument("invalid DOPRI step controls");
    double h=std::clamp(dt_next,min_dt,max_dt);
    for(int iter=0;iter<64;++iter){
        DopriTry trial=dopri_trial(q,v,a,t,h,accel,rtol,atol_pos_m,atol_vel_mps);
        if(!std::isfinite(trial.error))throw std::runtime_error("non-finite DOPRI error");
        double factor=std::clamp(0.90*std::pow(1.0/std::max(trial.error,1e-300),0.2),0.2,4.0);
        if(trial.error<=1.0){
            q=trial.q;v=trial.v;a=trial.a;t+=h;
            accepted_dense=std::move(trial.dense);
            if(accepted_error)*accepted_error=trial.error;
            dt_next=std::clamp(h*factor,min_dt,max_dt);
            return h;
        }
        if(h<=min_dt*(1+1e-12))
            throw std::runtime_error("DOPRI error tolerance unattainable at minimum timestep");
        h=std::clamp(h*std::min(0.9,factor),min_dt,max_dt);
    }
    throw std::runtime_error("DOPRI exceeded 64 rejected steps");
}
