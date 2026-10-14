// Low-lunar-orbit Monte Carlo in an inertial frame with a rotating lunar gravity
// field. The inertial acceleration is grad U(Rz(-omega*t)q) rotated back to
// inertial axes. The conserved integral is J = (v^2/2-U) - omega*(q x v)_z.
#include "vec3.hpp"
#include "sph_gravity.hpp"
#include "topography.hpp"
#include "integrator.hpp"
#include "third_bodies.hpp"

#include <algorithm>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdarg>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <optional>
#include <cstdint>
#include <functional>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;
static void tee(const char* fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    std::fflush(stdout);
}

struct Config {
    std::string grav = "data/gravity.bin", topo = "data/topo.bin";
    std::string out = "results/summary.csv", traj = "results/trajectories.csv";
    std::string log = "results/run.log";
    int n = 500, lmax = -1, steps_per_orbit = 480, traj_count = 0;
    int event_substeps = 4;
    // Adaptive clearance bracketing guards against short, grazing contacts.
    std::string event_search = "adaptive";
    int event_max_depth = 32;
    double alt_min = 25, alt_max = 500, days = 15;
    double omega = 2.661699e-6;
    double dense_envelope_km = 11, dense_factor = 1;
    double tol = 1e-10, traj_dt_min = 30, impact_tol_s = 1e-3;
    double atol_pos_m = 0.01, atol_vel_mps = 1e-5;
    double adaptive_max_dt_s = 90, adaptive_min_dt_s = 0.02;
    double adaptive_low_steps = 480, adaptive_high_steps = 120;
    double adaptive_scale_height_km = 75;
    std::string integrator = "adaptive";
    // Explicit master switch (off by default) for external tides/KL dynamics.
    // Individual bodies and ephemerides never activate external forces alone.
    bool kozai_lidov = false;
    bool third_bodies_specified = false;
    std::string third_bodies = "none", ephem_mode = "none", ephem_file;
    unsigned seed = 12345;
    std::string checkpoint_dir;
    bool resume = false;
    std::string hdf5;
    bool no_csv = false;
    // Explicit ensemble termination policy. A particle always terminates on impact;
    // in this independent-particle solver, all-impact termination is automatic.
    bool stop_when_all_impacted = false;
    std::string progress_dir;
    double progress_every_seconds = 600.0; // poll every 10 wall-clock minutes
    double progress_heartbeat_seconds = 1800.0; // unchanged count: every 30 minutes
    double progress_sim_hours = 6.0; // per-orbit heartbeat cadence
    double checkpoint_every_hours = 6;
};

struct Result {
    double lat0 = 0, lon0 = 0, alt0 = 0, heading = 0, inc0 = 0;
    int impacted = 0;
    double t_end_hr = 0, impact_lat = 0, impact_lon = 0;
    double impact_angle_deg = 0, radial_angle_deg = 0;
    double terminal_speed = 0, terrain_elev_m = 0;
    double vx_term = 0, vy_term = 0, vz_term = 0;
    double alt_min_seen = 0, alt_max_seen = 0, ecc_max = 0;
    int orbit_count = 0;
    long n_steps = 0;
    double inv0 = 0, invariant_drift_rel = 0;
    double lunar_jacobi_variation_rel = 0;
    double kl_incl0_deg = NAN, kl_incl_end_deg = NAN;
    double kl_jz0 = NAN, kl_jz_end = NAN;
};
struct TrajPoint {
    double t_hr, lat, lon, alt_km, invariant, ecc, x, y, z, vx, vy, vz;
    double kl_incl_deg, kl_jz;
};
struct State { Vec3 q, v, a; double t; };
struct Checkpoint {
    State state{};
    Result result{};
    double next_dt = 0, next_traj = 0, prev_z = 0, alt_start = 0;
    double next_checkpoint = 0;
    bool finished = false;
    std::vector<TrajPoint> history;
};

#include "progress.hpp"

// Native HDF5 restart snapshots retain full precision and a deterministic
// checksum; each particle has its current and previous atomic snapshot.
// HDF5 calls are mutex-serialized because the library may not be thread-safe.
static std::uint64_t fnv64(const std::string& data) {
    std::uint64_t h=14695981039346656037ULL;
    for(unsigned char b:data){h^=b;h*=1099511628211ULL;}
    return h;
}
static std::uint64_t file_hash(const std::string& path) {
    std::ifstream f(path,std::ios::binary);
    if(!f)throw std::runtime_error("cannot hash input file: "+path);
    std::uint64_t h=14695981039346656037ULL;
    char buffer[65536];
    while(f){f.read(buffer,sizeof buffer);
        for(std::streamsize i=0;i<f.gcount();++i){h^=static_cast<unsigned char>(buffer[i]);h*=1099511628211ULL;}
    }
    return h;
}
static std::string hex64(std::uint64_t h) {
    std::ostringstream x;x<<std::hex<<std::setw(16)<<std::setfill('0')<<h;return x.str();
}
static std::string run_fingerprint(const Config& c) {
    // The requested end date and CSV/log output paths may change on resume.
    // Everything affecting dynamics or trajectory/diagnostic sampling may not.
    std::ostringstream ss;ss<<std::scientific<<std::setprecision(17);
    ss<<"mascons_adaptive_checkpoint_v2 "<<c.n<<' '<<c.lmax<<' '<<c.steps_per_orbit<<' '
      <<c.traj_count<<' '<<c.event_substeps<<' '
      <<c.event_search<<' '<<c.event_max_depth<<' '<<c.alt_min<<' '<<c.alt_max<<' '
      <<c.omega<<' '<<c.dense_envelope_km<<' '<<c.dense_factor<<' '
      <<c.tol<<' '<<c.atol_pos_m<<' '<<c.atol_vel_mps<<' '
      <<c.adaptive_max_dt_s<<' '<<c.adaptive_min_dt_s<<' '
      <<c.adaptive_low_steps<<' '<<c.adaptive_high_steps<<' '
      <<c.adaptive_scale_height_km<<' '<<c.traj_dt_min<<' '<<c.impact_tol_s<<' '
      <<c.integrator<<' '<<c.seed<<' '<<c.kozai_lidov<<' '
      <<c.third_bodies<<' '<<c.ephem_mode<<' '
      <<hex64(file_hash(c.grav))<<' '<<hex64(file_hash(c.topo));
    if(!c.ephem_file.empty())ss<<' '<<hex64(file_hash(c.ephem_file));
    return hex64(fnv64(ss.str()));
}
static std::string checkpoint_path(const Config& c,int id){
    char name[64];std::snprintf(name,sizeof name,"particle_%08d.h5",id);
    return (fs::path(c.checkpoint_dir)/name).string();
}
static double read_real(std::istream& is) {
    std::string word;if(!(is>>word))throw std::runtime_error("short checkpoint float");
    std::size_t consumed=0;
    double v=std::stod(word,&consumed);
    if(consumed!=word.size())throw std::runtime_error("invalid checkpoint float: "+word);
    return v;
}
static void emit_real(std::ostream& os,double x){ os<<' '<<x; }
static void emit_result(std::ostream& os,const Result& r) {
    os<<"result";
    emit_real(os,r.lat0);emit_real(os,r.lon0);emit_real(os,r.alt0);
    emit_real(os,r.heading);emit_real(os,r.inc0);os<<' '<<r.impacted;
    emit_real(os,r.t_end_hr);emit_real(os,r.impact_lat);emit_real(os,r.impact_lon);
    emit_real(os,r.impact_angle_deg);emit_real(os,r.radial_angle_deg);
    emit_real(os,r.terminal_speed);emit_real(os,r.terrain_elev_m);
    emit_real(os,r.vx_term);emit_real(os,r.vy_term);emit_real(os,r.vz_term);
    emit_real(os,r.alt_min_seen);emit_real(os,r.alt_max_seen);emit_real(os,r.ecc_max);
    os<<' '<<r.orbit_count<<' '<<r.n_steps;
    emit_real(os,r.inv0);emit_real(os,r.invariant_drift_rel);
    emit_real(os,r.lunar_jacobi_variation_rel);emit_real(os,r.kl_incl0_deg);
    emit_real(os,r.kl_incl_end_deg);emit_real(os,r.kl_jz0);emit_real(os,r.kl_jz_end);
    os<<'\n';
}
static void read_result(std::istream& is,Result& r){
    r.lat0=read_real(is);r.lon0=read_real(is);r.alt0=read_real(is);
    r.heading=read_real(is);r.inc0=read_real(is);is>>r.impacted;
    r.t_end_hr=read_real(is);r.impact_lat=read_real(is);r.impact_lon=read_real(is);
    r.impact_angle_deg=read_real(is);r.radial_angle_deg=read_real(is);
    r.terminal_speed=read_real(is);r.terrain_elev_m=read_real(is);
    r.vx_term=read_real(is);r.vy_term=read_real(is);r.vz_term=read_real(is);
    r.alt_min_seen=read_real(is);r.alt_max_seen=read_real(is);r.ecc_max=read_real(is);
    is>>r.orbit_count>>r.n_steps;
    r.inv0=read_real(is);r.invariant_drift_rel=read_real(is);
    r.lunar_jacobi_variation_rel=read_real(is);r.kl_incl0_deg=read_real(is);
    r.kl_incl_end_deg=read_real(is);r.kl_jz0=read_real(is);r.kl_jz_end=read_real(is);
}
static void must_tag(std::istream& is,const char* expected){
    std::string tag;if(!(is>>tag)||tag!=expected)
        throw std::runtime_error(std::string("checkpoint expected '")+expected+"'");
}
static std::string checkpoint_payload(int id,const std::string& signature,const Checkpoint& p) {
    std::ostringstream os;os<<std::scientific<<std::setprecision(17);
    os<<"LUNARMC_CP_V1 "<<signature<<' '<<id<<' '<<int(p.finished)<<'\n';
    const auto& s=p.state;
    os<<"state";
    for(double d:{s.t,s.q.x,s.q.y,s.q.z,s.v.x,s.v.y,s.v.z,s.a.x,s.a.y,s.a.z})emit_real(os,d);
    os<<'\n';
    os<<"control";
    for(double d:{p.next_dt,p.next_traj,p.prev_z,p.alt_start,p.next_checkpoint})emit_real(os,d);
    os<<'\n';emit_result(os,p.result);
    os<<"history "<<p.history.size()<<'\n';
    for(const auto& t:p.history){
        os<<"sample";
        for(double d:{t.t_hr,t.lat,t.lon,t.alt_km,t.invariant,t.ecc,t.x,t.y,t.z,t.vx,t.vy,t.vz,t.kl_incl_deg,t.kl_jz})emit_real(os,d);
        os<<'\n';
    }
    os<<"END\n";
    return os.str();
}
static void atomic_write(const std::string& path,const std::string& content){
    // The temporary file belongs to exactly one orbit worker. On POSIX,
    // rename within the same directory replaces the old snapshot atomically.
    auto tmp=path+".tmp";
    {std::ofstream f(tmp,std::ios::binary|std::ios::trunc);
      if(!f)throw std::runtime_error("cannot open checkpoint temporary file: "+tmp);
      f.write(content.data(),static_cast<std::streamsize>(content.size()));f.flush();
      if(!f)throw std::runtime_error("cannot write checkpoint: "+tmp);
    }
    try {fs::rename(tmp,path);} catch(...) {fs::remove(tmp);throw;}
}
static void prepare_checkpoint_dir(const Config& c,const std::string& signature){
    if(c.checkpoint_dir.empty())return;
    fs::create_directories(c.checkpoint_dir);
    const auto manifest=(fs::path(c.checkpoint_dir)/"run.manifest").string();
    if(c.resume){
        std::ifstream f(manifest);std::string magic,stored;int count=-1;
        if(!(f>>magic>>stored>>count)||magic!="LUNARMC_RUN_V2"||stored!=signature||count!=c.n)
            throw std::runtime_error("missing/mismatched checkpoint manifest: "+manifest);
    } else {
        if(fs::exists(manifest))
            throw std::runtime_error("checkpoint run already exists; use --resume or a new --checkpoint-dir");
        for(const auto& e:fs::directory_iterator(c.checkpoint_dir)){
            if(e.path().extension()==".chk" || e.path().extension()==".h5" || e.path().extension()==".tmp")
                throw std::runtime_error("checkpoint directory contains existing snapshots: "+c.checkpoint_dir);
        }
        atomic_write(manifest,"LUNARMC_RUN_V2 "+signature+" "+std::to_string(c.n)+"\n");
    }
}


static Vec3 rot_z(const Vec3& v, double angle) {
    double c = std::cos(angle), s = std::sin(angle);
    return {v.x*c-v.y*s, v.x*s+v.y*c, v.z};
}
static Vec3 sphere_pos(double lat, double lon, double r) {
    double c = std::cos(lat);
    return {r*c*std::cos(lon), r*c*std::sin(lon), r*std::sin(lat)};
}
static void body_ll(const Vec3& x, double omega, double t, double& lat, double& lon) {
    Vec3 b = rot_z(x, -omega*t);
    lat = std::asin(std::clamp(b.z/b.norm(), -1.0, 1.0));
    lon = std::atan2(b.y,b.x);
}
static double clearance(const State& s, const Topography& topo, double omega) {
    double lat, lon; body_ll(s.q, omega, s.t, lat, lon);
    return s.q.norm() - topo.surface_radius(lat, lon);
}

// Reintegrate a fraction of the last ACCEPTED interval from its starting state.
// Used only for event probes/refinement. The step is always forward in time,
// even though Yoshida's internal composition contains a negative substep.
template <class Accel>
static State fraction_state(const State& start, double interval,
                             const Config& cfg, const Accel& accel) {
    State s = start;
    if (interval == 0) return s;
    if (cfg.integrator == "yoshida4")
        yoshida4_step(s.q,s.v,s.t,interval,accel,s.a);
    else if (cfg.integrator == "adaptive") {
        rk45_trial(s.q,s.v,s.t,interval,accel);
        s.t+=interval;
        s.a=accel(s.q,s.t);
    } else
        verlet_step(s.q,s.v,s.t,interval,accel,s.a);
    return s;
}

static Result propagate(const SphGravity& G, const Topography& T,
                        const thirdbody::Perturbations& P, const Config& cfg,
                        double lat0, double lon0, double alt0_km, double heading,
                        std::vector<TrajPoint>* traj,
                        const Checkpoint* old=nullptr,
                        const std::function<void(const Checkpoint&)>& save={},
                        const std::function<void(const State&,const Result&,double,double,double,double,double,bool,bool)>& telemetry={}) {
    Result R{};
    R.lat0=lat0*180/M_PI; R.lon0=lon0*180/M_PI;
    R.alt0=alt0_km; R.heading=heading*180/M_PI;
    double r0 = T.surface_radius(lat0,lon0) + alt0_km*1000.;
    Vec3 q0 = sphere_pos(lat0,lon0,r0);
    Vec3 east{-std::sin(lon0),std::cos(lon0),0};
    Vec3 north{-std::sin(lat0)*std::cos(lon0),-std::sin(lat0)*std::sin(lon0),std::cos(lat0)};
    Vec3 v0 = (north*std::cos(heading)+east*std::sin(heading)) * std::sqrt(G.GM/r0);
    Vec3 h0=q0.cross(v0);
    R.inc0=std::acos(std::clamp(h0.z/h0.norm(),-1.0,1.0))*180/M_PI;
    const double period=2*M_PI*std::sqrt(r0*r0*r0/G.GM);
    const double cruise_dt=period/cfg.steps_per_orbit;
    const double dense_dt=cruise_dt/cfg.dense_factor;
    const double max_t=cfg.days*86400.;
    const double omega=cfg.omega;
    const double envelope=cfg.dense_envelope_km*1000.;
    const bool adaptive=cfg.integrator=="adaptive";
    const double min_ad=cfg.adaptive_min_dt_s;
    double next_dt=cruise_dt, next_traj=0;
    DopriDense dense_last{};
    auto accel=[&](const Vec3& q,double t) {
        Vec3 qb=rot_z(q,-omega*t);
        return rot_z(G.acceleration(qb), omega*t) + P.acceleration(q,t);
    };
    auto energy=[&](const State& s){
        return 0.5*s.v.norm2()+G.potential_energy(rot_z(s.q,-omega*s.t));
    };
    const double initial_energy=0.5*v0.norm2()+G.potential_energy(q0);
    auto invariant=[&](const State& s) {
        return energy(s)-omega*(s.q.x*s.v.y-s.q.y*s.v.x);
    };
    auto ecc=[&](const State& s) {
        Vec3 hh=s.q.cross(s.v);
        return ((s.v.cross(hh))*(1.0/G.GM)-s.q*(1.0/s.q.norm())).norm();
    };
    // With third bodies, the lunar Jacobi value is only a reference quantity,
    // NOT a conserved integral. Report its physical variation separately.
    State s{q0,v0,accel(q0,0),0};
    if(P.enabled())R.invariant_drift_rel=std::numeric_limits<double>::quiet_NaN();
    auto kozai=[&](const State& state) {
        const double nan=std::numeric_limits<double>::quiet_NaN();
        if(!P.enabled())return std::pair<double,double>{nan,nan};
        Vec3 h=state.q.cross(state.v), n=P.reference_normal(state.t);
        double hi=h.norm();
        if(!(hi>0) || !std::isfinite(n.norm()))return std::pair<double,double>{nan,nan};
        double cosi=std::clamp(h.dot(n)/hi,-1.0,1.0);
        double e=ecc(state);
        // The classical dimensionless test-particle Kozai quantity H/L,
        // meaningful for bound Keplerian osculating orbits only.
        double jz=e<1.0?std::sqrt(std::max(0.0,1.0-e*e))*cosi:nan;
        return std::pair<double,double>{std::acos(cosi)*180./M_PI,jz};
    };
    {auto [inc,jz]=kozai(s);R.kl_incl0_deg=inc;R.kl_jz0=jz;}
    double alt_start=clearance(s,T,omega);
    double prev_z=s.q.z;
    double next_checkpoint=cfg.checkpoint_every_hours*3600.;
    double last_local_error=NAN;
    double next_progress_time=0.0;
    if(old){
        R=old->result;
        s=old->state;
        next_dt=old->next_dt;
        next_traj=old->next_traj;
        prev_z=old->prev_z;
        alt_start=old->alt_start;
        next_checkpoint=old->next_checkpoint;
        if(traj)*traj=old->history;
        if(s.t>max_t+1e-7)
            throw std::runtime_error("cannot shorten simulation below checkpoint time");
        if(old->finished && (R.impacted || s.t>=max_t-1e-7)){
            if(telemetry)telemetry(s,R,clearance(s,T,omega),energy(s),initial_energy,invariant(s),last_local_error,true,true);
            return R;
        }
        if(R.impacted)
            throw std::runtime_error("checkpoint has an unfinished but impacted orbit");
    } else {
        R.inv0=invariant(s);
        R.alt_min_seen=R.alt_max_seen=alt0_km;
        R.ecc_max=ecc(s);
    }

    auto store=[&](const State& x) {
        if (!traj) return;
        double lat,lon; body_ll(x.q,omega,x.t,lat,lon);
        traj->push_back({x.t/3600.,lat*180/M_PI,lon*180/M_PI,
                        clearance(x,T,omega)/1000.,invariant(x),ecc(x),
                        x.q.x,x.q.y,x.q.z,x.v.x,x.v.y,x.v.z,
                        kozai(x).first,kozai(x).second});
    };
    auto update_metrics=[&](const State& x,double alt) {
        R.alt_min_seen=std::min(R.alt_min_seen,alt/1000.);
        R.alt_max_seen=std::max(R.alt_max_seen,alt/1000.);
        R.ecc_max=std::max(R.ecc_max,ecc(x));
        double denom=std::max(std::fabs(R.inv0),1.0);
        const double difference=std::fabs(invariant(x)-R.inv0)/denom;
        R.lunar_jacobi_variation_rel=std::max(R.lunar_jacobi_variation_rel,difference);
        if(!P.enabled())R.invariant_drift_rel=R.lunar_jacobi_variation_rel;
    };
    auto terminate=[&](const State& impact) {
        double lat,lon; body_ll(impact.q,omega,impact.t,lat,lon);
        R.impacted=1;
        R.t_end_hr=impact.t/3600.;
        R.impact_lat=lat*180/M_PI; R.impact_lon=lon*180/M_PI;
        R.terrain_elev_m=T.surface_radius(lat,lon)-G.R;
        Vec3 b=rot_z(impact.q,-omega*impact.t);
        Vec3 ground_v=rot_z(impact.v,-omega*impact.t)-Vec3{-omega*b.y,omega*b.x,0};
        R.vx_term=ground_v.x; R.vy_term=ground_v.y; R.vz_term=ground_v.z;
        R.terminal_speed=ground_v.norm();
        Vec3 n=T.surface_normal(lat,lon);
        R.impact_angle_deg=std::asin(std::clamp(std::fabs(ground_v.dot(n))/ground_v.norm(),0.0,1.0))*180/M_PI;
        R.radial_angle_deg=std::asin(std::clamp(std::fabs(ground_v.dot(b))/(ground_v.norm()*b.norm()),0.0,1.0))*180/M_PI;
        update_metrics(impact,clearance(impact,T,omega));
        {auto [inc,jz]=kozai(impact);R.kl_incl_end_deg=inc;R.kl_jz_end=jz;}
    };

    auto progress=[&](bool terminal,bool restored=false){
        if(telemetry)telemetry(s,R,clearance(s,T,omega),energy(s),initial_energy,invariant(s),last_local_error,terminal,restored);
    };
    next_progress_time=s.t+cfg.progress_sim_hours*3600.0;
    auto checkpoint=[&](bool finished) {
        if(!save)return;
        Checkpoint snap;
        snap.state=s;snap.result=R;snap.next_dt=next_dt;
        snap.next_traj=next_traj;snap.prev_z=prev_z;snap.alt_start=alt_start;
        snap.next_checkpoint=next_checkpoint;snap.finished=finished;
        if(traj)snap.history=*traj;
        save(snap);
    };
    if (!old){
        if(traj){store(s);next_traj=cfg.traj_dt_min*60.;}
        checkpoint(false);  // initial snapshot even if killed during first orbit
    }
    // Starting on or beneath the triangulated/interpolated surface is an
    // immediate collision, not an event-search ambiguity at t=0.
    progress(false);
    if (alt_start<=1e-5) {terminate(s);checkpoint(true);progress(true);return R;}
    while (max_t-s.t > std::max(1e-8, 32*std::numeric_limits<double>::epsilon()*max_t)) {
        State start=s;
        double available=max_t-s.t;
        bool near_surface=alt_start<=envelope;
        double dt=std::min(available,near_surface?dense_dt:cruise_dt);
        if (adaptive) {
            // Soft frequency cap. The shortest lunar harmonics are attenuated
            // at high altitude, so allow longer accepted steps above 200 km.
            // Unlike the old broad terrain guard, this does not force ~2 s
            // steps for the entire near-surface orbital phase.
            const double hkm=std::max(0.0,alt_start/1000.);
            const double steps=cfg.adaptive_high_steps +
                (cfg.adaptive_low_steps-cfg.adaptive_high_steps)*
                std::exp(-hkm/cfg.adaptive_scale_height_km);
            double cap=std::min({cfg.adaptive_max_dt_s,period/steps,available});
            if (alt_start<10000.)
                cap=std::min(cap,T.equatorial_cell_width_m(s.q.norm())/
                                  (2.0*std::max(1.0,s.v.norm())));
            if (alt_start<1000.) cap=std::min(cap,2.0);
            // The final step may be shorter than the user-chosen minimum.
            const double floor=std::min(min_ad,cap);
            next_dt=std::min(next_dt,cap);
            dopri_adaptive_step(s.q,s.v,s.t,next_dt,accel,s.a,
                cfg.tol,cfg.atol_pos_m,cfg.atol_vel_mps,
                floor,cap,dense_last,&last_local_error);
            dt=s.t-start.t;
        } else if (cfg.integrator=="yoshida4") {
            yoshida4_step(s.q,s.v,s.t,dt,accel,s.a);
        } else {
            verlet_step(s.q,s.v,s.t,dt,accel,s.a);
        }
        if(!(s.t>start.t))throw std::runtime_error("time integration failed to advance");
        ++R.n_steps;
        double alt=clearance(s,T,omega);
        if (!std::isfinite(alt) || !std::isfinite(s.v.norm()))
            throw std::runtime_error("non-finite orbit state");

        // Search for the earliest terrain intersection within the accepted
        // interval. The legacy sampled method is retained for reproducibility;
        // the default adaptive mode recursively bounds the clearance between
        // samples, exposing unresolved tangencies instead of silently missing
        // one full revolution of flight. This numerical bound uses a generous
        // global bound on interpolated topographic gradients and on the speed.
        bool hit=false;
        // Use the accepted DOPRI continuous extension for all event probes:
        // it avoids seven expensive gravity evaluations per subdivision.
        auto event_state=[&](double fraction) {
            if (!adaptive) return fraction_state(start,dt*fraction,cfg,accel);
            State trial{};
            dense_last.evaluate(fraction,trial.q,trial.v);
            trial.t=start.t+fraction*dt;
            return trial;
        };
        if (cfg.event_search == "adaptive") {
            const double speed=std::max(start.v.norm(),s.v.norm()) + 20.0*dt;
            const double fastest_body=speed + std::fabs(omega)*std::max(start.q.norm(),s.q.norm());
            // Minimum-radius exclusion: a distant orbit cannot intersect any
            // terrain, irrespective of its longitude. 10m/s^2 comfortably
            // exceeds physical lunar gravity across these low orbits.
            const double distance_margin=0.5*speed*dt + 5.0*dt*dt;
            if (std::min(start.q.norm(),s.q.norm()) <= T.rmax()+distance_margin) {
                // A padded Lipschitz estimate for h=r-R(lat,lon) along the
                // trial orbit. The grid is piecewise bilinear, not C^1. Near
                // the excluded DH south pole the derivative bound is heuristic,
                // not a theorem (see numerical guide).
                const double slope=std::max(0.5,3.0*T.max_surface_slope());
                const double clearance_speed=(1.0+slope)*fastest_body;
                std::function<std::optional<double>(double,double,double,double,int)> locate;
                locate=[&](double left,double right,double hl,double hr,int depth)->std::optional<double> {
                    const double span_seconds=(right-left)*dt;
                    if (hl>0 && hr>0 && std::min(hl,hr)>0.5*span_seconds*clearance_speed)
                        return std::nullopt;  // clearance is bounded positive
                    if (span_seconds<=cfg.impact_tol_s) {
                        if (hr<=0) return right;
                        if (hl>0 && hr>0 && std::min(hl,hr)>0.5*span_seconds*clearance_speed)
                            return std::nullopt;
                    }
                    if (depth>=cfg.event_max_depth) {
                        if (hr<=0) return right;
                        // A potential sub-grid grazing contact could not be
                        // ruled out. Fail explicitly rather than report a
                        // potentially false survivor or later impact.
                        throw std::runtime_error("terrain event remains unresolved at --event-max-depth; increase resolution or depth");
                    }
                    double mid=0.5*(left+right);
                    State xm=event_state(mid);
                    double hm=clearance(xm,T,omega);
                    if (!std::isfinite(hm)) throw std::runtime_error("non-finite event clearance");
                    if (auto first=locate(left,mid,hl,hm,depth+1))return first;
                    return locate(mid,right,hm,hr,depth+1);
                };
                if (auto f=locate(0.0,1.0,alt_start,alt,0)) {
                    State impact=event_state(*f);
                    terminate(impact);s=impact;
                    if(traj)store(impact);
                    checkpoint(true);
                    hit=true;
                }
            }
        } else {
            // Original fixed-subdivision event sampling for exact comparisons.
            bool probe=(std::min(alt,alt_start) <= envelope + s.v.norm()*dt ||
                        (dt*s.v.norm() > 0.5*T.equatorial_cell_width_m(s.q.norm()) &&
                         std::min(alt,alt_start) <= 3*envelope));
            int divisions=probe?cfg.event_substeps:1;
            double prevf=0, prevh=alt_start;
            for (int k=1;k<=divisions;++k) {
                double f=static_cast<double>(k)/divisions;
                State current=k==divisions?s:event_state(f);
                double h=k==divisions?alt:clearance(current,T,omega);
                if (h<=0 && prevh>0) {
                    double left=prevf,right=f;
                    for (int it=0;it<64 && (right-left)*dt>cfg.impact_tol_s;++it) {
                        double mid=(left+right)*0.5;
                        State trial=event_state(mid);
                        if (clearance(trial,T,omega)<=0) right=mid;
                        else left=mid;
                    }
                    State impact=event_state(right);
                    terminate(impact);s=impact;
                    if (traj)store(impact);
                    checkpoint(true);hit=true;break;
                }
                prevf=f;prevh=h;
            }
        }
        if (hit) {progress(true);return R;}
        // A sampled interval can begin inside terrain only if a collision was
        // already missed; do not silently claim a valid survivor.
        if (alt<=0) throw std::runtime_error("unbracketed surface crossing");
        update_metrics(s,alt);
        if (prev_z<0 && s.q.z>=0) ++R.orbit_count;
        prev_z=s.q.z;
        alt_start=alt;
        if (traj && s.t+1e-8>=next_traj) {
            store(s);
            // Samples are at the first accepted point after the requested time.
            // Their actual timestamps appear in the trajectory CSV.
            while (next_traj<=s.t+1e-8) next_traj+=cfg.traj_dt_min*60.;
        }
        if(s.t+1e-8>=next_progress_time){
            while(next_progress_time<=s.t+1e-8)next_progress_time+=cfg.progress_sim_hours*3600.0;
            progress(false);
        }
        if(save && s.t+1e-8>=next_checkpoint){
            while(next_checkpoint<=s.t+1e-8)
                next_checkpoint+=cfg.checkpoint_every_hours*3600.;
            checkpoint(false);
        }
    }
    // Ignore the final sub-nanosecond rounding gap at the requested end time.
    s.t=max_t;
    R.t_end_hr=s.t/3600.;
    {auto [inc,jz]=kozai(s);R.kl_incl_end_deg=inc;R.kl_jz_end=jz;}
    R.terminal_speed=s.v.norm();
    checkpoint(true);
    progress(true);
    return R;
}

#include "hdf5_output.hpp"
#include "hdf5_checkpoint.hpp"

static void create_parent(const std::string& name) {
    fs::path p(name);
    if (p.has_parent_path()) fs::create_directories(p.parent_path());
}
static bool parse(int argc,char** argv,Config& c) {
    auto val=[&](int& i) {
        if (i+1>=argc) throw std::invalid_argument(std::string("missing value for ")+argv[i]);
        return std::string(argv[++i]);
    };
    for (int i=1;i<argc;++i) {
        std::string a=argv[i];
        if (a=="--grav") c.grav=val(i);
        else if (a=="--topo") c.topo=val(i);
        else if (a=="--out") c.out=val(i);
        else if (a=="--traj") c.traj=val(i);
        else if (a=="--log") c.log=val(i);
        else if (a=="--hdf5" || a=="--h5") c.hdf5=val(i);
        else if (a=="--no-csv") c.no_csv=true;
        else if (a=="--n") c.n=std::stoi(val(i));
        else if (a=="--lmax") {
            const std::string degree=val(i);
            c.lmax=(degree=="max")?-1:std::stoi(degree);
        }
        else if (a=="--alt-min") c.alt_min=std::stod(val(i));
        else if (a=="--alt-max") c.alt_max=std::stod(val(i));
        else if (a=="--days") c.days=std::stod(val(i));
        else if (a=="--steps-per-orbit") c.steps_per_orbit=std::stoi(val(i));
        else if (a=="--omega") c.omega=std::stod(val(i));
        else if (a=="--kozai-lidov") {
            const auto setting=val(i);
            if (setting=="on") c.kozai_lidov=true;
            else if (setting=="off") c.kozai_lidov=false;
            else throw std::invalid_argument("--kozai-lidov accepts on or off");
        }
        else if (a=="--enable-kozai-lidov") c.kozai_lidov=true;
        else if (a=="--disable-kozai-lidov") c.kozai_lidov=false;
        else if (a=="--third-bodies") {
            c.third_bodies=val(i);
            c.third_bodies_specified=true;
        }
        else if (a=="--ephem-mode") c.ephem_mode=val(i);
        else if (a=="--ephem-file") c.ephem_file=val(i);
        else if (a=="--dense-envelope-km") c.dense_envelope_km=std::stod(val(i));
        else if (a=="--dense-factor") c.dense_factor=std::stod(val(i));
        else if (a=="--integrator") c.integrator=val(i);
        else if (a=="--tol" || a=="--rtol") c.tol=std::stod(val(i));
        else if (a=="--atol-pos-m") c.atol_pos_m=std::stod(val(i));
        else if (a=="--atol-vel-mps") c.atol_vel_mps=std::stod(val(i));
        else if (a=="--adaptive-max-dt-s") c.adaptive_max_dt_s=std::stod(val(i));
        else if (a=="--adaptive-min-dt-s") c.adaptive_min_dt_s=std::stod(val(i));
        else if (a=="--adaptive-low-steps") c.adaptive_low_steps=std::stod(val(i));
        else if (a=="--adaptive-high-steps") c.adaptive_high_steps=std::stod(val(i));
        else if (a=="--adaptive-scale-height-km") c.adaptive_scale_height_km=std::stod(val(i));
        else if (a=="--traj-count") c.traj_count=std::stoi(val(i));
        else if (a=="--traj-dt-min") c.traj_dt_min=std::stod(val(i));
        else if (a=="--impact-tol-s") c.impact_tol_s=std::stod(val(i));
        else if (a=="--event-substeps") c.event_substeps=std::stoi(val(i));
        else if (a=="--event-search") c.event_search=val(i);
        else if (a=="--event-max-depth") c.event_max_depth=std::stoi(val(i));
        else if (a=="--seed") c.seed=static_cast<unsigned>(std::stoul(val(i)));
        else if (a=="--checkpoint-dir") c.checkpoint_dir=val(i);
        else if (a=="--checkpoint-every-hours") c.checkpoint_every_hours=std::stod(val(i));
        else if (a=="--resume") c.resume=true;
        else if (a=="--progress-dir") c.progress_dir=val(i);
        else if (a=="--progress-every-seconds") c.progress_every_seconds=std::stod(val(i));
        else if (a=="--progress-heartbeat-seconds") c.progress_heartbeat_seconds=std::stod(val(i));
        else if (a=="--progress-sim-hours") c.progress_sim_hours=std::stod(val(i));
        else if (a=="--stop-when-all-impacted") c.stop_when_all_impacted=true;
        else if (a=="--help" || a=="-h") return false;
        else throw std::invalid_argument("unknown argument: "+a);
    }
    if(!(c.progress_every_seconds>0) || !(c.progress_heartbeat_seconds>0) || !(c.progress_sim_hours>0) ||
       !std::isfinite(c.progress_every_seconds) || !std::isfinite(c.progress_heartbeat_seconds) ||
       !std::isfinite(c.progress_sim_hours))
        throw std::invalid_argument("progress cadence must be positive and finite");
    if(c.event_search!="sampled" && c.event_search!="adaptive")
        throw std::invalid_argument("--event-search must be adaptive or sampled");
    if(c.event_max_depth<1 || c.event_max_depth>36)
        throw std::invalid_argument("--event-max-depth must lie between 1 and 36");
    if (c.integrator=="yoshida") c.integrator="yoshida4"; // old script alias
    if (c.integrator!="verlet" && c.integrator!="yoshida4" && c.integrator!="adaptive")
        throw std::invalid_argument("unknown integrator: "+c.integrator);
    if (c.ephem_mode=="none" && !c.ephem_file.empty()) c.ephem_mode="file";
    // Explicit master gate: off means lunar-only, not merely hiding diagnostics.
    // Reject contradictory inputs instead of silently ignoring them.
    if (c.kozai_lidov) {
        if (!c.third_bodies_specified) c.third_bodies="all";
        const auto selected=thirdbody::select(c.third_bodies);
        if(!(selected[0] || selected[1] || selected[2]))
            throw std::invalid_argument("--kozai-lidov on needs at least one selected third body");
        // Single-switch demo mode. Requires an explicit dated SPICE file for
        // scientific interpretation; approximate phases are NOT physical epochs.
        if(c.ephem_mode=="none") c.ephem_mode="analytic";
    } else {
        if(c.third_bodies_specified && c.third_bodies!="none")
            throw std::invalid_argument("--third-bodies requires --kozai-lidov on");
        if(c.ephem_mode!="none" || !c.ephem_file.empty())
            throw std::invalid_argument("ephemeris supplied while --kozai-lidov is off");
        c.third_bodies="none";
    }
    if(c.ephem_mode!="none" && c.ephem_mode!="analytic" && c.ephem_mode!="file")
        throw std::invalid_argument("--ephem-mode must be analytic or file");
    if(c.ephem_mode=="analytic" && !c.ephem_file.empty())
        throw std::invalid_argument("--ephem-file cannot be combined with analytic mode");
    if(c.ephem_mode=="file" && c.ephem_file.empty())
        throw std::invalid_argument("file ephemeris mode needs --ephem-file FILE");
    if(c.no_csv && c.hdf5.empty())
        throw std::invalid_argument("--no-csv requires --hdf5 FILE");
#ifndef MASCONS_HAS_HDF5
    if(!c.hdf5.empty() || !c.checkpoint_dir.empty())
        throw std::invalid_argument("HDF5 output/checkpoint requested, but this build lacks HDF5; install libhdf5-dev and rebuild");
#endif
    if(c.resume && c.checkpoint_dir.empty())
        throw std::invalid_argument("--resume requires --checkpoint-dir DIR");
    if(c.checkpoint_dir.empty() && c.checkpoint_every_hours!=6)
        throw std::invalid_argument("--checkpoint-every-hours requires --checkpoint-dir DIR");
    if(!(c.checkpoint_every_hours>0) || !std::isfinite(c.checkpoint_every_hours))
        throw std::invalid_argument("checkpoint interval must be positive and finite");
    if (c.n<=0 || c.lmax<-1 || c.days<=0 || c.alt_min<0 || c.alt_max<c.alt_min ||
        c.steps_per_orbit<=0 || c.traj_count<0 || c.traj_dt_min<=0 ||
        c.dense_factor<1 || c.dense_envelope_km<0 || c.tol<=0 ||
        c.atol_pos_m<=0 || c.atol_vel_mps<=0 ||
        c.adaptive_max_dt_s<=0 || c.adaptive_min_dt_s<=0 ||
        c.adaptive_max_dt_s<c.adaptive_min_dt_s ||
        c.adaptive_low_steps<c.adaptive_high_steps || c.adaptive_high_steps<1 ||
        c.adaptive_scale_height_km<=0 ||
        !std::isfinite(c.atol_pos_m) || !std::isfinite(c.atol_vel_mps) ||
        !std::isfinite(c.adaptive_max_dt_s) || !std::isfinite(c.adaptive_min_dt_s) ||
        !std::isfinite(c.adaptive_low_steps) || !std::isfinite(c.adaptive_high_steps) ||
        !std::isfinite(c.adaptive_scale_height_km) ||
        c.impact_tol_s<=0 || c.event_substeps<1 || c.event_substeps>64 ||
        !std::isfinite(c.omega) || !std::isfinite(c.days) || !std::isfinite(c.tol) ||
        !std::isfinite(c.alt_min) || !std::isfinite(c.alt_max) ||
        !std::isfinite(c.dense_factor) || !std::isfinite(c.dense_envelope_km) ||
        !std::isfinite(c.impact_tol_s) || !std::isfinite(c.traj_dt_min))
        throw std::invalid_argument("invalid simulation parameter");
    return true;
}
static void usage() {
    std::puts("usage: lunarmc [--grav FILE][--topo FILE][--out FILE][--traj FILE][--log FILE]\n"
      " [--n N][--lmax L|max][--alt-min km][--alt-max km][--days D]\n"
      " [--steps-per-orbit S][--omega rad/s][--integrator verlet|yoshida4|adaptive]\n"
      " [--dense-envelope-km K][--dense-factor F][--rtol R]\n"
      " [--atol-pos-m A][--atol-vel-mps A][--adaptive-max-dt-s SEC]\n"
      " [--adaptive-min-dt-s SEC][--adaptive-low-steps N]\n"
      " [--adaptive-high-steps N][--adaptive-scale-height-km KM]\n"
      " [--impact-tol-s seconds][--event-substeps N]\n"
      " [--event-search adaptive|sampled][--event-max-depth 32]\n"
      " [--kozai-lidov on|off]  (default off; on activates external tides)\n"
      " [--enable-kozai-lidov|--disable-kozai-lidov] (aliases)\n"
      " [--third-bodies earth,sun,jupiter|all|none] (with switch on; default all)\n"
      " [--ephem-file moon_centered_positions.csv | --ephem-mode analytic]\n"
      " NOTE: 'on' without ephemeris selects UNDATED demonstration ephemerides.\n"
      " [--traj-count K][--traj-dt-min minutes][--seed S]\n"
      " [--checkpoint-dir DIR][--checkpoint-every-hours HOURS][--resume]\n"
      " [--stop-when-all-impacted] (explicit early-stop policy; all impacts already end the ensemble)\n"
      " [--hdf5 results/run.h5][--no-csv]  (compressed, self-describing output)\n"
      " [--progress-dir DIR][--progress-every-seconds 600][--progress-heartbeat-seconds 1800][--progress-sim-hours 6]\n"
      " Extend --days 30 to --days 60 with the SAME checkpoint directory + --resume.\n"
      " Resume requires identical physics/sampling settings and input file contents.");
}

int main(int argc,char** argv) {
  try {
    Config cfg;
    if (!parse(argc,argv,cfg)) { usage(); return 0; }
    if(cfg.lmax==-1){
        std::ifstream f(cfg.grav,std::ios::binary);
        std::int32_t n=-1;
        if(!f.read(reinterpret_cast<char*>(&n),sizeof(n)) || n<0 || n>2500)
            throw std::runtime_error("invalid gravity header when resolving --lmax max");
        cfg.lmax=int(n);
    }
    create_parent(cfg.log);
    if(!cfg.no_csv){create_parent(cfg.out);create_parent(cfg.traj);}
    if(!cfg.hdf5.empty())create_parent(cfg.hdf5);
    // The dedicated progress reporter owns run.log: stdout carries verbose diagnostics.
    auto wall0=std::chrono::steady_clock::now();
    std::time_t now=std::time(nullptr);
    char ts[64]; std::strftime(ts,sizeof(ts),"%Y-%m-%d %H:%M:%S",std::localtime(&now));
    tee("# LunarMasconsOrbit CORRECTED %s\n",ts);
    tee("config: n=%d alt=[%.2f,%.2f]km days=%.3f lmax=%d integrator=%s tol=%.2e\n",
        cfg.n,cfg.alt_min,cfg.alt_max,cfg.days,cfg.lmax,cfg.integrator.c_str(),cfg.tol);
    tee("adaptive controller: rtol=%.3e atol_pos=%.3g m atol_vel=%.3g m/s dt=[%.3g,%.3g]s\n",
        cfg.tol,cfg.atol_pos_m,cfg.atol_vel_mps,cfg.adaptive_min_dt_s,cfg.adaptive_max_dt_s);
    tee("adaptive soft cap: steps/orbit low=%.0f high=%.0f height_scale=%.1f km\n",
        cfg.adaptive_low_steps,cfg.adaptive_high_steps,cfg.adaptive_scale_height_km);
    tee("omega=%.9e steps_per_orbit=%d dense_factor=%.1f envelope=%.1fkm\n",
        cfg.omega,cfg.steps_per_orbit,cfg.dense_factor,cfg.dense_envelope_km);
    tee("impact_tol_s=%.3g event_substeps=%d search=%s max_depth=%d seed=%u\n",
        cfg.impact_tol_s,cfg.event_substeps,cfg.event_search.c_str(),cfg.event_max_depth,cfg.seed);
    const thirdbody::Perturbations P(cfg.third_bodies,cfg.ephem_mode,cfg.ephem_file);
    tee("kozai_lidov=%s (external-tide switch); third_bodies=%s ephemeris=%s %s\n",
        cfg.kozai_lidov?"on":"off",cfg.third_bodies.c_str(),
        cfg.ephem_mode.c_str(),cfg.ephem_file.c_str());
    if(P.enabled()) {
        tee("# WARNING: lunar Jacobi is not conserved with third bodies; "
            "invariant_drift_rel=nan; lunar_jacobi_variation_rel is physical + numerical.\n");
        if(P.mode()=="analytic")tee("# WARNING: approximate arbitrary-phase Kepler ephemeris; NOT for research claims.\n");
        else if(P.table_start()>-2*86400 || P.table_end()<cfg.days*86400+2*86400)
            tee("# Ephemeris edge caution: ensure coverage of Yoshida internal stage times.\n");
    }
    const std::string signature=cfg.checkpoint_dir.empty()?"":run_fingerprint(cfg);
    prepare_checkpoint_dir(cfg,signature);
    if(!cfg.checkpoint_dir.empty())
        tee("checkpoints: %s mode=%s every=%.6g simulated hours signature=%s\n",
            cfg.checkpoint_dir.c_str(),cfg.resume?"RESUME":"NEW",
            cfg.checkpoint_every_hours,signature.c_str());
    SphGravity G; G.load(cfg.grav,cfg.lmax);
    if (G.lmax!=cfg.lmax) throw std::runtime_error("requested lmax absent from gravity file");
    Topography T; T.load(cfg.topo);
    tee("gravity: lmax=%d GM=%.9e R=%.1f m | topo: %dx%d [%g,%g] km\n",
        G.lmax,G.GM,G.R,T.nlat,T.nlon,T.rmin()/1000.,T.rmax()/1000.);
    // Internal consistency test only; independent normalization checks in CTest.
    std::mt19937 test_rng(1);std::uniform_real_distribution<double> unit(0,1);
    double maxrel=0;
    for(int k=0;k<50;++k){
        double lat=std::asin(2*unit(test_rng)-1),lon=2*M_PI*unit(test_rng);
        double rr=G.R+(25+225*unit(test_rng))*1000;
        Vec3 p=sphere_pos(lat,lon,rr);
        maxrel=std::max(maxrel,(G.acceleration(p)-G.acceleration_fd(p)).norm()/G.acceleration_fd(p).norm());
    }
    tee("gradient consistency: maximum relative error %.3e %s\n",maxrel,
        maxrel<1e-5?"PASS":"FAIL");
    if(maxrel>=1e-5) throw std::runtime_error("gravity analytic-gradient check failed");
    std::mt19937 rng(cfg.seed);std::uniform_real_distribution<double> U(0,1);
    std::vector<double> lat0(cfg.n),lon0(cfg.n),alt0(cfg.n),head(cfg.n);
    for(int i=0;i<cfg.n;++i){
        lat0[i]=std::asin(2*U(rng)-1);lon0[i]=2*M_PI*U(rng);
        alt0[i]=cfg.alt_min+(cfg.alt_max-cfg.alt_min)*U(rng);
        head[i]=2*M_PI*U(rng);
    }
    std::vector<Result> res(cfg.n);
    std::vector<std::vector<TrajPoint>> trajs(cfg.n);
    std::vector<std::string> errors(cfg.n);
    ProgressReporter progress_reporter(cfg,P.enabled());
    tee("progress: poll=%.3g s heartbeat=%.3g s per-particle=%.3g simulated hours dir=%s\n",
        cfg.progress_every_seconds,cfg.progress_heartbeat_seconds,cfg.progress_sim_hours,
        (cfg.progress_dir.empty()?(fs::path(cfg.log).has_parent_path()?fs::path(cfg.log).parent_path().string():"."):cfg.progress_dir).c_str());
    progress_reporter.start();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for(int i=0;i<cfg.n;++i){
        try {
            std::vector<TrajPoint>* tp=(i<cfg.traj_count)?&trajs[i]:nullptr;
            std::optional<Checkpoint> previous;
            std::function<void(const Checkpoint&)> save;
            if(!cfg.checkpoint_dir.empty()){
                const auto path=checkpoint_path(cfg,i);
                if(cfg.resume && (fs::exists(path) || fs::exists(path.substr(0,path.size()-3)+".prev.h5"))){
                    #ifdef MASCONS_HAS_HDF5
                    previous=lunarh5::load_snapshot_pair(path,i,signature);
#endif
                    // Never accept an accidental checkpoint for a different orbit.
                    if(std::fabs(previous->result.lat0-lat0[i]*180/M_PI)>1e-9 ||
                       std::fabs(previous->result.lon0-lon0[i]*180/M_PI)>1e-9 ||
                       std::fabs(previous->result.alt0-alt0[i])>1e-9 ||
                       std::fabs(previous->result.heading-head[i]*180/M_PI)>1e-9)
                        throw std::runtime_error("checkpoint initial condition mismatch: "+path);
                }
                #ifdef MASCONS_HAS_HDF5
                save=[path,i,signature](const Checkpoint& cp){lunarh5::save_snapshot_pair(path,i,signature,cp);};
#endif
            }
            auto observe=[&progress_reporter,i](const State& s,const Result& r,double altitude,
                double energy,double initial_energy,double J,double local_error,bool terminal,bool restored){
                progress_reporter.update(i,s,r,altitude,energy,initial_energy,J,local_error,terminal,restored);
            };
            res[i]=propagate(G,T,P,cfg,lat0[i],lon0[i],alt0[i],head[i],tp,
                             previous?&*previous:nullptr,save,observe);
        } catch(const std::exception& e) { errors[i]=e.what();progress_reporter.error(i,e.what()); }
    }
    progress_reporter.finish();
    for(int i=0;i<cfg.n;++i) if(!errors[i].empty()) {
        throw std::runtime_error("trajectory "+std::to_string(i)+": "+errors[i]);
    }
    int impacts=0;long long steps=0;double drift=0;
    for(const auto& r:res){
        impacts+=r.impacted;steps+=r.n_steps;
        if(std::isfinite(r.invariant_drift_rel))drift=std::max(drift,r.invariant_drift_rel);
    }
    if(!cfg.no_csv){
    FILE* f=std::fopen(cfg.out.c_str(),"w");
    if(!f) throw std::runtime_error("cannot write output: "+cfg.out);
    std::fprintf(f,"id,lat0,lon0,alt0_km,heading,inc0,impacted,t_end_hr,impact_lat,impact_lon,"
                   "impact_angle_deg,terminal_speed_ms,vx_term,vy_term,vz_term,terrain_elev_m,"
                   "alt_min_km,alt_max_km,ecc_max,orbit_count,inv0,invariant_drift_rel,"
                   "radial_angle_deg,n_steps,lunar_jacobi_variation_rel,kl_incl0_deg,"
                   "kl_incl_end_deg,kl_jz0,kl_jz_end\n");
    for(int i=0;i<cfg.n;++i){
        const auto& r=res[i];
        std::fprintf(f,"%d,%.6f,%.6f,%.4f,%.4f,%.4f,%d,%.9f,%.6f,%.6f,%.6f,"
            "%.6f,%.6f,%.6f,%.6f,%.3f,%.6f,%.6f,%.8f,%d,%.12e,%.6e,%.6f,%ld,"
            "%.6e,%.8f,%.8f,%.9f,%.9f\n",
            i,r.lat0,r.lon0,r.alt0,r.heading,r.inc0,r.impacted,r.t_end_hr,
            r.impact_lat,r.impact_lon,r.impact_angle_deg,r.terminal_speed,
            r.vx_term,r.vy_term,r.vz_term,r.terrain_elev_m,r.alt_min_seen,
            r.alt_max_seen,r.ecc_max,r.orbit_count,r.inv0,r.invariant_drift_rel,
            r.radial_angle_deg,r.n_steps,r.lunar_jacobi_variation_rel,
            r.kl_incl0_deg,r.kl_incl_end_deg,r.kl_jz0,r.kl_jz_end);
    }
    std::fclose(f);
    FILE* g=std::fopen(cfg.traj.c_str(),"w");
    if(!g) throw std::runtime_error("cannot write trajectories: "+cfg.traj);
    std::fprintf(g,"id,t_hr,lat,lon,alt_km,invariant,ecc,x,y,z,vx,vy,vz,"
                   "kl_incl_deg,kl_jz\n");
    for(int i=0;i<std::min(cfg.n,cfg.traj_count);++i)
        for(const auto& p:trajs[i])
            std::fprintf(g,"%d,%.9f,%.6f,%.6f,%.6f,%.12e,%.9f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,"
                "%.8f,%.9f\n",
                i,p.t_hr,p.lat,p.lon,p.alt_km,p.invariant,p.ecc,
                p.x,p.y,p.z,p.vx,p.vy,p.vz,p.kl_incl_deg,p.kl_jz);
    std::fclose(g);
    } // optional legacy CSV export
#ifdef MASCONS_HAS_HDF5
    if(!cfg.hdf5.empty()){
        lunarh5::write(cfg.hdf5,cfg,res,trajs);
        tee("wrote HDF5: %s\n",cfg.hdf5.c_str());
    }
#endif
    double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-wall0).count();
    if(P.enabled()) {
        double variation=0;
        for(const auto& r:res) variation=std::max(variation,r.lunar_jacobi_variation_rel);
        tee("impacts: %d/%d (%.2f%%); lunar-only Jacobi variation=%.3e "
            "(NOT an integration-error estimate)\n",impacts,cfg.n,100.*impacts/cfg.n,variation);
    } else {
        tee("impacts: %d/%d (%.2f%%) worst %s drift=%.3e\n",impacts,cfg.n,
            100.*impacts/cfg.n,cfg.omega==0?"energy":"Jacobi",drift);
    }
    // Propagation is particle-wise: each impacted orbit stops immediately. If all
    // particles impact, no integrations remain before the requested horizon;
    // the option makes this termination policy explicit in the recorded log.
    if(cfg.stop_when_all_impacted && impacts==cfg.n){
        double last_impact_hr=0.;
        for(const auto& r:res) last_impact_hr=std::max(last_impact_hr,r.t_end_hr);
        tee("ensemble stop: all %d particles impacted; last impact at %.9f hours "
            "(%.9f simulated days)\n",cfg.n,last_impact_hr,last_impact_hr/24.);
    } else if(cfg.stop_when_all_impacted){
        tee("ensemble stop: requested %.6f days reached; %d/%d particles still active\n",
            cfg.days,cfg.n-impacts,cfg.n);
    }
    tee("steps: total=%lld average=%.1f walltime=%.2fs\n",steps,double(steps)/cfg.n,wall);
    if(cfg.n==1){
        const auto& r=res[0];
        tee("single: lat0=%.3f lon0=%.3f alt0=%.1fkm inc0=%.2f impacted=%d "
            "t_end=%.6fhr impact=(%.3f,%.3f) angle_local=%.5f angle_radial=%.5f "
            "vterm=%.2f orbits=%d steps=%ld\n",
            r.lat0,r.lon0,r.alt0,r.inc0,r.impacted,r.t_end_hr,r.impact_lat,
            r.impact_lon,r.impact_angle_deg,r.radial_angle_deg,r.terminal_speed,
            r.orbit_count,r.n_steps);
    }
    if(!cfg.no_csv)tee("wrote CSV: %s %s\n",cfg.out.c_str(),cfg.traj.c_str());
    return 0;
  } catch(const std::exception& e) {
    std::fprintf(stderr,"ERROR: %s\n",e.what());
    return 2;
  }
}
