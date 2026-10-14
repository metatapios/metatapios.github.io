#pragma once
// HDF5 schema v1 for portable, self-describing simulation results.
// Include this header from main.cpp after Config, Result and TrajPoint.
// The HDF5 writer runs serially after the parallel integration. Checkpoints
// are separate HDF5 files per orbit, with an intact rolling backup generation.

#ifdef MASCONS_HAS_HDF5
#include <hdf5.h>
#include <cstddef>
#include <type_traits>

namespace lunarh5 {
struct Handle {
    hid_t id = -1;
    herr_t (*close)(hid_t) = nullptr;
    Handle() = default;
    Handle(hid_t h, herr_t(*closer)(hid_t)) : id(h), close(closer) {
        if(id < 0) throw std::runtime_error("HDF5 failed to allocate object");
    }
    ~Handle(){ if(id >= 0) close(id); }
    Handle(const Handle&) = delete;
    Handle& operator=(const Handle&) = delete;
    Handle(Handle&& other) noexcept : id(other.id), close(other.close) {other.id=-1;}
    Handle& operator=(Handle&& other) noexcept {
        if(this!=&other){if(id>=0)close(id);id=other.id;close=other.close;other.id=-1;}
        return *this;
    }
    operator hid_t() const { return id; }
};
inline void check(herr_t error,const std::string& context){
    if(error<0) throw std::runtime_error("HDF5 "+context);
}
inline void attribute_string(hid_t obj,const char* name,const std::string& value){
    Handle type(H5Tcopy(H5T_C_S1),H5Tclose);
    check(H5Tset_size(type,value.size()+1),"setting text attribute size");
    check(H5Tset_strpad(type,H5T_STR_NULLTERM),"setting string padding");
    Handle space(H5Screate(H5S_SCALAR),H5Sclose);
    Handle attr(H5Acreate2(obj,name,type,space,H5P_DEFAULT,H5P_DEFAULT),H5Aclose);
    check(H5Awrite(attr,type,value.c_str()),"writing attribute "+std::string(name));
}
inline void attribute_double(hid_t obj,const char* name,double value){
    Handle space(H5Screate(H5S_SCALAR),H5Sclose);
    Handle attr(H5Acreate2(obj,name,H5T_NATIVE_DOUBLE,space,H5P_DEFAULT,H5P_DEFAULT),H5Aclose);
    check(H5Awrite(attr,H5T_NATIVE_DOUBLE,&value),"writing attribute "+std::string(name));
}
inline void attribute_int(hid_t obj,const char* name,int value){
    Handle space(H5Screate(H5S_SCALAR),H5Sclose);
    Handle attr(H5Acreate2(obj,name,H5T_NATIVE_INT,space,H5P_DEFAULT,H5P_DEFAULT),H5Aclose);
    check(H5Awrite(attr,H5T_NATIVE_INT,&value),"writing attribute "+std::string(name));
}

struct SummaryRow { int id; Result result; };
static_assert(std::is_standard_layout<SummaryRow>::value,"HDF5 rows must be standard layout");
static_assert(std::is_standard_layout<TrajPoint>::value,"HDF5 rows must be standard layout");
inline Handle summary_type(){
    Handle t(H5Tcreate(H5T_COMPOUND,sizeof(SummaryRow)),H5Tclose);
    check(H5Tinsert(t,"id",offsetof(SummaryRow,id),H5T_NATIVE_INT),"compound id");
#define SM(name, member, ty) check(H5Tinsert(t,name,offsetof(SummaryRow,result)+offsetof(Result,member),ty),"compound " name)
    SM("lat0_deg",lat0,H5T_NATIVE_DOUBLE); SM("lon0_deg",lon0,H5T_NATIVE_DOUBLE);
    SM("alt0_km",alt0,H5T_NATIVE_DOUBLE); SM("heading_deg",heading,H5T_NATIVE_DOUBLE);
    SM("inc0_deg",inc0,H5T_NATIVE_DOUBLE); SM("impacted",impacted,H5T_NATIVE_INT);
    SM("t_end_hr",t_end_hr,H5T_NATIVE_DOUBLE);
    SM("impact_lat_deg",impact_lat,H5T_NATIVE_DOUBLE); SM("impact_lon_deg",impact_lon,H5T_NATIVE_DOUBLE);
    SM("impact_angle_deg",impact_angle_deg,H5T_NATIVE_DOUBLE);
    SM("radial_angle_deg",radial_angle_deg,H5T_NATIVE_DOUBLE);
    SM("terminal_speed_ms",terminal_speed,H5T_NATIVE_DOUBLE);
    SM("terrain_elev_m",terrain_elev_m,H5T_NATIVE_DOUBLE);
    SM("vx_term_ms",vx_term,H5T_NATIVE_DOUBLE);SM("vy_term_ms",vy_term,H5T_NATIVE_DOUBLE);
    SM("vz_term_ms",vz_term,H5T_NATIVE_DOUBLE);
    SM("alt_min_km",alt_min_seen,H5T_NATIVE_DOUBLE);
    SM("alt_max_km",alt_max_seen,H5T_NATIVE_DOUBLE);
    SM("ecc_max",ecc_max,H5T_NATIVE_DOUBLE);
    SM("orbit_count",orbit_count,H5T_NATIVE_INT);
    SM("n_steps",n_steps,H5T_NATIVE_LONG);
    SM("inv0_si",inv0,H5T_NATIVE_DOUBLE);
    SM("invariant_drift_rel",invariant_drift_rel,H5T_NATIVE_DOUBLE);
    SM("lunar_jacobi_variation_rel",lunar_jacobi_variation_rel,H5T_NATIVE_DOUBLE);
    SM("kl_incl0_deg",kl_incl0_deg,H5T_NATIVE_DOUBLE);
    SM("kl_incl_end_deg",kl_incl_end_deg,H5T_NATIVE_DOUBLE);
    SM("kl_jz0",kl_jz0,H5T_NATIVE_DOUBLE);
    SM("kl_jz_end",kl_jz_end,H5T_NATIVE_DOUBLE);
#undef SM
    return t;
}
inline Handle trajectory_type(){
    Handle t(H5Tcreate(H5T_COMPOUND,sizeof(TrajPoint)),H5Tclose);
#define TP(name, member) check(H5Tinsert(t,name,offsetof(TrajPoint,member),H5T_NATIVE_DOUBLE),"compound " name)
    TP("t_hr",t_hr); TP("lat_deg",lat); TP("lon_deg",lon);
    TP("alt_km",alt_km); TP("jacobi_si",invariant); TP("ecc",ecc);
    TP("x_m",x);TP("y_m",y);TP("z_m",z);
    TP("vx_ms",vx);TP("vy_ms",vy);TP("vz_ms",vz);
    TP("kl_incl_deg",kl_incl_deg);TP("kl_jz",kl_jz);
#undef TP
    return t;
}
inline void dataset(hid_t parent,const char* name,hid_t type,hsize_t count,const void* data,bool compress=true){
    const hsize_t dims[1]={count};
    Handle space(H5Screate_simple(1,dims,nullptr),H5Sclose);
    Handle plist(H5Pcreate(H5P_DATASET_CREATE),H5Pclose);
    if(compress){
        const hsize_t chunk[1]={std::min<hsize_t>(4096,std::max<hsize_t>(1,count))};
        check(H5Pset_chunk(plist,1,chunk),"setting chunks");
        check(H5Pset_shuffle(plist),"setting byte shuffle");
        check(H5Pset_deflate(plist,4),"enabling gzip");
    }
    Handle d(H5Dcreate2(parent,name,type,space,H5P_DEFAULT,plist,H5P_DEFAULT),H5Dclose);
    if(count)check(H5Dwrite(d,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,data),"writing "+std::string(name));
}

inline void write(const std::string& dest,const Config& c,
                  const std::vector<Result>& summaries,
                  const std::vector<std::vector<TrajPoint>>& trajectories){
    // Write to a sibling staging file so an existing result is never replaced
    // by an interrupted HDF5 write. Rename only after H5Fclose succeeds.
    const std::string staging=dest+".tmp";
    try {
        {
            Handle f(H5Fcreate(staging.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT),H5Fclose);
            attribute_string(f,"schema","lunar_orbits_hdf5_v1");
            attribute_string(f,"description","Rotating lunar harmonics, Monte Carlo orbit integration");
            attribute_string(f,"coordinate_frames","trajectory x/y/z and v are lunar-centered inertial; trajectory lat/lon are Moon-body-fixed; longitudes east-positive");
            attribute_string(f,"units","SI unless field name gives km, hr, deg; specific Jacobi integral has units m^2 s^-2");
            attribute_string(f,"impact_status","summary.impacted=1 denotes a permanent surface impact; impacted orbits never restart");
            Handle meta(H5Gcreate2(f,"metadata",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Gclose);
            attribute_int(meta,"n_particles",c.n);
            attribute_int(meta,"gravity_lmax",c.lmax);
            attribute_int(meta,"steps_per_orbit",c.steps_per_orbit);
            attribute_int(meta,"traj_count",c.traj_count);
            attribute_int(meta,"seed",static_cast<int>(c.seed));
            attribute_int(meta,"event_substeps",c.event_substeps);
            attribute_int(meta,"event_max_depth",c.event_max_depth);
            attribute_string(meta,"event_search",c.event_search);
            attribute_double(meta,"days",c.days);
            attribute_double(meta,"alt_min_km",c.alt_min);
            attribute_double(meta,"alt_max_km",c.alt_max);
            attribute_double(meta,"omega_rad_s",c.omega);
            attribute_double(meta,"traj_dt_min",c.traj_dt_min);
            attribute_double(meta,"dense_envelope_km",c.dense_envelope_km);
            attribute_double(meta,"dense_factor",c.dense_factor);
            attribute_double(meta,"integrator_tolerance",c.tol);
            attribute_double(meta,"adaptive_rtol",c.tol);
            attribute_double(meta,"adaptive_atol_position_m",c.atol_pos_m);
            attribute_double(meta,"adaptive_atol_velocity_ms",c.atol_vel_mps);
            attribute_double(meta,"adaptive_min_step_s",c.adaptive_min_dt_s);
            attribute_double(meta,"adaptive_max_step_s",c.adaptive_max_dt_s);
            attribute_double(meta,"adaptive_low_alt_steps_per_orbit",c.adaptive_low_steps);
            attribute_double(meta,"adaptive_high_alt_steps_per_orbit",c.adaptive_high_steps);
            attribute_double(meta,"adaptive_scale_height_km",c.adaptive_scale_height_km);
            attribute_double(meta,"impact_tol_s",c.impact_tol_s);
            attribute_string(meta,"integrator",c.integrator);
            attribute_string(meta,"kozai_lidov",c.kozai_lidov?"on":"off");
            attribute_string(meta,"third_bodies",c.third_bodies);
            attribute_string(meta,"ephemeris_mode",c.ephem_mode);
            attribute_string(meta,"gravity_path",c.grav);
            attribute_string(meta,"topography_path",c.topo);
            attribute_string(meta,"gravity_fnv64",hex64(file_hash(c.grav)));
            attribute_string(meta,"topography_fnv64",hex64(file_hash(c.topo)));
            attribute_string(meta,"run_signature",run_fingerprint(c));
            attribute_string(meta,"checkpoint_directory",c.checkpoint_dir);
            attribute_string(meta,"progress_directory",c.progress_dir.empty()?
                (fs::path(c.log).has_parent_path()?fs::path(c.log).parent_path().string():"."):
                c.progress_dir);
            attribute_double(meta,"progress_wall_interval_s",c.progress_every_seconds);
            attribute_double(meta,"progress_particle_sim_interval_hr",c.progress_sim_hours);
            attribute_string(meta,"restart_mode",c.resume?"resumed":"new");
            std::vector<SummaryRow> summary_rows;
            summary_rows.reserve(summaries.size());
            for(size_t i=0;i<summaries.size();++i)
                summary_rows.push_back({static_cast<int>(i),summaries[i]});
            auto stype=summary_type();
            dataset(f,"summary",stype,summary_rows.size(),summary_rows.data());
            Handle tg(H5Gcreate2(f,"trajectories",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Gclose);
            std::vector<std::uint64_t> offsets(trajectories.size()+1,0);
            for(size_t i=0;i<trajectories.size();++i)
                offsets[i+1]=offsets[i]+trajectories[i].size();
            std::vector<TrajPoint> points;
            points.reserve(static_cast<size_t>(offsets.back()));
            for(const auto& tr:trajectories)points.insert(points.end(),tr.begin(),tr.end());
            auto ttype=trajectory_type();
            dataset(tg,"samples",ttype,points.size(),points.data());
            dataset(tg,"offsets",H5T_NATIVE_UINT64,offsets.size(),offsets.data());
            attribute_string(tg,"layout","Ragged array: samples[offsets[i]:offsets[i+1]] are particle i; empty for unsampled particles");
            attribute_string(tg,"sampling","Samples at the first integrator step following requested output time, plus initial and impact samples");
            check(H5Fflush(f,H5F_SCOPE_GLOBAL),"flushing simulation output");
        }
        fs::rename(staging,dest);
    } catch(...) {std::error_code ec;fs::remove(staging,ec);throw;}
}
} // namespace lunarh5
#endif
