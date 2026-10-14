#pragma once
// Version 2: self-describing, native HDF5 per-particle restart snapshots.
// Included after hdf5_output.hpp in main.cpp. One HDF5 file is written per
// particle, atomically rotated with the previous valid generation. The final
// analysis HDF5 is a separate archival product, produced at run completion.
#ifdef MASCONS_HAS_HDF5
#include <mutex>

namespace lunarh5 {
inline std::mutex& checkpoint_mutex(){static std::mutex m;return m;}
inline std::string read_string_attr(hid_t obj,const char* name){
    Handle a(H5Aopen(obj,name,H5P_DEFAULT),H5Aclose);
    Handle t(H5Aget_type(a),H5Tclose);
    const std::size_t n=H5Tget_size(t);
    if(!n || n>4096 || H5Tget_class(t)!=H5T_STRING || H5Tis_variable_str(t))
        throw std::runtime_error("invalid HDF5 checkpoint text attribute");
    std::vector<char> b(n+1,0);
    check(H5Aread(a,t,b.data()),"reading string attribute");
    if(b[n-1]!=0) throw std::runtime_error("unterminated HDF5 checkpoint attribute");
    return std::string(b.data());
}
inline int read_int_attr(hid_t obj,const char* name){
    Handle a(H5Aopen(obj,name,H5P_DEFAULT),H5Aclose);int result=0;
    check(H5Aread(a,H5T_NATIVE_INT,&result),"reading int attribute");return result;
}
inline void read_data(hid_t file,const char* name,hid_t type,hsize_t expected,void* target){
    Handle d(H5Dopen2(file,name,H5P_DEFAULT),H5Dclose);
    Handle space(H5Dget_space(d),H5Sclose);
    hsize_t dims[1]={0};
    if(H5Sget_simple_extent_ndims(space)!=1||H5Sget_simple_extent_dims(space,dims,nullptr)<0||dims[0]!=expected)
        throw std::runtime_error("invalid checkpoint dataset shape: "+std::string(name));
    if(expected)check(H5Dread(d,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,target),"reading "+std::string(name));
}
inline void write_snapshot_file(const std::string& path,int id,const std::string& signature,const Checkpoint& cp){
    // Canonical checkpoint hash protects all restart-critical fields, not just
    // HDF5 metadata. Hash recomputed from decoded native datasets on resume.
    const std::string checksum=hex64(fnv64(checkpoint_payload(id,signature,cp)));
    Handle f(H5Fcreate(path.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT),H5Fclose);
    attribute_string(f,"schema","lunar_particle_restart_hdf5_v2");
    attribute_string(f,"run_signature",signature);
    attribute_string(f,"state_checksum_fnv64",checksum);
    attribute_int(f,"particle_id",id);
    attribute_int(f,"finished",cp.finished?1:0);
    attribute_double(f,"simulation_time_s",cp.state.t);
    attribute_int(f,"history_count",static_cast<int>(cp.history.size()));
    const double s[10]={cp.state.t,cp.state.q.x,cp.state.q.y,cp.state.q.z,
        cp.state.v.x,cp.state.v.y,cp.state.v.z,
        cp.state.a.x,cp.state.a.y,cp.state.a.z};
    const double ctrl[5]={cp.next_dt,cp.next_traj,cp.prev_z,cp.alt_start,cp.next_checkpoint};
    dataset(f,"state_t_q_v_a",H5T_NATIVE_DOUBLE,10,s,false);
    dataset(f,"integrator_control",H5T_NATIVE_DOUBLE,5,ctrl,false);
    const SummaryRow summary{id,cp.result};auto stype=summary_type();
    dataset(f,"summary",stype,1,&summary,false);
    auto dtype=trajectory_type();
    dataset(f,"history",dtype,cp.history.size(),cp.history.data());
    check(H5Fflush(f,H5F_SCOPE_GLOBAL),"flushing checkpoint");
}
inline Checkpoint read_snapshot_file(const std::string& path,int id,const std::string& signature){
    // Corrupt HDF5 files may emit error stacks; exceptions are handled by caller.
    Handle f(H5Fopen(path.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT),H5Fclose);
    if(read_string_attr(f,"schema")!="lunar_particle_restart_hdf5_v2" ||
       read_string_attr(f,"run_signature")!=signature || read_int_attr(f,"particle_id")!=id)
        throw std::runtime_error("checkpoint run signature/id/schema mismatch: "+path);
    int finished=read_int_attr(f,"finished");
    if(finished!=0 && finished!=1)throw std::runtime_error("invalid checkpoint status");
    const int n=read_int_attr(f,"history_count");
    if(n<0 || n>10000000)throw std::runtime_error("invalid checkpoint trajectory length");
    Checkpoint cp;cp.finished=finished;
    double s[10],ctrl[5];
    read_data(f,"state_t_q_v_a",H5T_NATIVE_DOUBLE,10,s);
    read_data(f,"integrator_control",H5T_NATIVE_DOUBLE,5,ctrl);
    cp.state={Vec3{s[1],s[2],s[3]},Vec3{s[4],s[5],s[6]},Vec3{s[7],s[8],s[9]},s[0]};
    cp.next_dt=ctrl[0];cp.next_traj=ctrl[1];cp.prev_z=ctrl[2];
    cp.alt_start=ctrl[3];cp.next_checkpoint=ctrl[4];
    SummaryRow summary{};auto stype=summary_type();
    read_data(f,"summary",stype,1,&summary);
    if(summary.id!=id)throw std::runtime_error("mismatched checkpoint summary id");
    cp.result=summary.result;
    cp.history.resize(n);auto dtype=trajectory_type();
    read_data(f,"history",dtype,cp.history.size(),cp.history.data());
    if(cp.state.t<0 || !std::isfinite(cp.state.t)||!(cp.next_dt>0)||
       !(cp.next_checkpoint>0)||!std::isfinite(cp.alt_start))
        throw std::runtime_error("invalid/non-finite restart state");
    const std::string hash=hex64(fnv64(checkpoint_payload(id,signature,cp)));
    if(hash!=read_string_attr(f,"state_checksum_fnv64"))
        throw std::runtime_error("checkpoint checksum mismatch: "+path);
    return cp;
}
inline std::string prior_path(const std::string& latest){
    if(latest.size()<3 || latest.substr(latest.size()-3)!=".h5")
        throw std::runtime_error("checkpoint path must end with .h5");
    return latest.substr(0,latest.size()-3)+".prev.h5";
}
inline void save_snapshot_pair(const std::string& latest,int id,const std::string& signature,const Checkpoint& cp){
    std::lock_guard<std::mutex> lock(checkpoint_mutex());
    const auto older=prior_path(latest), temp=latest+".tmp";
    try {
        write_snapshot_file(temp,id,signature,cp);
        // Verify staged snapshot before risking either of the saved generations.
        (void)read_snapshot_file(temp,id,signature);
        if(fs::exists(latest)){
            bool latest_valid=false;
            try{(void)read_snapshot_file(latest,id,signature);latest_valid=true;}
            catch(const std::exception&){ /* keep existing previous generation */ }
            if(latest_valid){
                // POSIX rename atomically replaces the prior older generation.
                fs::rename(latest,older);
            }else{
                // A bad latest must never overwrite the intact fallback.
                fs::remove(latest);
            }
        }
        fs::rename(temp,latest);
    }catch(...){std::error_code ec;fs::remove(temp,ec);throw;}
}
inline Checkpoint load_snapshot_pair(const std::string& latest,int id,const std::string& signature){
    std::lock_guard<std::mutex> lock(checkpoint_mutex());
    const auto older=prior_path(latest);
    std::string errors;
    for(const auto& path:{latest,older}){
        if(!fs::exists(path))continue;
        try{
            Checkpoint cp=read_snapshot_file(path,id,signature);
            if(path==older){
                std::fprintf(stderr,"WARNING: particle %d using previous HDF5 checkpoint: %s\n",id,path.c_str());
            }
            return cp;
        }
        catch(const std::exception& e){errors+=path+": "+e.what()+"; ";}
    }
    if(errors.empty())throw std::runtime_error("missing both checkpoint generations for particle "+std::to_string(id));
    throw std::runtime_error("both checkpoint generations unavailable/corrupt: "+errors);
}
} // namespace lunarh5
#endif
