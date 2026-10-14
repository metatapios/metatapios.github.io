#pragma once
#include "progress_schedule.hpp"
// Thread-safe, low-frequency telemetry for independent lunar orbit workers.
// No HDF5 calls occur here; final scientific HDF5 output remains consolidated.
#include <condition_variable>
#include <mutex>
#include <thread>

struct ParticleTelemetry {
    enum Status : int {pending=0, running=1, impacted=2, survived=3, failed=4};
    Status status=pending;
    double t_days=0.0, altitude_km=NAN, eccentricity=NAN;
    double energy_change_rel=NAN, jacobi_change_rel=NAN, local_error_ratio=NAN;
    double max_lunar_jacobi_variation_rel=NAN;
    double numeric_invariant_drift_rel=NAN;
    long n_steps=0;
    double impact_lat_deg=NAN, impact_lon_deg=NAN, terminal_speed_mps=NAN;
    std::string error;
};

class ProgressReporter {
    const Config& cfg_;
    const bool perturbed_;
    const std::string root_;
    std::vector<ParticleTelemetry> particles_;
    std::ofstream events_, timeline_, simple_log_;
    int newly_completed_session_=0;
    int last_logged_completed_=-1;
    std::chrono::steady_clock::time_point last_log_time_{};
    std::mutex mutex_;
    std::condition_variable wake_;
    bool stop_=false;
    std::thread monitor_;
    const std::chrono::steady_clock::time_point started_=std::chrono::steady_clock::now();
    static const char* label(ParticleTelemetry::Status s) {
        switch(s){case ParticleTelemetry::pending:return "pending";
        case ParticleTelemetry::running:return "running";
        case ParticleTelemetry::impacted:return "impacted";
        case ParticleTelemetry::survived:return "survived";
        case ParticleTelemetry::failed:return "failed";}
        return "unknown";
    }
    static std::string iso_time() {
        std::time_t now=std::time(nullptr);
        std::tm tm{};
#if defined(_WIN32)
        localtime_s(&tm,&now);
#else
        localtime_r(&now,&tm);
#endif
        char buf[40];std::strftime(buf,sizeof(buf),"%Y-%m-%d %H:%M:%S",&tm);
        return buf;
    }
    // ETA is calculated from *new* completions this process has produced,
    // not restored terminal checkpoints. It is deliberately provisional: early
    // collisions and 25--500 km lifetime heterogeneity can bias this rate.
    static std::string eta_timestamp(double remaining_seconds) {
        if (!std::isfinite(remaining_seconds) || remaining_seconds<0 || remaining_seconds>1.0e11)
            return "unknown";
        const auto now=std::time(nullptr);
        auto eta=now+static_cast<std::time_t>(std::ceil(remaining_seconds));
        std::tm tm{};
#if defined(_WIN32)
        localtime_s(&tm,&eta);
#else
        localtime_r(&eta,&tm);
#endif
        char buf[40];
        if(!std::strftime(buf,sizeof(buf),"%Y-%m-%d %H:%M:%S",&tm))return "unknown";
        return buf;
    }
    void write_snapshot(const std::vector<ParticleTelemetry>& copy,const std::string& ts) {
        const auto path=(fs::path(root_)/"particle_status.csv").string();
        const auto temp=path+".tmp";
        {
            std::ofstream f(temp,std::ios::out|std::ios::trunc);
            if(!f)throw std::runtime_error("cannot write live particle status: "+temp);
            f<<"id,status,sim_day,altitude_km,max_eccentricity_so_far,energy_change_rel_signed,"
               "jacobi_change_rel_signed,max_lunar_jacobi_variation_rel,"
               "numeric_invariant_drift_rel,last_accepted_local_error_ratio,accepted_steps,impact_lat_deg,impact_lon_deg,"
               "terminal_speed_mps,error,updated_at\n";
            f<<std::scientific<<std::setprecision(15);
            for(size_t i=0;i<copy.size();++i) {
                const auto& x=copy[i];
                f<<i<<','<<label(x.status)<<','<<x.t_days<<','<<x.altitude_km<<','
                 <<x.eccentricity<<','<<x.energy_change_rel<<','<<x.jacobi_change_rel<<','
                 <<x.max_lunar_jacobi_variation_rel<<','<<x.numeric_invariant_drift_rel<<','
                 <<x.local_error_ratio<<','<<x.n_steps<<','<<x.impact_lat_deg<<','<<x.impact_lon_deg<<','
                 <<x.terminal_speed_mps<<',';
                // CSV-safe quoted error text (normally empty).
                f<<'"';for(char c:x.error){if(c=='"')f<<'"';if(c=='\n'||c=='\r')c=' ';f<<c;}
                f<<"\","<<ts<<"\n";
            }
            f.flush();if(!f)throw std::runtime_error("failed to flush live particle status");
        }
        fs::rename(temp,path);
    }
    void report(bool final,bool force=false) {
        std::vector<ParticleTelemetry> copy;
        int session_new=0, observed_done=0;
        const auto clock_now=std::chrono::steady_clock::now();
        {
            std::lock_guard<std::mutex> lock(mutex_);
            for(const auto& x:particles_) {
                if(x.status==ParticleTelemetry::impacted || x.status==ParticleTelemetry::survived)++observed_done;
            }
            if(!force && !should_write_progress(observed_done,last_logged_completed_,
                  std::chrono::duration<double>(clock_now-last_log_time_).count(),
                  cfg_.progress_heartbeat_seconds,final)) return;
            copy=particles_;
            session_new=newly_completed_session_;
        }
        int pending=0,running=0,impacted=0,survived=0,failed=0;
        double work=0, sum_day=0,min_day=INFINITY,max_day=0;
        double max_numeric=0,max_physical=0;
        for(const auto& x:copy){
            if(x.status==ParticleTelemetry::pending)++pending;
            if(x.status==ParticleTelemetry::running)++running;
            if(x.status==ParticleTelemetry::impacted)++impacted;
            if(x.status==ParticleTelemetry::survived)++survived;
            if(x.status==ParticleTelemetry::failed)++failed;
            if(x.status!=ParticleTelemetry::pending){
                sum_day+=x.t_days;min_day=std::min(min_day,x.t_days);max_day=std::max(max_day,x.t_days);
            }
            work+=(x.status==ParticleTelemetry::impacted||x.status==ParticleTelemetry::survived||x.status==ParticleTelemetry::failed)
                   ? 1.0 : std::min(1.0,x.t_days/cfg_.days);
            if(std::isfinite(x.numeric_invariant_drift_rel))
                max_numeric=std::max(max_numeric,x.numeric_invariant_drift_rel);
            if(std::isfinite(x.max_lunar_jacobi_variation_rel))
                max_physical=std::max(max_physical,x.max_lunar_jacobi_variation_rel);
        }
        const auto now=iso_time();
        const double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-started_).count();
        const int done=impacted+survived; // Failed is not a completed orbit.
        const double fraction=100.0*work/copy.size();
        timeline_<<now<<','<<wall<<','<<pending<<','<<running<<','<<done<<','<<impacted<<','
                 <<survived<<','<<failed<<','<<fraction<<','<<sum_day<<','
                 <<(std::isfinite(min_day)?min_day:NAN)<<','<<max_day<<','
                 <<(perturbed_?NAN:max_numeric)<<','<<max_physical<<'\n';
        timeline_.flush();
        write_snapshot(copy,now);
        std::string eta="unknown";
        if(failed>0) eta="unknown";
        else if(done==cfg_.n) eta=now;
        else if(session_new>0 && wall>0)
            eta=eta_timestamp(wall*double(cfg_.n-done)/session_new);
        simple_log_<<now<<" | "<<done<<'/'<<cfg_.n<<" completed | "
                   <<std::fixed<<std::setprecision(2)<<(100.0*done/cfg_.n)
                   <<"% | ETA: "<<eta<<'\n';
        simple_log_.flush();
        if(!simple_log_) throw std::runtime_error("failed to write simple progress log");
        last_log_time_=clock_now;
        last_logged_completed_=done;
    }
public:
    ProgressReporter(const Config& c,bool third_bodies):cfg_(c),perturbed_(third_bodies),
        root_(c.progress_dir.empty()?(fs::path(c.log).has_parent_path()?fs::path(c.log).parent_path().string():"."):c.progress_dir),
        particles_(static_cast<size_t>(c.n)) {
        fs::create_directories(root_);
        const auto events_path=(fs::path(root_)/"particle_events.csv").string();
        const auto progress_path=(fs::path(root_)/"progress.csv").string();
        const bool previous_events=c.resume&&fs::exists(events_path)&&fs::file_size(events_path)>0;
        const bool previous_progress=c.resume&&fs::exists(progress_path)&&fs::file_size(progress_path)>0;
        last_log_time_=started_;
        simple_log_.open(c.log,c.resume?std::ios::app:std::ios::trunc);
        events_.open(events_path,c.resume?std::ios::app:std::ios::trunc);
        timeline_.open(progress_path,c.resume?std::ios::app:std::ios::trunc);
        if(!events_||!timeline_||!simple_log_)throw std::runtime_error("cannot create progress output");
        if(!previous_events)events_<<"wall_timestamp,id,outcome,sim_day,steps,energy_change_rel_signed,jacobi_change_rel_signed,max_lunar_jacobi_variation_rel,numeric_invariant_drift_rel,last_accepted_local_error_ratio,impact_lat_deg,impact_lon_deg,terminal_speed_mps,error\n";
        if(!previous_progress)timeline_<<"wall_timestamp,elapsed_wall_s,pending,running,completed,impacted,survived,failed,work_fraction_percent,integrated_particle_days,min_started_sim_day,max_started_sim_day,max_numerical_invariant_drift_rel,max_lunar_jacobi_variation_rel\n";
        events_.flush();timeline_.flush();
    }
    void update(int i,const State& state,const Result& result,double alt,
                double E_now,double E_initial,double J_now,double local_error,bool terminal,bool restored=false) {
        ParticleTelemetry x;
        x.status=terminal?(result.impacted?ParticleTelemetry::impacted:ParticleTelemetry::survived):ParticleTelemetry::running;
        x.t_days=state.t/86400.0;
        x.altitude_km=alt/1000.0;
        x.n_steps=result.n_steps;
        x.local_error_ratio=local_error;
        x.energy_change_rel=(E_now-E_initial)/std::max(1.0,std::fabs(E_initial));
        x.jacobi_change_rel=(J_now-result.inv0)/std::max(1.0,std::fabs(result.inv0));
        x.max_lunar_jacobi_variation_rel=result.lunar_jacobi_variation_rel;
        x.numeric_invariant_drift_rel=perturbed_?NAN:result.invariant_drift_rel;
        if(result.impacted){x.impact_lat_deg=result.impact_lat;x.impact_lon_deg=result.impact_lon;x.terminal_speed_mps=result.terminal_speed;}
        // Eccentricity of the instantaneous osculating monopole orbit:
        // the caller supplies the already-computed maximum as a conservative
        // monitor value, not an instantaneous eccentricity.
        x.eccentricity=result.ecc_max;
        std::lock_guard<std::mutex> lock(mutex_);
        auto& old=particles_.at(static_cast<size_t>(i));
        const bool was_terminal=(old.status==ParticleTelemetry::impacted || old.status==ParticleTelemetry::survived);
        old=std::move(x);
        if(terminal && !was_terminal && !restored){
            ++newly_completed_session_;
            events_<<iso_time()<<','<<i<<','<<label(old.status)<<','<<std::setprecision(15)
                   <<old.t_days<<','<<old.n_steps<<','<<old.energy_change_rel<<','
                   <<old.jacobi_change_rel<<','<<old.max_lunar_jacobi_variation_rel<<','
                   <<old.numeric_invariant_drift_rel<<','<<old.local_error_ratio<<','<<old.impact_lat_deg<<','
                   <<old.impact_lon_deg<<','<<old.terminal_speed_mps<<",\n";
            events_.flush();
        }
    }
    void error(int i,const std::string& message){
        std::lock_guard<std::mutex> lock(mutex_);
        auto& x=particles_.at(static_cast<size_t>(i));
        x.status=ParticleTelemetry::failed;x.error=message;
        events_<<iso_time()<<','<<i<<",failed,"<<std::setprecision(15)<<x.t_days<<','<<x.n_steps
               <<','<<x.energy_change_rel<<','<<x.jacobi_change_rel<<','
               <<x.max_lunar_jacobi_variation_rel<<','<<x.numeric_invariant_drift_rel<<','<<x.local_error_ratio<<','
               <<x.impact_lat_deg<<','<<x.impact_lon_deg<<','<<x.terminal_speed_mps<<",\"";
        for(char c:message){if(c=='"')events_<<'"';if(c=='\n'||c=='\r')c=' ';events_<<c;}
        events_<<"\"\n";events_.flush();
    }
    void start(){
        if(!cfg_.resume)report(false,true);
        monitor_=std::thread([this]{
            std::unique_lock<std::mutex> lock(mutex_);
            while(!stop_){
                if(wake_.wait_for(lock,std::chrono::duration<double>(cfg_.progress_every_seconds),[&]{return stop_;}))break;
                lock.unlock();
                try{report(false);}catch(const std::exception& e){
                    std::fprintf(stderr,"WARNING: progress log failed: %s\n",e.what());
                }
                lock.lock();
            }
        });
    }
    void finish(){
        {std::lock_guard<std::mutex> lock(mutex_);stop_=true;}
        wake_.notify_all();
        if(monitor_.joinable())monitor_.join();
        report(true);
    }
    ~ProgressReporter(){
        {std::lock_guard<std::mutex> lock(mutex_);stop_=true;}
        wake_.notify_all();if(monitor_.joinable())monitor_.join();
    }
};
