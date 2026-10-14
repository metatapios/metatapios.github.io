// Third-body differential (direct-minus-indirect) gravity in lunar-centered
// inertial coordinates. The secular Kozai-Lidov dynamics, if present, arise
// from these full Newtonian forces; no ad-hoc extra "KL force" is applied.
#pragma once
#include "vec3.hpp"
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace thirdbody {
constexpr double GM_EARTH = 3.98600436e14;      // m^3/s^2
constexpr double GM_SUN = 1.32712440018e20;    // m^3/s^2
constexpr double GM_JUPITER = 1.266865319e17;  // m^3/s^2
constexpr double AU = 149597870700.0;          // IAU exact astronomical unit (m)
constexpr double DAY = 86400.0;
constexpr double PI = 3.14159265358979323846;

enum Body { EARTH=0, SUN=1, JUPITER=2 };
using Positions = std::array<Vec3,3>;   // Moon->body, in the simulation's inertial axes

inline std::string trim(std::string s) {
    const auto a=s.find_first_not_of(" \t\r\n");
    if(a==std::string::npos) return "";
    const auto b=s.find_last_not_of(" \t\r\n");
    return s.substr(a,b-a+1);
}
inline std::array<bool,3> select(const std::string& spec) {
    std::array<bool,3> result{false,false,false};
    if(spec=="none" || spec.empty()) return result;
    if(spec=="all") return {true,true,true};
    std::istringstream ss(spec);
    std::string token;
    while(std::getline(ss,token,',')) {
        token=trim(token);
        std::transform(token.begin(),token.end(),token.begin(),
                       [](unsigned char c){return static_cast<char>(std::tolower(c));});
        if(token=="earth")result[EARTH]=true;
        else if(token=="sun")result[SUN]=true;
        else if(token=="jupiter")result[JUPITER]=true;
        else throw std::invalid_argument("invalid third body: '"+token+"' (earth,sun,jupiter,all,none)");
    }
    if(!result[0] && !result[1] && !result[2])
        throw std::invalid_argument("empty --third-bodies list");
    return result;
}

// Test/demo ONLY: Keplerian Earth-Moon, Earth-Sun and Jupiter-Sun orbits with
// arbitrary initial phases. No date, planetary precession, lunar libration,
// or fully consistent many-body ephemeris. Never use for final science tables.
inline Vec3 kepler_orbit(double a, double e, double period_s,
                         double inclination, double mean_anomaly0, double t) {
    const double M=std::remainder(mean_anomaly0+2*PI*t/period_s,2*PI);
    double E=M;
    for(int k=0;k<16;++k) {
        const double correction=(E-e*std::sin(E)-M)/(1-e*std::cos(E));
        E-=correction;
        if(std::abs(correction)<2e-15) break;
    }
    const double x=a*(std::cos(E)-e);
    const double y=a*std::sqrt(1-e*e)*std::sin(E);
    return {x, y*std::cos(inclination), y*std::sin(inclination)};
}
inline Positions analytic_positions(double t) {
    const double d2r=PI/180.0;
    const Vec3 earth=kepler_orbit(384400e3,0.0549,27.321661*DAY,6.68*d2r,0,t);
    const Vec3 earth_to_sun=kepler_orbit(AU,0.0167,365.25636*DAY,
                                          1.54*d2r,1.3,t);
    const Vec3 sun_to_jupiter=kepler_orbit(5.2044*AU,0.0489,11.862615*365.25*DAY,
                                           2.84*d2r,2.2,t);
    return {earth, earth+earth_to_sun, earth+earth_to_sun+sun_to_jupiter};
}

// Optional dense ephemeris sampled from NASA NAIF/SPICE. File convention:
// time_s,earth_x_m,earth_y_m,earth_z_m,sun_x_m,sun_y_m,sun_z_m,
// jupiter_x_m,jupiter_y_m,jupiter_z_m
// All positions must be Moon-centered and expressed in the inertial axes of
// the satellite propagator. Time is measured from the run start epoch.
// Four-point (third-order) Lagrange interpolation for nonuniform samples.
class PositionTable {
    struct Row {double t; Positions pos;};
    std::vector<Row> rows_;
public:
    void load(const std::string& filename) {
        std::ifstream f(filename);
        if(!f) throw std::runtime_error("cannot open third-body ephemeris: "+filename);
        rows_.clear();
        std::string line;
        int line_no=0;
        while(std::getline(f,line)) {
            ++line_no;
            const auto comment=line.find('#');
            if(comment!=std::string::npos) line.resize(comment);
            line=trim(line);
            if(line.empty())continue;
            if(line.rfind("time_s,",0)==0)continue;
            std::array<double,10> values{};
            std::istringstream stream(line);
            std::string token;
            int n=0;
            while(std::getline(stream,token,',')) {
                if(n>=10)throw std::runtime_error("too many ephemeris columns on line "+std::to_string(line_no));
                token=trim(token);
                try {size_t used=0;values[n]=std::stod(token,&used);
                     if(used!=token.size() || !std::isfinite(values[n]))
                         throw std::invalid_argument("bad coordinate");}
                catch(...) {throw std::runtime_error("invalid ephemeris value on line "+std::to_string(line_no));}
                ++n;
            }
            if(n!=10)throw std::runtime_error("expected 10 ephemeris columns on line "+std::to_string(line_no));
            if(!rows_.empty() && !(values[0]>rows_.back().t))
                throw std::runtime_error("ephemeris timestamps must be strictly increasing");
            Row row{}; row.t=values[0];
            for(int j=0;j<3;++j) {
                row.pos[j]={values[1+3*j],values[2+3*j],values[3+3*j]};
                if(row.pos[j].norm()<1e6)
                    throw std::runtime_error("third-body position too small (incorrect units?)");
            }
            rows_.push_back(row);
        }
        if(rows_.size()<4) throw std::runtime_error("ephemeris needs at least four time samples");
    }
    double start()const{return rows_.front().t;}
    double end()const{return rows_.back().t;}
    size_t size()const{return rows_.size();}
    Positions at(double t)const {
        if(!std::isfinite(t) || t<start()-1e-8 || t>end()+1e-8)
            throw std::runtime_error("third-body ephemeris sampled outside time coverage at t="+std::to_string(t)+" s");
        t=std::clamp(t,start(),end());
        auto hi=std::lower_bound(rows_.begin(),rows_.end(),t,
                                 [](const Row& r,double time){return r.t<time;});
        int index=int(hi-rows_.begin());
        int i0=std::clamp(index-2,0,int(rows_.size())-4);
        Positions result{};
        for(int j=0;j<4;++j) {
            const Row& node=rows_[i0+j];
            double w=1.0;
            for(int k=0;k<4;++k)if(j!=k) {
                const Row& other=rows_[i0+k];
                w*=(t-other.t)/(node.t-other.t);
            }
            for(int b=0;b<3;++b)result[b]+=node.pos[b]*w;
        }
        return result;
    }
};

class Perturbations {
    std::array<bool,3> chosen_{};
    std::string mode_="none";
    PositionTable table_;
public:
    Perturbations()=default;
    Perturbations(const std::string& bodies,const std::string& mode,
                  const std::string& ephemeris_file):chosen_(select(bodies)),mode_(mode) {
        if(!enabled()) {
            if(mode!="none")throw std::invalid_argument("ephemeris specified without --third-bodies");
            return;
        }
        if(mode=="file") {
            if(ephemeris_file.empty())throw std::invalid_argument("--ephem-file is required for file mode");
            table_.load(ephemeris_file);
        } else if(mode!="analytic") {
            throw std::invalid_argument("select --ephem-file FILE (SPICE) or --ephem-mode analytic (toy)");
        }
    }
    bool enabled()const{return chosen_[0]||chosen_[1]||chosen_[2];}
    const std::array<bool,3>& enabled_bodies()const{return chosen_;}
    const std::string& mode()const{return mode_;}
    double table_start()const{return table_.start();}
    double table_end()const{return table_.end();}
    Positions positions(double t)const {
        if(!enabled())return {};
        return mode_=="file"?table_.at(t):analytic_positions(t);
    }
    Vec3 acceleration(const Vec3& q,double t)const {
        if(!enabled())return {};
        const Positions p=positions(t);
        constexpr double mu[3]={GM_EARTH,GM_SUN,GM_JUPITER};
        Vec3 result{};
        for(int b=0;b<3;++b)if(chosen_[b]) {
            const Vec3 r=p[b],d=r-q;
            const double r2=r.norm2(), d2=d.norm2();
            if(!(r2>1e12 && d2>1e12))throw std::runtime_error("orbit too close to a perturbing body's center");
            result += (d*(1.0/(d2*std::sqrt(d2)))-
                       r*(1.0/(r2*std::sqrt(r2))))*mu[b];
        }
        return result;
    }
    // Instantaneous third-body orbital plane for the classical quadrupole
    // Kozai diagnostic. This is NOT conserved with lunar harmonics, multiple
    // perturbers, or a precessing outer orbit.
    Vec3 reference_normal(double t)const {
        int b=EARTH;
        while(b<3 && !chosen_[b])++b;
        if(b==3)return {0,0,1};
        constexpr double delta=30.;
        const Vec3 R=positions(t)[b];
        const Vec3 v=(positions(t+delta)[b]-positions(t-delta)[b])*(0.5/delta);
        const Vec3 h=R.cross(v);
        if(h.norm()<1e-20)return {NAN,NAN,NAN};
        return h*(1.0/h.norm());
    }
};

} // namespace thirdbody
