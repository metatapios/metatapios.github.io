// topography.hpp - Lunar surface radius R_surf(lat,lon) from a LOLA/topo grid.
//
// Used for:
//   * initial satellite placement: altitude measured above local terrain,
//   * impact detection: satellite falls when r < R_surf(lat,lon) (v0 used a
//     flat r - R_MOON <= 0, ignoring topography; craters/highlands shift impact
//     altitude by up to ~ +/-9 km).
//
// Grid convention MATCHES pyshtools .expand(grid='DH'):
//   binary: int32 nlat, int32 nlon, then nlat*nlon doubles = absolute radius [m]
//   row  j (0..nlat-1) -> latitude  lat = 90 - j*(180/nlat)   [deg], N to S,
//                         (note: excludes the south pole, like DH sampling)
//   col  i (0..nlon-1) -> longitude lon = i*(360/nlon)        [deg], wraps.
#pragma once
#include <vector>
#include <cstdint>
#include <cmath>
#include <string>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include "vec3.hpp"

class Topography {
public:
    int nlat = 0, nlon = 0;

    void load(const std::string& path) {
        std::ifstream f(path, std::ios::binary);
        if (!f) throw std::runtime_error("cannot open topo file: " + path);
        int32_t a, b;
        f.read(reinterpret_cast<char*>(&a), sizeof(a));
        f.read(reinterpret_cast<char*>(&b), sizeof(b));
        if (!f || a < 2 || b < 2 || a > 50000 || b > 100000 ||
            (uint64_t)a * (uint64_t)b > 200000000ULL)
            throw std::runtime_error("invalid topo dimensions: " + path);
        nlat = a; nlon = b;
        data_.resize((size_t)nlat * nlon);
        f.read(reinterpret_cast<char*>(data_.data()),
               (std::streamsize)(sizeof(double) * data_.size()));
        if (!f) throw std::runtime_error("topo file truncated: " + path);
        rmin_ = rmax_ = data_[0];
        for (double v : data_) {
            if (!std::isfinite(v) || v <= 0.0) throw std::runtime_error("invalid topo radius: " + path);
            rmin_ = std::min(rmin_, v); rmax_ = std::max(rmax_, v);
        }
        dlat_ = 180.0 / nlat;   // DH latitude spacing
        dlon_ = 360.0 / nlon;
        // Conservative mesh-gradient estimate, including wrapped longitude.
        // Near the DH south pole, longitude derivatives are coordinate-singular
        // under clamped extrapolation; the event search documents this limit.
        max_surface_slope_=0.0;
        const double delta_lat_rad=dlat_*M_PI/180.,delta_lon_rad=dlon_*M_PI/180.;
        for(int j=0;j<nlat;++j){
            double coslat=std::fabs(std::cos((90.0-j*dlat_)*M_PI/180.));
            double largest_lat=0,largest_lon=0;
            for(int i=0;i<nlon;++i){
                if(j+1<nlat)largest_lat=std::max(largest_lat,std::fabs(at(j+1,i)-at(j,i)));
                largest_lon=std::max(largest_lon,std::fabs(at(j,(i+1)%nlon)-at(j,i)));
            }
            double lat_slope=largest_lat/(delta_lat_rad*rmin_);
            double lon_slope=(largest_lon==0)?0.0:largest_lon/(delta_lon_rad*rmin_*std::max(coslat,0.5*delta_lat_rad));
            max_surface_slope_=std::max(max_surface_slope_,std::hypot(lat_slope,lon_slope));
        }
    }

    double rmin() const { return rmin_; } // lowest surface radius (safety floor)
    double rmax() const { return rmax_; }
    double max_surface_slope() const { return max_surface_slope_; }

    struct SurfaceSample {
        double radius, dlat, dlon;  // partials of radius wrt lat and lon [m/rad]
    };

    // The DH grid samples +90 deg but not -90 deg. Longitude wraps; latitude
    // is clamped at its sampled endpoints to avoid out-of-range reads.
    SurfaceSample surface_sample(double lat_rad, double lon_rad) const {
        double latd = lat_rad * 180.0 / M_PI;
        double lond = std::fmod(lon_rad * 180.0 / M_PI, 360.0);
        if (lond < 0) lond += 360.0;
        double fj = (90.0 - latd) / dlat_;
        double fi = lond / dlon_;
        int j0 = static_cast<int>(std::floor(fj));
        double tj = fj - j0;
        if (j0 < 0) { j0 = 0; tj = 0; }
        if (j0 >= nlat - 1) { j0 = nlat - 1; tj = 0; }
        int j1 = std::min(j0 + 1, nlat - 1);
        int i0 = static_cast<int>(std::floor(fi)) % nlon;
        if (i0 < 0) i0 += nlon;
        int i1 = (i0 + 1) % nlon;
        double ti = fi - std::floor(fi);
        double r00 = at(j0, i0), r01 = at(j0, i1);
        double r10 = at(j1, i0), r11 = at(j1, i1);
        double r0 = (1 - ti) * r00 + ti * r01;
        double r1 = (1 - ti) * r10 + ti * r11;
        double r = (1 - tj) * r0 + tj * r1;
        double dlat = j1 == j0 ? 0.0 : -(r1 - r0) / (dlat_ * M_PI / 180.0);
        double dlon = ((1 - tj) * (r01 - r00) + tj * (r11 - r10)) /
                      (dlon_ * M_PI / 180.0);
        return {r, dlat, dlon};
    }

    double surface_radius(double lat_rad, double lon_rad) const {
        return surface_sample(lat_rad, lon_rad).radius;
    }

    // Normal of radial surface r = R(lat,lon), not the spherical radial vector.
    // Locally piecewise differentiable; grid edges can have slope discontinuities.
    // At coordinate poles a regular topographic gradient is not defined by DH;
    // use a spherical normal and report this approximation in the documentation.
    Vec3 surface_normal(double lat, double lon) const {
        const auto g = surface_sample(lat, lon);
        double cl = std::cos(lat), sl = std::sin(lat);
        double co = std::cos(lon), so = std::sin(lon);
        Vec3 er{cl * co, cl * so, sl};
        if (std::fabs(cl) < 1e-7) return er;
        Vec3 elat{-sl * co, -sl * so, cl};
        Vec3 elon{-so, co, 0};
        Vec3 n = er - elat * (g.dlat / g.radius) - elon * (g.dlon / (g.radius * cl));
        return n * (1.0 / n.norm());
    }

    // Equatorial terrain sampling scale, not the actual map's spatial fidelity.
    double equatorial_cell_width_m(double radius) const {
        return radius * dlon_ * M_PI / 180.0;
    }

private:
    std::vector<double> data_;
    double rmin_ = 0, rmax_ = 0, dlat_ = 0, dlon_ = 0;
    double max_surface_slope_ = 0;
    double at(int j, int i) const { return data_[(size_t)j * nlon + i]; }
};
