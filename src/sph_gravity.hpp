// sph_gravity.hpp - Lunar gravity from fully-normalized spherical harmonics.
//
// Fixes vs moon_orbit_v0.py:
//   * NO order truncation. Full inner sum m = 0..n (v0 used only m in {0,1,2,l}).
//   * NO radial-only gradient hack. Acceleration is the true gradient of the
//     SAME potential used for the energy audit -> force/energy are consistent by
//     construction, so the symplectic integrator conserves total energy.
//   * Fully-normalized (4pi / geodesy) associated Legendre functions via the
//     Holmes & Featherstone (2002) forward-column recursion, numerically stable
//     to high degree. GRGM1200B coefficients are 4pi-normalized -> match.
//
// Potential (positive convention, acceleration a = grad U):
//   U = (GM/r) * sum_{n=0}^{N} (R/r)^n * sum_{m=0}^{n}
//         Pbar_{nm}(sin phi) * [ Cbar_{nm} cos(m lam) + Sbar_{nm} sin(m lam) ]
// Specific potential ENERGY = -U ; total specific energy E = 1/2 v^2 - U.
//
// Thread-safe: no mutable state; all scratch is thread_local. One shared
// const SphGravity instance may be evaluated from many OpenMP threads.
#pragma once
#include "vec3.hpp"
#include <vector>
#include <cstdint>
#include <cmath>
#include <string>
#include <fstream>
#include <stdexcept>

class SphGravity {
public:
    int lmax = 0;
    double GM = 0.0;   // m^3 / s^2
    double R  = 0.0;   // reference radius, m

    void load(const std::string& path, int lmax_use = -1) {
        std::ifstream f(path, std::ios::binary);
        if (!f) throw std::runtime_error("cannot open gravity file: " + path);
        int32_t lf; f.read(reinterpret_cast<char*>(&lf), sizeof(lf));
        f.read(reinterpret_cast<char*>(&GM), sizeof(double));
        f.read(reinterpret_cast<char*>(&R),  sizeof(double));
        if (!f || lf < 0 || lf > 2500 || !std::isfinite(GM) || GM <= 0.0 || !std::isfinite(R) || R <= 0.0)
            throw std::runtime_error("invalid gravity file header: " + path);
        int Lf = lf;
        if (lmax_use > Lf) throw std::runtime_error("requested harmonic degree exceeds data file: " + path);
        lmax = (lmax_use < 0) ? Lf : lmax_use;
        alloc(lmax);
        for (int n = 0; n <= Lf; ++n)
            for (int m = 0; m <= n; ++m) {
                double c, s;
                f.read(reinterpret_cast<char*>(&c), sizeof(double));
                f.read(reinterpret_cast<char*>(&s), sizeof(double));
                if (n <= lmax) { C[idx(n, m)] = c; S[idx(n, m)] = s; }
            }
        if (!f) throw std::runtime_error("gravity file truncated: " + path);
        if (std::fabs(C[idx(0, 0)] - 1.0) > 1e-6)
            throw std::runtime_error("gravity coefficients require fully normalized C00 = 1");
        precompute();
    }

    // Fully-normalized potential U (positive convention). Thread-safe.
    double potential(const Vec3& p) const {
        Scratch& sc = scratch();
        const double r = p.norm();
        if (!(r > 0.0) || !std::isfinite(r)) throw std::runtime_error("invalid gravity evaluation position");
        const double sphi = p.z / r;                                  // sin(lat)
        const double cphi = std::sqrt(std::max(0.0, 1.0 - sphi * sphi));
        const double lam = std::atan2(p.y, p.x);
        const double Rr = R / r;

        // cos(m lam), sin(m lam) and (R/r)^n by recurrence
        double cl = std::cos(lam), sl = std::sin(lam);
        sc.cosm[0] = 1.0; sc.sinm[0] = 0.0; sc.Rrn[0] = 1.0;
        for (int m = 1; m <= lmax; ++m) {
            sc.cosm[m] = sc.cosm[m - 1] * cl - sc.sinm[m - 1] * sl;
            sc.sinm[m] = sc.sinm[m - 1] * cl + sc.cosm[m - 1] * sl;
            sc.Rrn[m] = sc.Rrn[m - 1] * Rr;
        }

        // Sectoral seeds Pbar_{m,m}
        sc.diag[0] = 1.0;
        // 4pi-normalized geodesy functions (csphase=+1): P11 = sqrt(3) cos(phi).
        // m=1 is exceptional; using the m>=2 recurrence for it introduces
        // a spurious 1/sqrt(2) factor in ALL tesseral/sectoral orders.
        if (lmax >= 1) sc.diag[1] = std::sqrt(3.0) * cphi;
        for (int m = 2; m <= lmax; ++m)
            sc.diag[m] = cphi * std::sqrt((2.0 * m + 1.0) / (2.0 * m)) * sc.diag[m - 1];

        double U = 0.0;
        for (int m = 0; m <= lmax; ++m) {
            double Pnm2 = 0.0, Pnm1 = 0.0;
            for (int n = m; n <= lmax; ++n) {
                double Pbar;
                if (n == m)          Pbar = sc.diag[m];
                else if (n == m + 1) Pbar = sphi * std::sqrt(2.0 * m + 3.0) * sc.diag[m];
                else                 Pbar = A_[idx(n, m)] * sphi * Pnm1 - B_[idx(n, m)] * Pnm2;
                Pnm2 = Pnm1; Pnm1 = Pbar;
                double harm = C[idx(n, m)] * sc.cosm[m] + S[idx(n, m)] * sc.sinm[m];
                U += sc.Rrn[n] * Pbar * harm;
            }
        }
        return (GM / r) * U;
    }

    // Acceleration a = grad U, ANALYTIC (spherical-coordinate gradient of the same
    // normalized potential -> consistent with the energy audit, no finite-diff
    // floor, and one field evaluation instead of six). Non-singular except within
    // ~1e-12 rad of a pole (cos(lat) floored). Thread-safe.
    Vec3 acceleration(const Vec3& p) const {
        Scratch& sc = scratch();
        const double x = p.x, y = p.y, z = p.z;
        const double r = std::sqrt(x * x + y * y + z * z);
        if (!(r > 0.0) || !std::isfinite(r)) throw std::runtime_error("invalid gravity evaluation position");
        const double rho = std::sqrt(x * x + y * y);
        const double sphi = z / r, cphi = rho / r;
        // Spherical-coordinate derivative has a removable pole singularity.
        // Near the poles use the Cartesian gradient of the same potential,
        // avoiding cancellation and incorrect m=1 contributions at rho=0.
        if (cphi < 1e-7) return acceleration_fd(p);
        const double cphi_s = cphi;
        const double clam = (rho > 0.0) ? x / rho : 1.0;
        const double slam = (rho > 0.0) ? y / rho : 0.0;
        const double Rr = R / r;

        sc.cosm[0] = 1.0; sc.sinm[0] = 0.0; sc.Rrn[0] = 1.0;
        for (int m = 1; m <= lmax; ++m) {
            sc.cosm[m] = sc.cosm[m - 1] * clam - sc.sinm[m - 1] * slam;
            sc.sinm[m] = sc.sinm[m - 1] * clam + sc.cosm[m - 1] * slam;
            sc.Rrn[m] = sc.Rrn[m - 1] * Rr;
        }
        sc.diag[0] = 1.0;
        // 4pi-normalized geodesy functions (csphase=+1): P11 = sqrt(3) cos(phi).
        // m=1 is exceptional; using the m>=2 recurrence for it introduces
        // a spurious 1/sqrt(2) factor in ALL tesseral/sectoral orders.
        if (lmax >= 1) sc.diag[1] = std::sqrt(3.0) * cphi;
        for (int m = 2; m <= lmax; ++m)
            sc.diag[m] = cphi * std::sqrt((2.0 * m + 1.0) / (2.0 * m)) * sc.diag[m - 1];

        double dUdr_s = 0.0, dUdphi_s = 0.0, dUdlam_s = 0.0;
        for (int m = 0; m <= lmax; ++m) {
            double Pnm2 = 0.0, Pnm1 = 0.0;
            for (int n = m; n <= lmax; ++n) {
                double Pbar;
                if (n == m)          Pbar = sc.diag[m];
                else if (n == m + 1) Pbar = sphi * std::sqrt(2.0 * m + 3.0) * sc.diag[m];
                else                 Pbar = A_[idx(n, m)] * sphi * Pnm1 - B_[idx(n, m)] * Pnm2;
                // dPbar/dphi = (e_nm * Pbar_{n-1,m} - n*sin(phi)*Pbar) / cos(phi)
                double dP = (E_[idx(n, m)] * Pnm1 - n * sphi * Pbar) / cphi_s;
                Pnm2 = Pnm1; Pnm1 = Pbar;
                double cc = C[idx(n, m)] * sc.cosm[m] + S[idx(n, m)] * sc.sinm[m];
                double sc2 = S[idx(n, m)] * sc.cosm[m] - C[idx(n, m)] * sc.sinm[m];
                double Rn = sc.Rrn[n];
                dUdr_s   += (n + 1) * Rn * Pbar * cc;
                dUdphi_s += Rn * dP * cc;
                dUdlam_s += Rn * Pbar * m * sc2;
            }
        }
        const double GMr = GM / r;
        const double ar   = -(GMr / r) * dUdr_s;         // dU/dr
        const double aphi = (GMr * dUdphi_s) / r;        // (1/r) dU/dphi
        const double alam = (GMr * dUdlam_s) / (r * cphi_s);  // (1/(r cos phi)) dU/dlam
        // rotate spherical -> Cartesian (rhat, phihat, lamhat)
        return {ar * cphi * clam - aphi * sphi * clam - alam * slam,
                ar * cphi * slam - aphi * sphi * slam + alam * clam,
                ar * sphi + aphi * cphi};
    }

    // Acceleration via central finite difference of the potential (reference /
    // self-test only). h ~ 1 m.
    Vec3 acceleration_fd(const Vec3& p) const {
        const double h = 1.0, inv = 1.0 / (2.0 * h);
        return {(potential({p.x + h, p.y, p.z}) - potential({p.x - h, p.y, p.z})) * inv,
                (potential({p.x, p.y + h, p.z}) - potential({p.x, p.y - h, p.z})) * inv,
                (potential({p.x, p.y, p.z + h}) - potential({p.x, p.y, p.z - h})) * inv};
    }

    // Specific potential energy (negative convention) for the energy audit.
    double potential_energy(const Vec3& p) const { return -potential(p); }

private:
    std::vector<double> C, S, A_, B_, E_;

    struct Scratch { std::vector<double> cosm, sinm, diag, Rrn; };
    // one Scratch per thread, sized to lmax on first use
    Scratch& scratch() const {
        thread_local Scratch sc;
        if ((int)sc.cosm.size() != lmax + 1) {
            sc.cosm.assign(lmax + 1, 0.0); sc.sinm.assign(lmax + 1, 0.0);
            sc.diag.assign(lmax + 1, 0.0); sc.Rrn.assign(lmax + 1, 0.0);
        }
        return sc;
    }

    static inline int tri(int n) { return n * (n + 1) / 2; }
    inline int idx(int n, int m) const { return tri(n) + m; }

    void alloc(int L) {
        int sz = tri(L + 1);
        C.assign(sz, 0.0); S.assign(sz, 0.0);
        A_.assign(sz, 0.0); B_.assign(sz, 0.0); E_.assign(sz, 0.0);
    }

    void precompute() {
        for (int m = 0; m <= lmax; ++m)
            for (int n = m + 2; n <= lmax; ++n) {
                double nm = double(n - m), np = double(n + m);
                A_[idx(n, m)] = std::sqrt((2.0 * n - 1.0) * (2.0 * n + 1.0) / (nm * np));
                B_[idx(n, m)] = std::sqrt((2.0 * n + 1.0) * (np - 1.0) * (nm - 1.0) /
                                          (nm * np * (2.0 * n - 3.0)));
            }
        // latitude-derivative coefficient e_nm = sqrt((2n+1)(n-m)(n+m)/(2n-1))
        for (int n = 1; n <= lmax; ++n)
            for (int m = 0; m <= n; ++m)
                E_[idx(n, m)] = std::sqrt((2.0 * n + 1.0) * double(n - m) * double(n + m) /
                                          (2.0 * n - 1.0));
    }
};
