#ifndef HYPSYS1D_LIMITERS_HPP
#define HYPSYS1D_LIMITERS_HPP

#include <cmath>
#include <algorithm>

inline double sign(double a) { return copysign(1.0, a); }

inline double minmod(double a, double b) {
    return 0.5 * (sign(a) + sign(b)) * std::min(std::abs(a), std::abs(b));
}

inline double maxmod(double a, double b) {
    return 0.5 * (sign(a) + sign(b)) * std::max(std::abs(a), std::abs(b));
}

inline double minmod(double a, double b, double c) {
    return minmod(a, minmod(b, c));
}

/// FVM slope limiters

//----------------LimitersBegin----------------  
struct MinMod {
    inline double operator()(double sL, double sR) const
    {
        return minmod(sL,sR);
    }
};

struct SuperBee {
    inline double operator()(double sL, double sR) const {
        double SL = minmod(2.0*sL, sR);
        double SR = minmod(sL, 2.0*sR);
        return maxmod(SL, SR);
    }
};

struct MonotonizedCentral {
    inline double operator()(double sL, double sR) const {
        return minmod(2.0*sR, minmod((sL+sR)*0.5, 2.0*sL));
    }
};
//----------------LimitersEnd----------------  


/// DG limiters
struct VanLeer {
    inline double operator()(double s, double sm, double sp) const
    {
        return 0.;
    }
};

struct Shu {
    explicit Shu (const double dx_) : dx (dx_) {}

    inline double operator()(double s, double sm, double sp) const {
        if (std::abs(s) < M*dx*dx) {
            return s;
        } else {
            return minmod(s, sm, sp);
        }
    }

  private:
    double dx;
    double M = 50;
};


#endif // HYPSYS1D_LIMITERS_HPP
