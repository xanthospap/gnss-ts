#ifndef __NGPT_VINCENTY__HPP
#define __NGPT_VINCENTY_HPP

#include <cmath>
#include <stdexcept>
#include "ellipsoid.hpp"

namespace ngpt {

/// Given the (ellipsoidal) coordinates of two points (1 and 2) in radians,
/// calculate the forward azimouths from 1->2 (i.e. a12) and fro 2->1 (i.e. a21)
/// and the ellipsoidal distance between the two points. The inverse Vincenty's
/// formula is used (https://en.wikipedia.org/wiki/Vincenty's_formulae)
///
template<ellipsoid E = ellipsoid::grs80>
    double inverse_vincenty(double f1, double l1, double f2, double l2,
        double& a12, double& a21, double convergence_limit = 1e-8)
{
    const int MAX_ITERATIONS = 100;
    int iteration = 0;

    double a = ellipsoid_traits<E>::a;
    double f = ellipsoid_traits<E>::f;
    double b = semi_minor<E>();

    double  U1 { std::atan((1-f)*std::tan(f1)) },
            U2 { std::atan((1-f)*std::tan(f2)) },
            L  { l2 -l1 },
            l  {L}, new_l {l};
    double sins, coss, s, sina, cosa2, sinl, cosl, C, sinU1, cosU1, sinU2,
        cosU2, cos2sm;
    sinU1 = std::sin(U1);
    sinU2 = std::sin(U2);
    cosU1 = std::cos(U1);
    cosU2 = std::cos(U2);

    do {
        if (++iteration > MAX_ITERATIONS) {
            throw std::out_of_range("Inverse Vincenty cannot converge after 1000 iterations!");
        }
        l     = new_l;
        sinl  = std::sin(l);
        cosl  = std::cos(l);
        sins  = std::sqrt( (cosU2*sinl)*(cosU2*sinl) +
                (cosU1*sinU2-sinU1*cosU2*cosl)*(cosU1*sinU2-sinU1*cosU2*cosl) );
        coss  = sinU1*sinU2 + cosU1*cosU2*cosl;
        s     = atan2(sins, coss);
        sina  = cosU1*cosU2*sinl / sins;
        cosa2 = 1 - sina*sina;
        cos2sm= coss - (2*sinU1*sinU2)/cosa2;
        C     = (f/16)*cosa2*(4+f*(4-3*cosa2));
        new_l = L + (1-C)*f*sina*(s+C*sins*(cos2sm+C*coss*(-1+s*cos2sm*cos2sm)));
    } while ( std::abs(l-new_l) < convergence_limit );

    l = new_l;
    sinl  = std::sin(l);
    cosl  = std::cos(l);
    double u2 { cosa2*(a*a-b*b)/(b*b) };
    double A  { 1 + (u2*(4096+u2*(-768+u2*(320-175*u2))))/16384 };
    double B  { u2*(256+u2*(-128+u2*(74-47*u2))) / 1024 };
    double Ds { B*sins*(cos2sm + 0.25*B*(coss*(-1+2*cos2sm*cos2sm)-
            (1/6.0)*B*cos2sm*(-3+4*sins*sins)*(-3+4*cos2sm*cos2sm))) };
    double distance { b*A*(s-Ds) };
    
    a12 = atan2(cosU2*sinl, cosU1*sinU2-sinU1*cosU2*cosl);
    a21 = atan2(cosU1*sinl, -sinU1*cosU2+cosU1*sinU2*cosl);
    return distance;
}

} // end namespace ngpt
#endif
