#ifndef __NGPT_VINCENTY__HPP
#define __NGPT_VINCENTY_HPP

#include <cmath>
#include <stdexcept>
#include "ellipsoid.hpp"
#ifdef DEBUG
#include <iostream>
#include <iomanip>
#endif

namespace ngpt {

/// Compute the haversine formula.
double haversine(double angle) noexcept
{
    double sinHalfTheta { std::sin(angle/2.0) };
    return sinHalfTheta * sinHalfTheta;
}

/// Given the (ellipsoidal) coordinates of two points (1 and 2) in radians,
/// compute the great circle distance between them, using the haversine formula
/// see https://en.wikipedia.org/wiki/Haversine_formula. The Earth's radius is
/// computed as an average of the Normal radious of curvature at the two points.
template<ellipsoid E = ellipsoid::wgs84>
    double haversine(double lat1, double lon1, double lat2, double lon2)
{
    // double r1 { N<E>(lat1) };
    // double r2 { N<E>(lat2) };
    double h  { haversine(lat2-lat1) + std::cos(lat1) * std::cos(lat2) * 
        haversine(lon2-lon1) };
    // double EarthRadius { r1/2 + r2/2 };
    double EarthRadius { (2*ellipsoid_traits<E>::a + semi_minor<E>())/3 };

    return 2 * EarthRadius * std::asin(std::sqrt(h));
}


/// Given the (ellipsoidal) coordinates of two points (1 and 2) in radians,
/// calculate the forward azimouths from 1->2 (i.e. a12) and fro 2->1 (i.e. a21)
/// and the ellipsoidal distance between the two points. The inverse Vincenty's
/// formula is used (https://en.wikipedia.org/wiki/Vincenty's_formulae)
///
template<ellipsoid E = ellipsoid::wgs84>
    double inverse_vincenty(double lat1, double lon1, double lat2, double lon2,
        double& a12, double& a21, double convergence_limit = 1e-10)
{
    const int MAX_ITERATIONS = 100;
    int iteration = 0;

    double a = ellipsoid_traits<E>::a;
    double f = ellipsoid_traits<E>::f;
    double b = semi_minor<E>();

    double U1     { std::atan((1-f)*std::tan(lat1)) };
    double U2     { std::atan((1-f)*std::tan(lat2)) };
    double L      { lon2 -lon1 };
    double sinU1  { std::sin(U1) };
    double sinU2  { std::sin(U2) };
    double cosU1  { std::cos(U1) };
    double cosU2  { std::cos(U2) };
    double lambda {L};
    double sinSigma, cosSigma, sigma, sinAlpha, cosSqAlpha, C, cos2SigmaM,
           lambdaP, sinLambda, cosLambda;

    do {
        if (++iteration > MAX_ITERATIONS) {
            throw std::out_of_range("Inverse Vincenty cannot converge after 1000 iterations!");
        }
        sinLambda  = std::sin(lambda);
        cosLambda  = std::cos(lambda);
        sinSigma   = std::sqrt( (cosU2*sinLambda) * (cosU2*sinLambda) +
                (cosU1*sinU2-sinU1*cosU2*cosLambda) * (cosU1*sinU2-sinU1*cosU2*cosLambda) );
        cosSigma   = sinU1*sinU2 + cosU1*cosU2*cosLambda;
        sigma      = atan2(sinSigma, cosSigma);
        sinAlpha   = cosU1*cosU2*sinLambda / sinSigma;
        cosSqAlpha = 1.0 - sinAlpha*sinAlpha;
        cos2SigmaM = cosSigma - 2.0*sinU1*sinU2/cosSqAlpha;
        C          = (f/16.0)*cosSqAlpha*(4.0+f*(4.0-3.0*cosSqAlpha));
        lambdaP    = lambda;
        lambda     = L + (1.0-C)*f*sinAlpha*
                    (sigma+C*sinSigma*(cos2SigmaM+C*cosSigma*
                    (-1.0+2.0*cos2SigmaM*cos2SigmaM)));
    } while ( std::abs(lambda-lambdaP) > convergence_limit );
#ifdef DEBUG
    std::cout<<"\tVincenty Inverse converged after "<<iteration<<" iterations\n";
#endif

    double uSq { cosSqAlpha*(a*a-b*b)/(b*b) };
    //double A  { 1 + (u2*(4096+u2*(-768+u2*(320-175*u2))))/16384 };
    //double B  { u2*(256+u2*(-128+u2*(74-47*u2))) / 1024 };
    double k1 { (std::sqrt(1.0+uSq)-1.0)/(std::sqrt(1.0+uSq)+1.0) };
    double A  { (1+0.25*k1*k1)/(1.0-k1) };
    double B  { k1*(1.0-(3.0/8.0)*k1*k1) };
    double deltaSigma { B*sinSigma*(cos2SigmaM+B/4.0*(cosSigma*
            (-1.0+2.0*cos2SigmaM*cos2SigmaM)-B/6.0*cos2SigmaM*
            (-3.0+4.0*sinSigma*sinSigma)*(-3.0+4.0*cos2SigmaM*cos2SigmaM))) };
    double distance { b*A*(sigma-deltaSigma) };
    
    a12 = std::atan2(cosU2*sinLambda, cosU1*sinU2-sinU1*cosU2*cosLambda);
    a21 = std::atan2(cosU1*sinLambda, -sinU1*cosU2+cosU1*sinU2*cosLambda);
    return distance;
}

} // end namespace ngpt
#endif
