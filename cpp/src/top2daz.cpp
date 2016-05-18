#include <cmath>
#include <stdexcept>
#include "geodesy.hpp"

/** \details  Compute the Distance, Azimouth and Zenith distance given a
 *            vector expressed in a local, topocentric system (i.e. given
 *            the north, east and up components of the vector).
 *
 *  \param[in]  north    Vector north component, meters.
 *  \param[in]  east     Vector east component, meters.
 *  \param[in]  up       Vector up component, meters.
 *  \param[out] distance The length of the vector, meters.
 *  \param[out] azimouth The azimouth of the vector, radians [0,2*pi).
 *  \param[out] zenith   The zenith distance, radians [0,pi)
 * 
 *  \throw               std::runtime_error if zero division encountered.
 *
 *  Reference: Physical Geodesy, p. 210
 */
void 
ngpt::top2daz(double north, double east, double up,
              double& distance, double& azimouth, double& zenith)
{

  // spatial distance of vector
  distance  = std::sqrt (north*north + east*east + up*up);

  // check if zero distance or north are zero
  if ( (!distance) || (!north) ) {
      throw std::runtime_error("geodesy::top2daz -> Zero Division !!");
  }

  // azimouth
  double a { std::atan2(east, north) };

  // normalize to range [0-2pi)
  azimouth  = std::fmod(a, ngpt::D2PI);
  while (azimouth < .0) azimouth += ngpt::D2PI;

  // zenith angle [0-pi)
  zenith = std::acos (up / distance);

  // finished
  return;
}
