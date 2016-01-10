// #include "ellipsoid.hpp"
#include "car2ell.hpp"
#include "ell2car.hpp"
#include "car2top.hpp"

#include <stdio.h>
#include <cmath>

struct Point {
    double x,y,z;
};

using namespace geodesy;

int main ()
{
    Point p1 {4595220.02261,2039434.13622,3912625.96044};
    Point p2,p3;

    car2ell<ELLIPSOID::GRS80>(p1.x,p1.y,p1.z,p2.x,p2.y,p2.z);
    ell2car<ELLIPSOID::GRS80>(p2.x,p2.y,p2.z,p3.x,p3.y,p3.z);

    printf ("To ellipsoidal and back:");
    printf ("\ndx=%8.5f dy=%8.5F dz=%8.5f\n",
            std::abs(p1.x-p3.x),std::abs(p1.y-p3.y),std::abs(p1.z-p3.z));

    car2top<ELLIPSOID::GRS80>(p1.x,p1.y,p1.z,p3.x,p3.y,p3.z,p2.x,p2.y,p2.z);
    printf ("Topocentric vector (same as above)");
    printf ("\ndn=%8.5f de=%8.5F du=%8.5f\n",p2.x,p2.y,p2.z);

    return 0;
}
