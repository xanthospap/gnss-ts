#include <iostream>
#include <iomanip>
#include <tuple>
#include "vincenty.hpp"
#include "geoconst.hpp"

constexpr double deg2rad = ngpt::DPI / 180.0;

auto rad2hdeg(double rad)
{
    int sign {1};
    if (rad < 0) { rad *= -1.0; sign = -1;}
    double ddeg = rad * 180 / ngpt::DPI;
    int deg = static_cast<int>(ddeg);
    int min = static_cast<int>( (ddeg-static_cast<double>(deg))*60.0 );
    double sec = (ddeg - (deg+min/60.0)) * 3600.0;
    return std::make_tuple(deg*sign, min, sec);
}

void
print_hex_deg(std::tuple<int,int,double>& angle)
{
    std::cout << std::get<0>(angle) << " " << std::get<1>(angle) << " " << std::get<2>(angle);
    return;
}

int main(int argc, char* argv[])
{
    if (argc != 5 && argc != 13) {
        std::cerr<<"Usage: vincenty <lon1> <lat1> <lon2> <lat2> in decimal degres.\n";
        std::cerr<<"       Coordinates can also be given in format HH MM SS.SSSS\n";
        return 1;
    }
    
    double lon1, lat1, lon2, lat2;

    if (argc == 5) {
        lon1 = std::atof(argv[1]);
        lat1 = std::atof(argv[2]);
        lon2 = std::atof(argv[3]);
        lat2 = std::atof(argv[4]);
    } else {
        lon1 = std::abs(std::atof(argv[1]))  + (std::atof(argv[2])  + std::atof(argv[3])/60.0)/60.0;
        if (std::atof(argv[1]) < 0) lon1 *= -1;
        lat1 = std::abs(std::atof(argv[4]))  + (std::atof(argv[5])  + std::atof(argv[6])/60.0)/60.0;
        if (std::atof(argv[4]) < 0) lat1 *= -1;
        lon2 = std::abs(std::atof(argv[7]))  + (std::atof(argv[8])  + std::atof(argv[9])/60.0)/60.0;
        if (std::atof(argv[7]) < 0) lon2 *= -1;
        lat2 = std::abs(std::atof(argv[10])) + (std::atof(argv[11]) + std::atof(argv[12])/60.0)/60.0;
        if (std::atof(argv[10]) < 0) lat2 *= -1;
    }

    double a12, a21, distance;
    distance = ngpt::inverse_vincenty<ngpt::ellipsoid::wgs84>(lat1*deg2rad, lon1*deg2rad,
        lat2*deg2rad, lon2*deg2rad, a12, a21, 1e-12);

    auto a12d = rad2hdeg(a12);
    auto a21d = rad2hdeg(a21);

    std::cout << "Point 1      : (" << lon1 << ", " << lat1 << ")\n";
    std::cout << "Point 2      : (" << lon2 << ", " << lat2 << ")\n";
    std::cout << "Distance     : "  << std::fixed << std::setprecision(3) << distance << "\n";
    std::cout << "Azimouth 1->2: "  << a12 << " or "; print_hex_deg(a12d);
    std::cout << "\nAzimouth 2->1: "<< a21 << " or "; print_hex_deg(a21d);
    // std::cout << "\nPoint 1      : (" << lon1*deg2rad << ", " << lat1*deg2rad << ") rad.\n";
    // std::cout << "Point 2      : (" << lon2*deg2rad << ", " << lat2*deg2rad << ") rad.\n";
    
    distance = ngpt::haversine<ngpt::ellipsoid::wgs84>(lat1*deg2rad, lon1*deg2rad,
        lat2*deg2rad, lon2*deg2rad);
    std::cout << "\nDistance Haversine (wgs84): "  << std::fixed << std::setprecision(3) << distance << "\n";
    
    distance = ngpt::haversine(lat1*deg2rad, lon1*deg2rad, lat2*deg2rad, lon2*deg2rad);
    std::cout << "\nDistance Haversine (using Earth's radius): "  << std::fixed << std::setprecision(3) << distance << "\n";
    
    std::cout<<"\n";

    return 0;
}
