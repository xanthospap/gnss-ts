#include <iostream>
#include <fstream>
#include "cts_read.hpp"

void help()
{
    std::cout<<"\nProgram lomb-scargle";
    std::cout<<"\nPurpose Read in a (coordinate) time-series file, and estimate";
    std::cout<<"\n        the Lomb-Scargle periodogram.";
    std::cout<<"\nUsage   lomb_scargle ts-file";
    return;
}

std::string
split_path(std::string s)
{
    auto pos = s.find_last_of('/');
    if (pos == std::string::npos ) {
        if ( s.size() < 4) return std::string("xxxx");
        return s.substr(0, 4);
    }
    return s.substr(pos+1, 4);
}

int
main(int argc, char* argv[])
{
    if (argc != 2) {
        help();
        return 1;
    }

    std::string ctsf  = std::string(argv[1]);
    std::string sname = split_path(ctsf);
    std::string ofile = sname + ".prd";
     
    ngpt::crdts<ngpt::milliseconds> ts = ngpt::readin<ngpt::milliseconds>(ctsf, sname);

    std::vector<ngpt::timeseries<ngpt::milliseconds, ngpt::pt_marker>*> components;
    components.push_back(&ts.x_component());
    components.push_back(&ts.y_component());
    components.push_back(&ts.z_component());

    std::ofstream fout {ofile.c_str()};

    auto it = components.begin();
    for (; it != components.end(); ++it) {
        std::size_t N = (*it)->data_pts() - (*it)->skipped_pts();
        double          ofac{4}, hifac{.9};
        int             nout = 0.5*ofac*hifac*N + 1;
        double          *px, *py, prob;
        int             jmax;
        double          days_in_year = 365.25e0;

        px = new double[nout];
        py = new double[nout];
        ngpt::lomb_scargle_period( **it, ofac, hifac, px, py, nout, nout, jmax, prob );
        // std::cout<<"\nDominant frequency in time-series: "<<py[jmax]<<" (at: "<<jmax<<")";
        // std::cout<<"\nAs yearly fraction: 1/"<<days_in_year*px[jmax]<<"\n";
        fout << "\nNew Component";
        for (int i=0;i<nout;i++) {
            fout << "\n" << px[i] << " " << py[i];
        }

        delete[] px;
        delete[] py;
    }

    fout.close();

    return 0;
}
