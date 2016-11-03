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
    std::cout<<"\nAnalysis output written to file: "<< (sname + ".prd");
    std::string ofile = sname + ".prd";
     
    ngpt::crdts<ngpt::milliseconds> ts = ngpt::readin<ngpt::milliseconds>(ctsf, sname);

    std::vector<ngpt::timeseries<ngpt::milliseconds, ngpt::pt_marker>*> components;
    components.push_back(&ts.x_component());
    components.push_back(&ts.y_component());
    components.push_back(&ts.z_component());

    std::vector<std::string> cmp_names = {std::string("North"), std::string("East"), std::string("Up")};

    std::ofstream fout {ofile.c_str()};

    auto   tdif = ts.last_epoch().delta_date( ts.first_epoch() );
    double ddif = tdif.days().as_underlying_type() / 365.25;
    double div  = 1e0/ddif;
    std::cout<<"\nDiv = " << div;

    auto it  = components.begin();
    auto cit = cmp_names.cbegin();
    for (; it != components.end(); ++it) {
        std::size_t N = (*it)->data_pts() - (*it)->skipped_pts();
        double          ofac{4}, hifac{div};
        int             nout = 0.5*ofac*hifac*N + 1;
        double          *px, *py, prob;
        int             jmax;
        double          days_in_year = 365.25e0;

        px = new double[nout];
        py = new double[nout];
        ngpt::lomb_scargle_period( **it, ofac, hifac, px, py, nout, nout, jmax, prob );
        fout << "\n#New Component : " << *cit;
        for (int i=0;i<nout;i++) {
            fout << "\n" << px[i] << " " << py[i] << " " << 1e0/px[i];
        }
        std::cout<<"\nComponent: "<<*cit;
        std::cout<<"\n---------------------------------------------------------";
        std::cout<<"\n\tDominant frequency in time-series: "<<px[jmax]<<" (at: "<<jmax<<")";
        std::cout<<"; this is a period of "<<1e0/px[jmax]<<" days";
        std::cout<<"\n\tMinimum frequency examined is: "<<px[0]
            <<", i.e. a period of "<<1e0/px[0]<<" days";
        std::cout<<"\n\tMaximum frequency examined is: "<<px[nout-1]
            <<", i.e. a period of "<<1e0/px[nout-1]<<" days\n";

        ++cit;
        delete[] px;
        delete[] py;
    }

    fout.close();

    return 0;
}
