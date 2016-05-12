#ifndef _TSTOOLS_HPP
#define _TSTOOLS_HPP

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <fstream>

#include "vartools.hpp"
#include "cppbasic.hpp"

/* approximately one day in fractional years */
#define DD2Y (1.0 / 365.25)

/* date transformations */
int iauCal2jd(const int& iy,const int& im,const int& id,double& djm);
int iauJd2cal(const double& dj,int& iy,int& im,int& id,double& fd);
double jd2fy(const double& jd);
void fy2cal(const double& fy,int& iy,int& im,int& id);
std::string m2m(const int& im);

/* 
** a specific point in time-series; mainly used to read in
** a ts file 
*/
struct tspoint
{
    double time;
    double f,sf;
    double l,sl;
    double u,su;
    /* flags:
     * 0 -> rapid
     * 1 -> final
     */
    short int flag;
};

/*
** flags to mark ts data points
*/
enum FLAG {
  OK=0,
  OUTLIER=1,
  SKIP=2,
  DISC=3,/* discontinuity */
  VELC=4,/* velocity change  */
  ERQT=5/* earthquacke */
};

/*
** quick weight computation
*/
inline double weight(const double& var0,const double& sigma)
{
  return var0 / (sigma*sigma);
}

/*
** class to hold a ts of a specific component
** time information is external
*/
class tscomponent {
  public:
    tscomponent(int,double* p=NULL);
    tscomponent();
    ~tscomponent();
    /* resize the ts; may cause data loss */
    void resize(int);
    /* get size skipping all but the OK-flagged data points */
    int getrealsize() const;
    /* mark a data point */
    inline void mark(const FLAG& f,const int& i) {*(flags+i)=f;}
    /* mark a series of data points */
    void mark(const FLAG&,const int&,const int&);
    /* compute weighted average and its std deviation */
    int average (double& mean,double& smean) const;
    /* get size skipping all but the OK-flagged data points */
    int cleansize() const;
    /* mark all data points as OK */
    void cleanoutliers();
    /* members */
    int     size;
    double* time;
    double* vals;
    double* sigmas;
    FLAG*   flags;
};

/* read in a .g.cts file */
int readingcts(std::ifstream&,std::vector<tspoint>&,bool userapid=true);

/*
** ts class; it holds 3 ts component instants
** and a time array
*/
class timeseries
{
  public:
    /* initialize from a .g.cts file; sypply also statin name */
    timeseries(const std::string& file,const std::string& name,bool userapid=true);
    timeseries(std::vector<tspoint>&,int);
    ~timeseries();
    /* return a specific date */
    inline double date(const int& i) const {return *(time+i);} 
    /* transform from Julian dates to fractional year; will also transform
    the discontinuity vector */
    void time2fy();
    /* read in discontinuities from a .dis file. Corresponding dates
    will be marked with DISC in the time array */
    void readindisc(const std::string&);
    /* apply discontinuities */
    void applydisc();
    /* clean sizes for all components (only OK-flaged) */
    void cleansize(int&,int&,int&) const;
    /* cleans flags OUTLIER and SKIP */
    void cleanoutliers(int c=0);
    /* full size; equal to the size of each component */
    inline int size()  const {return x.size;}
    /* number of discontinuities */
    inline int dsize() const {return discontinuities.size();}
    inline int vsize() const {return velchanges.size();}
    inline int esize() const {return earthquakes.size();}
    /* return index to the closest element (for time component) */
    int  locate (const double&) const;
    /* given a window in fractiona years and an index to a time ptr,
     * this function will return the indexes to index - window/2
     * and index + window/2. Flags are not checked (!!)
     * */
    int window(const int& index,const double& window,int& start,int& stop);
    /* compute average value and variance of average */
    void average(double*,int c=0) const;
    /* ellipsoidal to topocentric */
    void ell2top(const double& lat0,const double& lon0,const double& hgt0);
    void ell2top();
    /* mark outliers using 3-sigma test */
    int markoutliers(double* res,const double& sigma,const int& c);
    /* simple print to stdout or file */
    void print (std::string="",int c=0,bool gmtformat=false) const;
    /* smouth print with 3-sigma; only works for printing in file with gmt format */
    void smprint(std::string="",int c=0,bool gmtformat=false) const;
    /* make matrices for least squares */
    bool makemat(double* A,double* l,const std::vector<double>& freqs,
      const double& sigma0,const int& rows,const int& cols,double& t0,
      int c=0,bool use_weights=true);
    /* solve least squares via QR */
    int lsqrsolve(double* A,double* b,double* x,const int& rows,
      const int& cols,const std::vector<double>& freqs,double& var,
      const double& sigma0,double* sigmas,int c=0,bool m=1);
    /* fit model to ts */
    void modelfitc(double* estimates,const int& c,const int& model_pars,
      const std::vector<double>& freqs,const double& s0,
      bool iterative=true,std::string rf="",bool usegmt=false);
    /* report fit results */
    void printfit(const int& parameters,double* estimates,double* sigmas,
      const double& apriori_sigma,const double& var0,
      const std::vector<double>& freqs,const int& outliers,
      const double& percent,const int& c,const int& iterations);
    /* make model line */
    void makemodel(const double& tl,const double& tu,const double& dt,
      double* estimates,const int& pars,const std::vector<double>& freqs,
      const double& t0,int c,std::string& modelfile,bool usegmt=false);
/* private members */
  private:
    double*             time;
    tscomponent         x,y,z;
    std::vector<double> discontinuities;
    std::vector<double> earthquakes;
    std::vector<double> velchanges;
    std::string         station;
    short int*          soltype;
};

#endif
