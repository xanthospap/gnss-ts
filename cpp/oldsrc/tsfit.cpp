#include "tstools.hpp"
#include "vartools.hpp"
#include <stdexcept>

using std::string;
using std::vector;

void help();

int main(int argc,char *argv[])
{
  int i,parameters;
  bool gmtdates=true;

  /* check command line arguments */
  if (argc<=1) {
    help();
    return 2;
  }

  /* assign name of ts file */
  std::string tsfile = argv[1];
  std::string station = "unknown";

  /* assign station name */
  station=argv[2];

  /* assign vector of frequency multipliers */
  vector<double> frequencies;
  for (i=3;i<argc;i++) {
    double f = ::atof(argv[i]);
    frequencies.push_back(f);
    printf("\nfrequency multiplier %1i is %5.2f",i,f);
  }

  /* create a timeseries instant */
  timeseries ts (tsfile,station);

  /* make sure it has enough size */
  if (ts.size()<10) {
    printf("\n### Too few epochs in time-series file! (%2i)",ts.size());
    printf("\n### Exiting with status 1");
    return 1;
  }

  /* read in discontinuities (if any) */
  std::string disfile = tsfile;
  try {
  	std::string cts = "g.cts"; 
  	disfile.replace(disfile.find(cts),cts.length(),"dis");
  	ts.readindisc(disfile);
  	printf("\nIntroduced %2i discontinuities",ts.dsize());
  } catch (std::out_of_range& e) {
	  std::cerr<<"\n*** Could not find discontinuity file!";
	  std::cerr<<"\n*** No discontinuity introduced!";
  }

  /* transform to topocentric */
  ts.ell2top();

  /* transform dates to decimal years */
  ts.time2fy();

  /* check if time-series is less than one year; if so
   * delete frequencies */
  if (ts.date(ts.size()-1)-ts.date(0)<1.0) {
	printf("\n$ Time-Series spans less than a year (%01.2f)",ts.date(ts.size()-1)-ts.date(0));
	printf("\n$ Frequency multipliers removed!");
	frequencies.clear();
  }

  /* compute number of parameters to be estimated */
  parameters = 2 + ts.dsize() + 2*frequencies.size() + ts.vsize();

  /* array of estimated parameters */
  double X[parameters];

  /* perform LS for each component */
  printf("\nCOMPONENT North\n================================");
  ts.modelfitc(X,1,parameters,frequencies,0.001,true,station+".north.res.dat",gmtdates);
  ts.print(station+".north.raw.dat",1,gmtdates);
  printf("\nCOMPONENT East\n================================");
  ts.modelfitc(X,2,parameters,frequencies,0.001,true,station+".east.res.dat",gmtdates);
  ts.print(station+".east.raw.dat",2,gmtdates);
  printf("\nCOMPONENT Up\n================================");
  ts.modelfitc(X,3,parameters,frequencies,0.001,true,station+".up.res.dat",gmtdates);
  ts.print(station+".up.raw.dat",3,gmtdates);
  

  /* DEBUGING */
  /*
  double window = 6.0 * DD2Y;
  int wsize,start,stop;
  for (int ii=0;ii<ts.size();ii++) {
	wsize = ts.window(ii,window,start,stop);
	printf("\n position %5i window size: %2i [",ii,wsize);
	//printf("%10.5f < %10.5f < %10.5f]",ts.date(start),ts.date(ii),ts.date(stop));
  }
  */
  /* done */
  printf("\n");
  return 0;
}

void help()
{
  printf ("\nHELP MESSAGE MPLA MPLA MLPA\n");
}
