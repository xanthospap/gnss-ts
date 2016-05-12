#include "tstools.hpp"
#include "matlib.hpp"

/* INCLUDE CBLAS AND LAPACKE HEADERS */
 #include "/usr/include/lapacke.h"
//#include "/home/xanthos/include/lapacke.h"
#include "/usr/include/cblas.h"
//#include "/home/xanthos/include/cblas.h"

/*
** =============================================================================
** DATE TRANSFORMATIONS
** =============================================================================
*/
std::string m2m(const int& im)
{
  switch (im) {
    case (1):
      return "Jan";
    case (2):
      return "Feb";
    case (3):
      return "Mar";
    case (4):
      return "Apr";
    case (5):
      return "May";
    case (6):
      return "Jun";
    case (7):
      return "Jul";
    case (8):
      return "Aug";
    case (9):
      return "Sep";
    case (10):
      return "Oct";
    case (11):
      return "Nov";
    case (12):
      return "Dec";
    default:
      std::cerr<<"\n*** Invalid month!!!";
      return "NAN";
  }
}
void fy2cal(const double& fy,int& iy,int& im,int& id)
{
  double iyd,fdec,daysinyear;
  long leap,yd,guess,more;

  static long month_day[2][13] = {
    {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365},
    {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366}
  };

  /* compute day of year */
  fdec = ::modf(fy,&iyd);
  iy = (int)iyd;
  leap = (iy%4==0);
  daysinyear = 365.0e00 + (double)leap;
  yd = fdec * daysinyear;

  /* transform to decimal */
  guess = yd*0.032;
  more = ((yd-month_day[leap][guess+1])>0);
  im = guess + more + 1;
  id = yd - month_day[leap][guess+more];

  /* done */
  return;
}

double jd2fy(const double& jd)
{
  /* julian date to franctional year */
  int iy,im,id;
  long leap;
  double fd;
  double yd;

  static long month_day[2][12] = {
    {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334},
    {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335}
  };

  if (iauJd2cal(jd,iy,im,id,fd)) {
    std::cerr<<"\n*** Invalid date!";
    return -999.9;
  }

  leap = (iy%4==0);

  /* day of year */
  yd = month_day[leap][im-1] + id + 0.2;

  return (double)iy + yd / (365.0e00+leap);
}

void timeseries::time2fy()
{
	int SIZE=size(),i;

	for (i=0;i<SIZE;i++) time[i]=jd2fy(time[i]);

	for (i=0;i<dsize();i++) discontinuities[i]=jd2fy(discontinuities[i]);

	for (i=0;i<vsize();i++) velchanges[i]=jd2fy(velchanges[i]);

	for (i=0;i<esize();i++) earthquakes[i]=jd2fy(earthquakes[i]);

	return;
}

/*
** =============================================================================
** TSCOMPONENT CLASS
** =============================================================================
*/
tscomponent::tscomponent()
  :size(0),vals(NULL),sigmas(NULL),flags(NULL)
{
}

tscomponent::tscomponent(int s,double* t)
  : size(s)
{
  time=t;

  if (size) {
    try {
      vals   = new double[size];
      sigmas = new double[size];
      flags  = new FLAG[size];
    } catch (std::bad_alloc& e) {
      std::cerr<<"\n*** Cannot allocate memory!";
    }
  }

  for (int i=0;i<size;i++) *(flags+i)=OK;  
}

tscomponent::~tscomponent()
{
  if (size) {
    delete[] vals;
    delete[] sigmas;
    delete[] flags;
  }
}

void tscomponent::resize(int s)
{
  if (s!=size && size) {
      delete[] vals;
      delete[] sigmas;
      delete[] flags;
  }
  
  size=s;

  if (size) {
    try {
      vals   = new double[size];
      sigmas = new double[size];
      flags  = new FLAG[size];
    } catch (std::bad_alloc& e) {
      std::cerr<<"\n*** Cannot allocate memory!";
    }
  }

  for (int i=0;i<size;i++) *(flags+i)=OK; 
  
}
/*
int tscomponent::getrealsize() const
{
  int rsize=0;

  for (int i=0;i<size;i++)
    if (*(flags+i)==0||*(flags+i)==3)
      rsize++;

  return rsize;  
}
*/
void tscomponent::mark(const FLAG& f,const int& srt,const int& sop)
{
  for (int i=srt;i<sop;i++)
    *(flags+i)=f;
  return;
}

void tscomponent::cleanoutliers()
{
  for (int i=0;i<size;i++)
    if (flags[i]==OUTLIER) *(flags+i)=OK;
  return;
}

int tscomponent::average(double& mean,double& smean) const
{
/* compute average value and variance of average */
  int i=-1,j;
  double k,gain;

  do {
    i++;
    mean  = vals[i];
    smean = sigmas[i];
  } while (flags[i]!=OK);

  j=1;
  for (i=i;i<size;i++) {
    if (flags[i]!=OUTLIER) {
      j++;
      k=1.0e00/((double)j+1.0e00);
      gain=(vals[i]-mean);
      mean+=k*gain;
      smean+=k*gain*gain;
      smean*=(1.0-k);
    }
  } 

  return j;
}

int tscomponent::cleansize() const
{
  int cs=0;

  for (int i=0;i<size;i++)
     if (flags[i]!=OUTLIER) cs++;

  return cs;
}

/*
** =============================================================================
** TIMESERIES CLASS
** =============================================================================
*/

timeseries::timeseries(std::vector<tspoint>& p,int s)
  : x(s),y(s),z(s)
{
   try {
    time   = new double[s];
    soltype = new short int[s];
  } catch (std::bad_alloc& e) {
    std::cerr<<"\n*** Cannot allocate memory!";
  }

  for (int i=0;i<s;i++) {
    time[i] = p[i].time;
    soltype[i] = 1;
    x.vals[i] = p[i].f;
    x.sigmas[i] = p[i].sf;
    y.vals[i] = p[i].l;
    y.sigmas[i] = p[i].sl;
    z.vals[i] = p[i].u;
    z.sigmas[i] = p[i].su;
  }

  x.time=time;
  y.time=time;
  z.time=time;
}

timeseries::timeseries(const std::string& name,const std::string& sta,bool userapid)
{
  station=sta;
  std::vector<tspoint> p;
  int tssize;

  std::ifstream fin (name.c_str());
  if (!fin.is_open()) {
    std::cerr<<"\nCould not open time-series file: "<<name;
    return;
  } else {
    tssize = readingcts(fin,p,userapid);
  }

  if (tssize) {
    x.resize(tssize);
    y.resize(tssize);
    z.resize(tssize);
    try {
      time   = new double[tssize];
      soltype = new short int[tssize];
    } catch (std::bad_alloc& e) {
      std::cerr<<"\n*** Cannot allocate memory!";
    }
    for (int i=0;i<tssize;i++) {
      time[i] = p[i].time;
      soltype[i]=p[i].flag;
      x.vals[i] = p[i].f;
      x.sigmas[i] = p[i].sf;
      y.vals[i] = p[i].l;
      y.sigmas[i] = p[i].sl;
      z.vals[i] = p[i].u;
      z.sigmas[i] = p[i].su;
    }
    x.time=time;
    y.time=time;
    z.time=time;
  } else {
    std::cerr<<"\nTime-series file: "<<name<<" has size =0";
  }

}

timeseries::~timeseries()
{
  delete[] time;
  delete[] soltype;
}

int timeseries::locate(const double& tm) const
{
  int s(size()),i;

  i=qslocate(time,s,tm);

  return i;
}

void timeseries::print(std::string nf,int c,bool usegmt) const
{
  int i,s=size();

  switch (usegmt) {
    /* DO NOT USE GMT FORMAT; DATES ARE PRINTED AS DECIMAL YEARS */
    case (0):
      if (nf=="") {/* print to stream */
        switch (c) {
          case (0):
            for (i=0;i<s;i++) {
              printf("\n%10.5f %8.5f %8.6f %1i %1i",time[i],x.vals[i],x.sigmas[i],x.flags[i],soltype[i]);
              printf(" %8.5f %8.6f %1i",y.vals[i],y.sigmas[i],y.flags[i]);
              printf(" %8.5f %8.6f %1i",z.vals[i],z.sigmas[i],z.flags[i]);
            }
            break;
          case (1):
            for (i=0;i<s;i++) printf("\n%10.5f %8.5f %8.6f %1i %1i",time[i],x.vals[i],x.sigmas[i],x.flags[i],soltype[i]);
            break;
          case (2):
            for (i=0;i<s;i++) printf("\n%10.5f %8.5f %8.6f %1i %1i",time[i],y.vals[i],y.sigmas[i],y.flags[i],soltype[i]);
            break;
          case (3):
            for (i=0;i<s;i++) printf("\n%10.5f %8.5f %8.6f %1i %1i",time[i],z.vals[i],z.sigmas[i],z.flags[i],soltype[i]);
            break;
        }
      } else {/* print to file */
        FILE* pFile;
        pFile = fopen(nf.c_str(),"w");
	        if (!pFile) {
		        std::cerr<<"\n*** Could not write time-series file!";
	        } else {
                switch (c) {
                  case (0):
                    for (i=0;i<s;i++) {
                      printf("\n%10.5f %8.5f %8.6f %1i %1i",time[i],x.vals[i],x.sigmas[i],x.flags[i],soltype[i]);
                      printf(" %8.5f %8.6f %1i",y.vals[i],y.sigmas[i],y.flags[i]);
                      printf(" %8.5f %8.6f %1i",z.vals[i],z.sigmas[i],z.flags[i]);
                    }
                    break;
                  case (1):
                    for (i=0;i<s;i++) fprintf(pFile,"\n%10.5f %8.5f %8.6f %1i %1i",time[i],x.vals[i],x.sigmas[i],x.flags[i],soltype[i]);
                    break;
                  case (2):
                    for (i=0;i<s;i++) fprintf(pFile,"\n%10.5f %8.5f %8.6f %1i %1i",time[i],y.vals[i],y.sigmas[i],y.flags[i],soltype[i]);
                    break;
                  case (3):
                    for (i=0;i<s;i++) fprintf(pFile,"\n%10.5f %8.5f %8.6f %1i %1i",time[i],z.vals[i],z.sigmas[i],z.flags[i],soltype[i]);
                    break;
              }
	        	fclose(pFile);
		        printf("\n$$ t-series file written: %20s",nf.c_str());
	        }
      }
      break;
    /* USE GMT FORMAT; DATES ARE PRINTED AS DD-MM-YYYY */
    case (1):
      int iy,im,id;
      std::string sm;
      if (nf=="") {/* print to stream */
        switch (c) {
          case (0):
            for (i=0;i<s;i++) {
              fy2cal(time[i],iy,im,id);
              sm = m2m(im);
              printf("\n%02i-%3s-%4i %8.5f %8.6f %1i %1i",id,sm.c_str(),iy,x.vals[i],x.sigmas[i],x.flags[i],soltype[i]);
              printf(" %8.5f %8.6f %1i",y.vals[i],y.sigmas[i],y.flags[i]);
              printf(" %8.5f %8.6f %1i",z.vals[i],z.sigmas[i],z.flags[i]);
            }
            break;
          case (1):
            for (i=0;i<s;i++) {
              fy2cal(time[i],iy,im,id);
              sm = m2m(im);
              printf("\n%02i-%3s-%4i %8.5f %8.6f %1i %1i",id,sm.c_str(),iy,x.vals[i],x.sigmas[i],x.flags[i],soltype[i]);
            }
            break;
          case (2):
            for (i=0;i<s;i++) {
              fy2cal(time[i],iy,im,id);
              sm = m2m(im);
              printf("\n%02i-%3s-%4i %8.5f %8.6f %1i %1i",id,sm.c_str(),iy,y.vals[i],y.sigmas[i],y.flags[i],soltype[i]);
            }
            break;
          case (3):
            for (i=0;i<s;i++) {
              fy2cal(time[i],iy,im,id);
              sm = m2m(im);
              printf("\n%02i-%3s-%4i %8.5f %8.6f %1i %1i",id,sm.c_str(),iy,z.vals[i],z.sigmas[i],z.flags[i],soltype[i]);
            }
            break;
        }
      } else {/* print to file */
        FILE* pFile;
        pFile = fopen(nf.c_str(),"w");
	if (!pFile) {
	        std::cerr<<"\n*** Could not write time-series file!";
	} else {
                switch (c) {
          case (0):
            for (i=0;i<s;i++) {
              fy2cal(time[i],iy,im,id);
              sm = m2m(im);
              fprintf(pFile,"\n%02i-%3s-%4i %8.5f %8.6f %1i",id,sm.c_str(),iy,x.vals[i],x.sigmas[i],x.flags[i]);
              fprintf(pFile," %8.5f %8.6f %1i",y.vals[i],y.sigmas[i],y.flags[i]);
              fprintf(pFile," %8.5f %8.6f %1i",z.vals[i],z.sigmas[i],z.flags[i]);
            }
            break;
          case (1):
            for (i=0;i<s;i++) {
              fy2cal(time[i],iy,im,id);
              sm = m2m(im);
              fprintf(pFile,"\n%02i-%3s-%4i %8.5f %8.6f %1i %1i",id,sm.c_str(),iy,x.vals[i],x.sigmas[i],x.flags[i],soltype[i]);
            }
            break;
          case (2):
            for (i=0;i<s;i++) {
              fy2cal(time[i],iy,im,id);
              sm = m2m(im);
              fprintf(pFile,"\n%02i-%3s-%4i %8.5f %8.6f %1i %1i",id,sm.c_str(),iy,y.vals[i],y.sigmas[i],y.flags[i],soltype[i]);
            }
            break;
          case (3):
            for (i=0;i<s;i++) {
              fy2cal(time[i],iy,im,id);
              sm = m2m(im);
              fprintf(pFile,"\n%02i-%3s-%4i %8.5f %8.6f %1i %1i",id,sm.c_str(),iy,z.vals[i],z.sigmas[i],z.flags[i],soltype[i]);
            }
            break;
              }
	        	fclose(pFile);
		        printf("\n$$ t-series file written: %20s",nf.c_str());
	        }
      }
      break;
  }

  return;
}

void timeseries::average(double* m,int c) const
{
  /* compute average value and variance of average */
  switch (c) {
    case (0):
      x.average(m[0],m[1]);
      y.average(m[2],m[3]);
      z.average(m[4],m[5]);
      break;
    case (1):
      x.average(m[0],m[1]);
      break;
    case (2):
      y.average(m[0],m[1]);
      break;
    case (3):
      z.average(m[0],m[1]);
      break;
  }

  return;
}
void timeseries::ell2top(const double& lat0,const double& lon0,const double& hgt0)
{
/*
ellispsoidal to topocentric
*/
  int SIZE = size();

  double dlat,dlon,dhgt;
  double p=mesitomi<double>(lat0);
  double r=parkik<double>(lat0);

  for (int i=0;i<SIZE;i++) {
    dlat=x.vals[i]-lat0;/* radians */
    dlon=y.vals[i]-lon0;/* radians */
    dhgt=z.vals[i]-hgt0;/* meters */
    x.vals[i]=dlat*p;
    y.vals[i]=dlon*r;
    z.vals[i]=dhgt;
  }

  return;
}

void timeseries::ell2top()
{
/*
ellispsoidal to topocentric
*/
  double m[6];
  average(m,0);

  ell2top(m[0],m[2],m[4]);

  return;
}

void timeseries::readindisc(const std::string& name)
{
/*
** read in and apply discontinuities. The file must be in the form:
** YEAR MONTH DAY [FLAG]
** lines starting with '#' are considered as comments
**  DISC=3
**  VELC=4
**  ERQT=5
*/
  std::ifstream fin (name.c_str());

  if (!fin.is_open()) {
    std::cerr<<"\nCould not open discontinuity file: "<<name;
    return;
  }
  
  std::string line;
  std::vector<std::string> strv;
  int iy,im,id;
  double t;
  
  getline(fin,line);
  if (!line.size()||fin.eof()) {
    fin.close();
    return;
  }

  while (!fin.eof()) {
    if (line[0]!='#') {/* if not comment line */
      strv = lineToStringVector(line,' ');
      if (strv.size()<3) {/* must have at least 3 valid fields */
        std::cout<<"\n### skipping line:"<<line;
        fin.close();
        applydisc();
        return;
      } else if (strv.size()==3) {/* fileds are ok in number; assign date */
        iy  = ::atoi(strv[0].c_str());
        im  = ::atoi(strv[1].c_str());
        id  = ::atoi(strv[2].c_str());
        if (iauCal2jd(iy,im,id,t)) {/* convert to Julian date */
          std::cout<<"\n### invalid line:"<<line;/* error; apply and return */
          applydisc();
          fin.close();
          return;
        } else {/* add the date to the discontinuity vector */
          discontinuities.push_back(t);
        }
      } else { //* fileds are ok in number, flag exists; assign date */
	iy  = ::atoi(strv[0].c_str());
	im  = ::atoi(strv[1].c_str());
	id  = ::atoi(strv[2].c_str());
	short int flg = ::atoi(strv[3].c_str());
	if (iauCal2jd(iy,im,id,t)) {/* convert to Julian date */
		std::cout<<"\n### invalid line:"<<line;/* error; apply and return */
		applydisc();
		fin.close();
		return;
	} else {/* add the date to the discontinuity vector */
		switch (flg) {
			case (3):
        printf("\n$$ Introduced discontinuity at: %02i-%02i-%04i",id,im,iy);
				discontinuities.push_back(t);
				break;
			case (4):
        printf("\n$$ Introduced velocity change at: %02i-%02i-%04i",id,im,iy);
				velchanges.push_back(t);
				break;
			case (5):
        printf("\n$$ Introduced earthquake event at: %02i-%02i-%04i",id,im,iy);
				earthquakes.push_back(t);
				break;
			default:
				std::cerr<<"\n*** Unknown event flag";
				break;
			}
		}
	}
    }
    getline(fin,line);
  }

  /* now apply them (i.e. mark the data-points).
  WARNING 
  discontinuities and time arrays must have the same date format (
  both decimal year or noth fractional year)!!
  */
  applydisc();

  /* done */
  return;
}

void timeseries::applydisc()
{
/*
apply the discontinuities contained in the objects
discontinuity vector
!!.  WARNING .!!
discontinuities and time arrays must have the same date format (
both decimal year or noth fractional year)!!
*/

  int index;

  /* apply discontinuities */
  for (int i=0;i<dsize();i++) {
    index=locate(discontinuities[i]);
    if (index>=0) {
      x.flags[index]=DISC;
      y.flags[index]=DISC;
      z.flags[index]=DISC;
    }
  }

  /* apply velocity changes */
  for (int i=0;i<vsize();i++) {
	  index=locate(velchanges[i]);
	  if (index>=0) {
		  x.flags[index]=VELC;
		  y.flags[index]=VELC;
		  z.flags[index]=VELC;
	  }
  }

  /* apply earthquackes */
  for (int i=0;i<esize();i++) {
	  index=locate(earthquakes[i]);
	  if (index>=0) {
		  x.flags[index]=ERQT;
		  y.flags[index]=ERQT;
		  z.flags[index]=ERQT;
	  }
  }

  return;
}

int timeseries::window(const int& index,const double& window,int& start,int& stop)
{
	int SIZE = size();
	int i;
	double limitl (time[index]-window/2.0);
	double limitu (time[index]+window/2.0);
	//printf("\n[%10.5f < %10.5f < %10.5f] ~",limitl,time[index],limitu);

	if (limitu<limitl)
		return -1;

	start = stop = index;

	for (i=start;i>0;i--)
		if (limitl<=time[i]) start=i;
		else break;

	for (i=stop;i<SIZE;i++)
		if (limitu>time[i]) stop=i;
		else break;


	//printf(" [%10.5f < %10.5f < %10.5f]",time[start],time[index],time[stop]);
	return stop-start;
}

void timeseries::cleanoutliers(int c)
{
  switch (c) {
    case (0):
      for (int i=0;i<size();i++){
        if (x.flags[i]==OUTLIER) x.flags[i]=OK;
        if (y.flags[i]==OUTLIER) y.flags[i]=OK;
        if (z.flags[i]==OUTLIER) z.flags[i]=OK;
      }
      break;
    case (1):
      x.cleanoutliers(); break;
    case (2):
      y.cleanoutliers(); break;
    case (3):
      z.cleanoutliers(); break;
  }

  return;
}

void timeseries::cleansize(int& sx,int& sy,int& sz) const
{
/*
count sizes of all non-flaged elements
for each component
*/
	int SIZE = size(),i;
	sx=sy=sz=0;

	for (i=0;i<SIZE;i++){
		if (x.flags[i]!=OUTLIER) sx++;
		if (y.flags[i]!=OUTLIER) sy++;
		if (z.flags[i]!=OUTLIER) sz++;
	}
	return;
}

bool timeseries::makemat
  (double* A/*design matrix*/,double* l/*array of observations*/,
  const std::vector<double>& freqs/*array with frequency multipliers*/,
  const double& sigma0/*a-priori std. deviation*/,const int& rows/*number of rows (cleaned)*/
  ,const int& cols/*number of columns=parameters*/,double& t0/*central epoch*/,
  int c/*component:1=n,2=e,3=u*/,bool use_weights/*use weights of data points*/)
{
/*
this function will form the following two matrices:
A*sqrt(P) and l*sqrt(P), with Pi=sigma0/sigma[i],
if use_weights is on, else
A, l
All non-flaged elements of the timeseries will be
used.
rows must be the clean-size of element, cols is
the number of parameters.
*/
  int SIZE = size();
  int fsize = freqs.size();
  int dcsize = dsize();
  int vcsize = vsize();
  int j,i,k;
  double* v,*s;
  FLAG* f;
  
  /* check that number of columns is correct */
  if (cols!=2+dcsize+fsize*2+vcsize) {
    std::cerr<<"\n*** Wrong number of columns given!";
    std::cerr<<"\n***cols="<<cols<<", parameters="<<2<<"+"<<dcsize<<"+"<<(fsize*2)<<"+"<<vcsize;
    return false;
  }

  /* if central epoch given = 0, compute the average */
  if (!t0) t0 = qaverage(time,SIZE);

  /* assign pointers, depending on component */
  switch (c) {
    case (1):
      v=x.vals;
      s=x.sigmas;
      f=x.flags;
      break;
    case (2):
      v=y.vals;
      s=y.sigmas;
      f=y.flags;
      break;
    case (3):
      v=z.vals;
      s=z.sigmas;
      f=z.flags;
      break;
    default:
      std::cerr<<"\n*** Invalid component!";
      return false;
  }

  if (use_weights) {/* USE WEIGHTS OF OBSERVATIONS */
  /* the following will only work for a diagonal weight matrix */
    j=0;
    for (i=0;i<SIZE;i++) {// from the whole list of elements ...
      if (f[i]!=OUTLIER) { // only take the ones flagged as OK
	double dT (time[i]-t0);// delta time
	double weight (sigma0/s[i]);// weight of observation p(i) = sigma0**2 / sigmai**2
        l[j]=v[i]*weight;// (weighted) observed value 
        A[j]=1.0e00*weight;// first term (constant)
        A[j+rows]=dT*weight;// second term, for velocity
        for (k=0;k<dcsize;k++) {// add discoontinuity terms ...
          A[j+(2+k)*rows]=(time[i]>=discontinuities[k])?1.0e00*weight:0.0e00;
        }
        for (k=0;k<fsize;k++) {// add harmonic terms ...
          A[j+(2+dcsize+2*k)*rows]  =weight* ::sin(freqs[k]*DPI*dT);
          A[j+(2+dcsize+2*k+1)*rows]=weight* ::cos(freqs[k]*DPI*dT);
        }
        for (k=0;k<vsize();k++) {// add velocity change terms ...
		double delta = time[i]-velchanges[k];
		if (delta>=0.0) {
          		A[j+(2+dcsize+2*fsize+k)*rows]=delta*weight;
		} else {
			A[j+(2+dcsize+2*fsize+k)*rows] = 0.0e00;
		}
        }
        j++;
      }
    }
  } else {
    j=0;
    for (i=0;i<SIZE;i++) {// from the whole list of elements ...
      if (f[i]!=OUTLIER) { // only take the ones flagged as OK
	      double dT (time[i]-t0);
        l[j]=v[i];// observed value 
        A[j]=1.0e00;// first term (constant)
        A[j+rows]=dT;// second term, for velocity
        for (k=0;k<dcsize;k++) {// add discoontinuity terms ...
          A[j+(2+k)*rows]=(time[i]>=discontinuities[k])?1.0e00:0.0e00;
        }
        for (k=0;k<fsize;k++) {// add harmonic terms ...
          A[j+(2+dcsize+2*k)*rows]  =::sin(freqs[k]*DPI*dT);
          A[j+(2+dcsize+2*k+1)*rows]=::cos(freqs[k]*DPI*dT);
        }
        for (k=0;k<vsize();k++) {// add velocity change terms ...
		double delta = time[i]-velchanges[k];
		if (delta>=0.0) {
			A[j+(2+dcsize+2*fsize+k)*rows]=delta;
		} else {
			A[j+(2+dcsize+2*fsize+k)*rows] = 0.0e00;
		}
        }
        j++;
      }
    }
  }

/* DEBUG MATRIX PRINTING */
/* ===================== */
//int jj=0;
//for (jj=0;jj<j;jj++) printf("\n%15.10f",l[jj]);
/*
for (i=0;i<j;i++){
	printf("\n");
	for (int k=0;k<cols;k++){
		printf(" %15.10f",A[j*k+i]);
	}
}
*/

  /* rows must be the same as j */
  return (j==rows && j+(2+dcsize+2*fsize+vsize()-1)*rows==rows*cols);
}

int timeseries::lsqrsolve(
    double* A/*design matrix*/,double* b,/*array of observations*/
    double* est/*array of etimated parameters*/,
    const int& rows,/*number of rows (cleaned)*/
    const int& cols,/*number of columns=parameters*/
    const std::vector<double>& freqs/*array with frequency multipliers*/,
    double& var/*estimated variance of fit*/,
    const double& sigma0,/* a-priori std. deviation of unit weight*/
    double* vcv,/* upper triangular part of the vcv matrix size of vcv = cols * (cols+1) / 2 */
    int c/*component:1=n,2=e,3=u*/,bool mark/*mark outliers*/)
{
	int i,j,ORSIZE,marked;
	double limit;
	double *bcopy, *acopy, *s;

  /* underdetermined system; quick exit */
  if (rows<=cols) {
    printf("\n*** Obsrvations <= Parameters!");
    return -2;
  }

	/* make a copy of the b and A matrices */
  try {
   bcopy = new double[rows];
   acopy = new double[rows*cols];
  } catch (std::bad_alloc& e) {
    std::cerr<<"\n*** Cannot allocate memory!";
    return -1;
  }
	for (i=0;i<rows;i++) bcopy[i]=b[i];
	for (i=0;i<rows*cols;i++) acopy[i]=A[i];

	/* solve the system via QR decomposition */
	int info = LAPACKE_dgels(LAPACK_COL_MAJOR,'N',rows,cols,1,acopy,rows,bcopy,rows);
	if (info>0) {
		std::cerr<<"\n*** Unable to perform least squares estimation";
		delete[] bcopy;
    		delete[] acopy;
		return -1;
	}

	/* make x vector (estimates) */
	for (i=0;i<cols;i++) est[i]=bcopy[i];

	/* compute residuals stored in b (i.e. res = 1.0*A*x-1.0*b ) */
	/* WARNIGN!! these residuals contain the square root of the weight matrix 
	 * meaning that res(i) = res(i) * sqrt(P(i,i))
	 * */
	cblas_dgemv(CblasColMajor,CblasNoTrans,rows,cols,1.0e00,A,rows,est,1,-1.0e00,b,1);

	/* 
	 * compute var-covar matrix (only upper triangular part) 
	 * and assign std. deviations of parameters
	 */
	ata(cols,rows,A,vcv);
	info = posdefinv(cols,vcv);
	if (!info) {
		std::cerr<<"\n*** Invalid computation of VCV matrix!";
	}

	/* compute a-posteriori variance (i.e. res**T * res / (rows-parameters)) */
	var = cblas_ddot(rows,b,1,b,1) / (double)(rows-cols);

  /* limit for outlier detection (3-sigma) */
  limit = 3.0e00*::sqrt(var);

  /* short report */
  #ifdef DEBUG
	  printf("\n+ StdDev of unit weight is %10.3f",::sqrt(var));
  #endif

	/* assigned needed pointers */
	FLAG* fs;
	switch (c) {
		case (1):
			fs=x.flags;
			ORSIZE=x.size;
			s=x.sigmas;
			break;
		case (2):
			fs=y.flags;
			ORSIZE=y.size;
			s=y.sigmas;
			break;
		case (3):
			fs=z.flags;
			ORSIZE=z.size;
			s=z.sigmas;
			break;
		default:
			std::cerr<<"\ninvalid component";
			delete[] bcopy;
			return -1;
	}

  	/* mark residuals */
	//marked = markoutliers(b,::sqrt(var),c);
	j=marked=0;
	if (mark) {
	  for (i=0;i<ORSIZE;i++) {
		  if (fs[i]!=OUTLIER) {
			  b[j] *= s[i]/sigma0;
			  if (::fabs(b[j])>limit) {
				  fs[i]=OUTLIER;
				  marked++;
			  }
			  // also de-weight residuals
			  //b[j] *= s[i]/sigma0;
			  j++;
		  }
	  }
  	}

	/* deallocate memory */
	delete[] bcopy;
	delete[] acopy;

	/* all done */
	return marked;
}

void timeseries::modelfitc(
  double* estimates,/*array of etimated parameters*/
  const int& c/*component:1=n,2=e,3=u*/,
  const int& model_pars,/*number of columns=parameters*/
  const std::vector<double>& freqs/*array with frequency multipliers*/,
  const double& s0/*a-priori std. deviation*/,bool iterative/*perform iterations*/,
  std::string resfile/*name of residual file*/,
  bool usegmt/*print gmt formated dates (DD-MMM-YYYY)for model and residuals*/)
{
  int outliers=0,deleted=0;
  int latsize,lonsize,hgtsize;
  int SIZE = size();
  int ROWS;
  int COLS=model_pars;
  int PARS=model_pars;
  double percent=0;
  int iteration=0;
  double* A, *b;
  double fit_variance;
  double sigma0=s0;
  bool keepon;
  double T0=0.0e00;/*if starting with 0, T0 will be assigned in makemat function*/
  FLAG* flg;

  /* allocate memmory for vcv matrix */
  double* vcv = new double[PARS * (PARS+1) / 2];

  /* clean all previously marked data points, apart from discontinuities */
  cleanoutliers(c);

  do { /* ITERATIVELY SOLVE */

    /* get real sizes, and assign flags */
    cleansize(latsize,lonsize,hgtsize);
    switch (c) {
      case (1):
        ROWS=latsize;
        flg=x.flags;
        break;
      case (2):
        ROWS=lonsize;
        flg=y.flags;
        break;
      case (3):
        ROWS=hgtsize;
        flg=z.flags;
        break;
      default:
        std::cerr<<"\n*** Invalid component!";
        return;
    }

    /* allocate memmory */
    try {
      A = new double[ROWS*COLS];
      b = new double[ROWS];
    } catch (std::bad_alloc& e) {
      std::cerr<<"\n*** Cannot allocate memory!";
      break;
    }

    /* form matrices (T0 assigned in first oteration!!) */
    deleted = makemat(A,b,freqs,sigma0,ROWS,COLS,T0,c);
    if (!deleted) {
      std::cerr<<"\n*** Something went wrong in matrix filling!";
      std::cerr<<"\n*** Model parameters are "<<model_pars<<", COLS="<<COLS;
      delete[] A;
      delete[] b;
      break;
    }

    /* solve least squares, mark outliers */
    deleted=lsqrsolve(A,b,estimates,ROWS,COLS,freqs,fit_variance,s0,vcv,c,true);
    if (deleted<0) {
      delete[] A;
      delete[] b;
      std::cerr<<"\n*** Could not solve the system!";
      break;
    }
    /* !. IMPORTANT .! b vector now contains residuals */

    /* update outliers, percentge, iteration counter */
    outliers+=deleted; 
    percent=outliers*100.0/SIZE;
    iteration++;

    /* new sigma 0 for next iteration */
    //sigma0=::sqrt(fit_variance);
    sigma0=0.001;

    /* decide if to perform more iterations */
    keepon = (deleted!=0&&percent<20&&iterative);

    /* if this was the last iteration, write the residuals and model! */
    if ((!keepon)&&resfile!="") {
	    FILE* pFile;
	    pFile = fopen(resfile.c_str(),"w");
	    if (!pFile) {
		    std::cerr<<"\n*** Could not write residual file!";
	    } else {
	    	int _size_ = size();
	    	int j=0,i;
	    	for (i=0;i<_size_;i++) {
			    if (flg[i]!=OUTLIER) {
				    if (!usegmt) {
		    			fprintf(pFile,"\n%10.5f %+10.6f %1i",time[i],b[j],soltype[i]);
				    } else {
					int iy,im,id;
					std::string sm;
					fy2cal(time[i],iy,im,id);
					sm = m2m(im);
					fprintf(pFile,"\n%02i-%3s-%4i %8.5f %1i",id,sm.c_str(),iy,b[j],soltype[i]);
				    }
				j++;
		    	}
	    	}
	    	fclose(pFile);
		printf("\n$$ residual file written: %20s",resfile.c_str());
	    }

      /* also print model */
      std::string modelfile = resfile;
      std::string resstr = "res"; 
      modelfile.replace(modelfile.find(resstr),resstr.length(),"mod");
      makemodel(time[0]-0.01,time[SIZE-1]+0.01,0.005,estimates,PARS,freqs,T0,c,modelfile,usegmt);
    }

  } while (keepon);/* TEST ITERATION CONDITION*/

  /* report results */
  printfit(PARS,estimates,vcv,s0,fit_variance,freqs,outliers,percent,c,iteration);

  /* deallocate memory */
  delete[] vcv;

  /* done */
  return;
}

void timeseries::printfit(
  const int& parameters/*number of etimated parameters*/,
  double* estimates/*array of estimated parameters*/,
  double* vcv/**/,
  const double& aprsigma/*a-priori std. deviation*/,
  const double& var0,/*a-posteriori variance of fit*/
  const std::vector<double>& freqs,/*array with frequency multipliers*/
  const int& outliers/*number of outliers and*/,const double& percent,/*their percentage*/
  const int& c/*component*/,const int& iterations/*number of iterations performed*/)
{
  /* variables and central epoch */
  double T0 = qaverage(time,size());
  int iy,im,id,i;
  //double fd;
  int fsize = freqs.size();
  std::string comp;

  /* Topocentric component name*/
  switch (c) {
    case (1):
      comp="north";break;
    case (2):
      comp="east";break;
    case (3):
      comp="up";break;
    default:
      comp="unknown";break;
  }

  /* assign std deviations of parameters */
  double* stdevs = new double[parameters];
  stdevs[0]=vcv[0];
  int j=0,add=2;
  i=0;
  do {
	  j+=add;
	  stdevs[i++] = ::sqrt(vcv[j])*var0;
	  add+=1;
  } while (i+1<parameters);

  /* central epoch to Gregorian Calender */
  // iauJd2cal(T0,iy,im,id,fd);

  printf("\n--------------------------------------------------");
  printf("\n| STATION NAME : %10s                      |",station.c_str());
  printf("\n| COMPONENT    : %5s                           |",comp.c_str());
  printf("\n+------------------------------------------------+");
  printf("\n| Number of data points : %10i             |",size());
  printf("\n| Number of model params: %10i             |",parameters);
  printf("\n| Number of iterations  : %10i             |",iterations);
  printf("\n| Reference sigma of fit: %10.4f (m)         |",sqrt(var0));
  printf("\n| A-Priori sigma        : %10.4f (m)         |",aprsigma);
  printf("\n| Number of outliers    : %5i (~%4.1f%%)         |",outliers,percent);
  printf("\n+------------------------------------------------+");
  printf("\n| Model Parameters                               |");
 /*
  printf("\n| Central Epoch : %10.3f (Julian Date)       |",T0);
  printf("\n|            or : %4i-%02i-%02i, %03.1f (hours)        |",iy,im,id,fd);
  */
  printf("\n| Central Epoch : %10.3f (Years)             |",T0);
  printf("\n| Velocity (initial):  %+6.2f (mm/yr) +- %6.3f  |",estimates[1]*1000.0e00,stdevs[1]*1000.0);

  /* print discontinuities */
  if (!dsize()) {
	  printf("\n| Discontinuities  (None)                        |");
  } else {
	  printf("\n| Discontinuities                                |");
	  for (i=0;i<dsize();i++) {
      fy2cal(discontinuities[i],iy,im,id);
      printf("\n| Epoch : %10.3f (%4i-%02i-%02i) ->             |",discontinuities[i],iy,im,id);
      printf("\n|        -> %+7.4f (m) +- %7.5f               |",estimates[2+i],stdevs[2+i]);
    }
  }

  /* print harmonics */
  if (!fsize) {
    printf("\n| Frequency Harmonics (None)                     |");
  } else {
    printf("\n| Frequency Harmonics                            |");
    for (i=0;i<fsize;i++) {
      printf("\n| Period: %+3.2f * PI * DeltaT                    |",freqs[i]);
      printf("\n|        -> in-phase     :%+5.4f (mm) +- %7.4f|",estimates[2+dsize()+i]*1000.0,stdevs[2+dsize()+i]*1000.0);
      printf("\n|        -> out-of-phase :%+5.4f (mm) +- %7.4f|",estimates[2+dsize()+i+1]*1000.0,stdevs[2+dsize()+i+1]*1000.0);
    }
  }

  /* print velocity changes */
  if (!vsize()) {
    printf("\n| Velocity changes  (None)                       |");
  } else {
    printf("\n| Velocity changes                               |");
    for (i=0;i<vsize();i++) {
      fy2cal(velchanges[i],iy,im,id);
      printf("\n| Epoch : %10.3f (%4i-%02i-%02i)                |",velchanges[i],iy,im,id);
      printf("\n|        -> %+6.2f (mm/yr) +- %6.3f             |",estimates[2+dsize()+fsize*2+i]*1000.0,stdevs[2+i]*1000.0);
    }
  }
  
  printf("\n--------------------------------------------------");

  /* deallocate */
  delete[] stdevs;

  return;
}

void timeseries::makemodel(
  const double& tl/*lower time limit*/,const double& tu/*upper time limit*/,
  const double& dt/*time increment*/,
  double* estimates,/*array of estimated parameters*/
  const int& pars/*number of etimated parameters*/,
  const std::vector<double>& freqs,/*array with frequency multipliers*/
  const double& t0/*central epoch*/,int c/*component*/,
  std::string& modelfile/*name of output file*/,
  bool usegmt/*print date as DD-MMM-YYYY*/)
{
  int size = (tu-tl)/dt;
  int i,k;
  double temp=tl;
  int fsize = freqs.size();
  int dcsize = dsize();

  /* allocate memory A*x=l*/
  double* t = new double[size];
  double* A = new double[size*pars];
  double* l = new double[size];

  /* assign time vector */
  for (i=0;i<size;i++) {
    t[i] = temp;
    temp+=dt;
  }

  /* form A matrix */
  for (i=0;i<size;i++) {
    double dT (t[i]-t0);
    A[i]=1.0e00;// first term (constant)
    A[i+size]=dT;// second term, for velocity
    for (k=0;k<dcsize;k++) {// add discoontinuity terms ...
      A[i+(2+k)*size]=(t[i]>=discontinuities[k])?1.0e00:0.0e00;
    }
    for (k=0;k<fsize;k++) {// add harmonic terms ...
      A[i+(2+dcsize+2*k)*size]  =::sin(freqs[k]*DPI*dT);
      A[i+(2+dcsize+2*k+1)*size]=::cos(freqs[k]*DPI*dT);
    }
    for (k=0;k<vsize();k++) {// add velocity change terms ...
	    double delta = t[i]-velchanges[k];
	    if (delta>=0.0) {
			A[i+(2+dcsize+2*fsize+k)*size]=delta;
	    } else {
		    A[i+(2+dcsize+2*fsize+k)*size]=0.0e00;
	    }
    }
  }

  /* form the model: A * x = y , where y-->l */
  cblas_dgemv(CblasColMajor,CblasNoTrans,size,pars,1.0e00,A,size,estimates,1,0.0e00,l,1);

  /* write model to file */
  FILE* pFile;
  pFile = fopen(modelfile.c_str(),"w");
  if (!pFile) {
    std::cerr<<"\n*** Could not write residual file!";
  } else {
    if (!usegmt) {
      for (i=0;i<size;i++) fprintf(pFile,"\n%10.5f %+10.6f %1i",t[i],l[i],soltype[i]);
    } else {
      int iy,im,id;
      std::string sm;
      for (i=0;i<size;i++) {
        fy2cal(t[i],iy,im,id);
        sm = m2m(im);
        fprintf(pFile,"\n%02i-%3s-%4i %8.5f %1i",id,sm.c_str(),iy,l[i],soltype[i]);
      }
    }
    fclose(pFile);
		printf("\n$$ model file written   : %20s",modelfile.c_str());
  }

  /* deallocate */
  delete[] l;
  delete[] A;
  delete[] t;

  /* all done */
  return; 
}
void timeseries::smprint(std::string nf,int c,bool usegmt) const
{
  int i,s=size();

  double* ave = new double[6];
  average(ave,0);
  printf("\nAverage x = %10.4f\ny = %10.4f\nz = %10.4f",ave[0],ave[2],ave[4]);

  switch (usegmt) {
    /* DO NOT USE GMT FORMAT; DATES ARE PRINTED AS DECIMAL YEARS */
    case (0):
      if (nf=="") {/* print to stream */
        switch (c) {
          case (0):
            for (i=0;i<s;i++) {
              printf("\n%10.5f %8.5f %8.6f %1i %1i",time[i],x.vals[i],x.sigmas[i],x.flags[i],soltype[i]);
              printf(" %8.5f %8.6f %1i",y.vals[i],y.sigmas[i],y.flags[i]);
              printf(" %8.5f %8.6f %1i",z.vals[i],z.sigmas[i],z.flags[i]);
            }
            break;
          case (1):
            for (i=0;i<s;i++) printf("\n%10.5f %8.5f %8.6f %1i %1i",time[i],x.vals[i],x.sigmas[i],x.flags[i],soltype[i]);
            break;
          case (2):
            for (i=0;i<s;i++) printf("\n%10.5f %8.5f %8.6f %1i %1i",time[i],y.vals[i],y.sigmas[i],y.flags[i],soltype[i]);
            break;
          case (3):
            for (i=0;i<s;i++) printf("\n%10.5f %8.5f %8.6f %1i %1i",time[i],z.vals[i],z.sigmas[i],z.flags[i],soltype[i]);
            break;
        }
      } else {/* print to file */
        FILE* pFile;
        pFile = fopen(nf.c_str(),"w");
	        if (!pFile) {
		        std::cerr<<"\n*** Could not write time-series file!";
	        } else {
                switch (c) {
                  case (0):
                    for (i=0;i<s;i++) {
                      printf("\n%10.5f %8.5f %8.6f %1i %1i",time[i],x.vals[i],x.sigmas[i],x.flags[i],soltype[i]);
                      printf(" %8.5f %8.6f %1i",y.vals[i],y.sigmas[i],y.flags[i]);
                      printf(" %8.5f %8.6f %1i",z.vals[i],z.sigmas[i],z.flags[i]);
                    }
                    break;
                  case (1):
                    for (i=0;i<s;i++) fprintf(pFile,"\n%10.5f %8.5f %8.6f %1i %1i",time[i],x.vals[i],x.sigmas[i],x.flags[i],soltype[i]);
                    break;
                  case (2):
                    for (i=0;i<s;i++) fprintf(pFile,"\n%10.5f %8.5f %8.6f %1i %1i",time[i],y.vals[i],y.sigmas[i],y.flags[i],soltype[i]);
                    break;
                  case (3):
                    for (i=0;i<s;i++) fprintf(pFile,"\n%10.5f %8.5f %8.6f %1i %1i",time[i],z.vals[i],z.sigmas[i],z.flags[i],soltype[i]);
                    break;
              }
	        	fclose(pFile);
		        printf("\n$$ t-series file written: %20s",nf.c_str());
	        }
      }
      break;
    /* USE GMT FORMAT; DATES ARE PRINTED AS DD-MM-YYYY */
    case (1):
      int iy,im,id;
      std::string sm;
      if (nf=="") {/* print to stream */
        switch (c) {
          case (0):
            for (i=0;i<s;i++) {
              fy2cal(time[i],iy,im,id);
              sm = m2m(im);
              printf("\n%02i-%3s-%4i %8.5f %8.6f %1i %1i",id,sm.c_str(),iy,x.vals[i],x.sigmas[i],x.flags[i],soltype[i]);
              printf(" %8.5f %8.6f %1i",y.vals[i],y.sigmas[i],y.flags[i]);
              printf(" %8.5f %8.6f %1i",z.vals[i],z.sigmas[i],z.flags[i]);
            }
            break;
          case (1):
            for (i=0;i<s;i++) {
              fy2cal(time[i],iy,im,id);
              sm = m2m(im);
              printf("\n%02i-%3s-%4i %8.5f %8.6f %1i %1i",id,sm.c_str(),iy,x.vals[i],x.sigmas[i],x.flags[i],soltype[i]);
            }
            break;
          case (2):
            for (i=0;i<s;i++) {
              fy2cal(time[i],iy,im,id);
              sm = m2m(im);
              printf("\n%02i-%3s-%4i %8.5f %8.6f %1i %1i",id,sm.c_str(),iy,y.vals[i],y.sigmas[i],y.flags[i],soltype[i]);
            }
            break;
          case (3):
            for (i=0;i<s;i++) {
              fy2cal(time[i],iy,im,id);
              sm = m2m(im);
              printf("\n%02i-%3s-%4i %8.5f %8.6f %1i %1i",id,sm.c_str(),iy,z.vals[i],z.sigmas[i],z.flags[i],soltype[i]);
            }
            break;
        }
      } else {/* print to file */
        FILE* pFile;
        pFile = fopen(nf.c_str(),"w");
	if (!pFile) {
	        std::cerr<<"\n*** Could not write time-series file!";
	} else {
                switch (c) {
          case (0):
            for (i=0;i<s;i++) {
              fy2cal(time[i],iy,im,id);
              sm = m2m(im);
              fprintf(pFile,"\n%02i-%3s-%4i %8.5f %8.6f %1i",id,sm.c_str(),iy,x.vals[i],x.sigmas[i],x.flags[i]);
              fprintf(pFile," %8.5f %8.6f %1i",y.vals[i],y.sigmas[i],y.flags[i]);
              fprintf(pFile," %8.5f %8.6f %1i",z.vals[i],z.sigmas[i],z.flags[i]);
            }
            break;
          case (1):
            for (i=0;i<s;i++) {
              fy2cal(time[i],iy,im,id);
              sm = m2m(im);
              double limit = 3.0 * ave[1];
              if (fabs(ave[0]-x.vals[i])<limit) {
                fprintf(pFile,"\n%02i-%3s-%4i %8.5f %8.6f %1i %1i",id,sm.c_str(),iy,x.vals[i],x.sigmas[i],x.flags[i],soltype[i]);
              } else {
		printf("\nSkipping point, mean=%7.4f, value=%7.4f, std.dev=%7.4f",ave[0],x.vals[i],ave[1]);
	      }
            }
            break;
          case (2):
            for (i=0;i<s;i++) {
              fy2cal(time[i],iy,im,id);
              sm = m2m(im);
              double limit = 3.0 * ave[3];
              if (fabs(ave[2]-y.vals[i])<limit) {
                fprintf(pFile,"\n%02i-%3s-%4i %8.5f %8.6f %1i %1i",id,sm.c_str(),iy,y.vals[i],y.sigmas[i],y.flags[i],soltype[i]);
              }
            }
            break;
          case (3):
            for (i=0;i<s;i++) {
              fy2cal(time[i],iy,im,id);
              sm = m2m(im);
              double limit = 3.0 * ave[5];
              if (fabs(ave[4]-z.vals[i])<limit) {
              fprintf(pFile,"\n%02i-%3s-%4i %8.5f %8.6f %1i %1i",id,sm.c_str(),iy,z.vals[i],z.sigmas[i],z.flags[i],soltype[i]);
              }
            }
            break;
              }
	        	fclose(pFile);
		        printf("\n$$ t-series file written: %20s",nf.c_str());
	        }
      }
      break;
  }

  return;
}
