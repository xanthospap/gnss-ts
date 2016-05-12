#include "tstools.hpp"
#include "vartools.hpp"
#include "cppbasic.hpp"

//# year(0) month(1) day(2) hour(3) doy(4) fractional_year(5)  n(6-8), e(9-11) u(12) sn(13) se(14) su(15) [1|0] #
// 1 -> FINAL
// 0 -> RAPID

bool readinline(const std::string& line,tspoint& ts,bool userapid=false);
void add(std::vector<tspoint>& vec,const tspoint& p);

int readingcts(std::ifstream& fin,std::vector<tspoint>& ts,bool userapid)
{
 	/* read in rad, rad, meters */
	std::string line;
	tspoint tsp;
	bool ok;

	getline(fin,line);
	while (!fin.eof()) {
		ok = readinline(line,tsp,userapid);// try to read line
		if (ok) {// if line read ok
			add (ts,tsp);
		}
		getline(fin,line);
	}// end of while

  #ifdef DEBUG
    printf("\n--> Read in time series of %6i points",(int)ts.size());
  #endif

  return ts.size();
}

bool readinline(const std::string& line,tspoint& ts,bool userapid)
{
	ts.flag=1;
	// check if line is not comment
	if (line[0]=='#')
		return false;
	// split line
	std::vector<std::string> svec=lineToStringVector(line,' ');
	if (svec.size()!=16&&svec.size()!=17) {
		#ifdef DEBUG
			std::cout<<"\n### skipping line:"<<line;
		#endif
		return false;
	}
	// resolve final/rapid flag if any
	if (svec.size()==17) {
		ts.flag=::atoi(svec[16].c_str());
	}
	if (userapid||(!userapid&&ts.flag)) {// is the solution type valid ?
		int deg,min,iy,im,id;
		double sec;
		// read in time
		//ts.time = ::atof(svec[5].c_str());
		iy  = ::atoi(svec[0].c_str());
		im  = ::atoi(svec[1].c_str());
		id  = ::atof(svec[2].c_str());
  		if (iauCal2jd(iy,im,id,ts.time)) {
    			std::cout<<"\n### skipping line:"<<line;
    			return false;
  		} 
		// read in latitude
		deg  = ::atoi(svec[6].c_str());
		min  = ::atoi(svec[7].c_str());
		sec  = ::atof(svec[8].c_str());
		ts.f=hexdegreestoradians<double>(deg,min,sec);
		// read in longtitude
		deg  = ::atoi(svec[9].c_str());
		min  = ::atoi(svec[10].c_str());
		sec  = ::atof(svec[11].c_str());
		ts.l=hexdegreestoradians<double>(deg,min,sec);
		// read in height
		ts.u = ::atof(svec[12].c_str());
		// read in sigmas
		ts.sf = ::atof(svec[13].c_str());
		ts.sl = ::atof(svec[14].c_str());
		ts.su = ::atof(svec[15].c_str());
		// done!
		return 1;
	} else {
		std::cout<<"\n### skipping line:"<<line;
		return false;
	}

}

void add(std::vector<tspoint>& vec,const tspoint& p)
{
	int size = (int)vec.size();
	std::vector<tspoint>::iterator it;

	/* list is empty */
	if (!size) vec.push_back(p);

	/* non-empty list: */
	if (vec[size-1].time<p.time) {/* in this case add point */
		vec.push_back(p);
	} else if (vec[size-1].time==p.time) {/* in this case replace point */
		vec[size-1] = p;
	} else {/* in this case search for the right position */
		/* quck case, insert/replace first element */
		it=vec.begin();
		if (vec[0].time>p.time) {/* insert at first position or ... */
			vec.insert(it,p);
		} else if (vec[0].time==p.time) {/* replace first position */
			vec[0]=p;
		} else {/* search list */
			unsigned int j(0);
		for (it=vec.begin();it!=vec.end();it++) {
			if (it->time<=p.time && (it+1)->time>p.time) {/* found right position ... */
				if (it->time==p.time) {
					*(it)=p;
					break;
				}/* insert */
				else {
					vec.insert(it+1,p);
					break;
				}     /* replace */
        		} else {
				j++;
			}
		}
		#ifdef DEBUG
			if ((int)j>=size-1) {
          			printf("\n*** Invalid replace");
				printf("\nTrying to fit point  : %9.3f",p.time);
				printf("\nFirst point in series: %9.3f",vec[0].time);
				printf("\nLast  point in series: %9.3f",vec[size-1].time);
			}
		#endif
		}
	}

	return;
}
