#ifndef _VARTOOLS_HPP
#define _VARTOOLS_HPP

#include <cmath>

// grs-80 parameters
#define grs80_a   6378137.
#define grs80_b   6356752.31414
#define grs80_f   0.003352810681225
#define grs80_e2  0.006694380022903
#define grs80_e2t 0.006739496775592
#define grs80_ba  0.996647189318816

static const double DPI (3.141592653589793238462643);
static const double D2PI(6.283185307179586476925287);
static const double DD2R(DPI/180.0e00);

template <class T>
T hexdegreestoradians(const int& id,const int& im,const T& sec)
{
  return (
    static_cast<T>(id) +
    static_cast<T>(im)/60.0e00 +
    sec/3600.0e00 ) * DD2R;
}

template <class T>
int qslocate(const T* data,const int& size,const T& x)
/*
return the index of the element <= x;
*/
{
	if (x<data[0]) {
		return -1;
	} else if (x>data[size-1]) {
		return size-1;
	} else {
    int ju(size-1);
		int jl(0),jm;
		while (ju-jl>1) {
			jm = (ju+jl)>>1;
			if (x>=data[jm])
				jl=jm;
			else
				ju=jm;
		}
		return jl;
	}

}

template <class T>
T qaverage(const T* data,const int& size)
{
  double m (0.0e00);
  for (int i=0;i<size;i++)
    m+=data[i];
  return m / (T)size;
}

template<class T>
inline T mesitomi (const T& lat)
{
  T sinl2 = ::sin(lat);
  sinl2*=sinl2;
  return grs80_a*(1.0e00-grs80_e2) / ::pow(1.0e00-grs80_e2*sinl2,3/2);
}

template<class T>
inline T parkik (const T& lat)
{
  T sinl2 = ::sin(lat);
  sinl2*=sinl2;
  return ::cos(lat)*(grs80_a / ::sqrt(1.0e00-grs80_e2*sinl2));
}

template<class T>
int minindex(T* data,const int& start,const int& stop)
{
	int i,m(start);
	T minv = data[start];

	for (i=start;i<=stop;i++)
		if (data[i]<minv)
			m=i;
	return m;
}

template<class T>
int median(T* data,const int& start,const int& stop)
{
        int middle = (stop-start) / 2;
	int i,j;
	T val;
	
	T* copy =  new T[stop-start+1];
	
	for (i=0;i<middle;i++) { 
		j=minindex(data,start+i,stop);
		val = data[j];
		data[j] = data[start+i];
		data[start+i] = val;
	}
	
	delete[] copy;
	return j;
}

#endif
