#include "matlib.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "/usr/include/lapacke.h"
//#include "/home/xanthos/include/lapacke.h"

double innerproduct(const int& r,double* v1,double* v2) 
{
	double s (v1[0]*v2[0]);
	for (int i=1;i<r;i++) s+= v1[i]*v2[i];
	return s;
}

inline double* getcolumn(const int& col,const int& ROWS,double* A)
{
	return A+col*ROWS;
}

void ata(const int& COLS,const int& ROWS,double* A,double* ATA)
{
	int i,j,k;
	double *rowv, *colv;
	
	k=0;
	for (i=0;i<COLS;i++) {
		colv = getcolumn(i,ROWS,A);
		for (j=0;j<=i;j++) {
			rowv = getcolumn(j,ROWS,A);
			ATA[k++] = innerproduct(ROWS,rowv,colv);
		}
	}
	
	return;
}

bool posdefinv(const int& dim,double* A)
{
 int info = LAPACKE_dpptrf(LAPACK_COL_MAJOR,'U',dim,A);

 if (info) return false;

 info = LAPACKE_dpptri(LAPACK_COL_MAJOR,'U',dim,A);

 return (!info);
}

void printsym(double* A,const int& dim)
{
	int i,j,k=0;
	for (i=0;i<dim;i++) {
		printf("\n");
		for (j=0;j<dim;j++) {
			if (j<i) {
				printf("      ");
			} else {
				printf(" %5.3f",A[k++]);
			}
		}
	}
	
	return;
}

double iqr(double* data,const int& start,const int& stop)
{
	int comp(const void* a, const void* b);

	int size = stop-start;
	int middle = size/2;
	double q1,q3,q4;
	double* sorted;
	int i;

	try {
		sorted = new double[size];
	} catch (std::bad_alloc& e) {
		std::cerr<<"\n*** Cannot allocate memory";
		return -999;
	}

	for (i=0;i<size;i++) *(sorted+i) = data[start+i];

	qsort (sorted,size,sizeof(double),comp);

	q1 = sorted[middle/2];
	q3 = sorted[middle+middle/2];

	q4 = q3-q1;

	delete[] sorted;

	return q4;
}

bool iqrstatistics(double* data,const int& start,const int& stop,
		double& q1,double& q2,double& q3,double& q4)
{
	int comp(const void* a, const void* b);
	
	int size = stop-start;
	int middle = size/2;
	double* sorted;
	int i;
	
	try {
		sorted = new double[size];
	} catch (std::bad_alloc& e) { 
		std::cerr<<"\n*** Cannot allocate memory";
		return false;
	}
	
	for (i=0;i<size;i++) *(sorted+i) = data[start+i];
	
	qsort (sorted,size,sizeof(double),comp);

	q1 = sorted[middle/2];
	q2 = sorted[middle];
	q3 = sorted[middle+middle/2];
	q4 = q3-q1;

	return true;
}


int comp(const void* a, const void* b)
{
	if (*(double*)a>*(double*)b)
		return 1;
	else if (*(double*)a<*(double*)b)
		return -1;
	return 0;
}
