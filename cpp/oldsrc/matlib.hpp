#ifndef _MATLIB_HPP
#define _MATLIB_HPP

/* compute the inner product of vectors v1 and v2 */
double innerproduct(const int& r,double* v1,double* v2);

/* 
 * get the pointer to the first element of column <col>
 * for a matrix <A> stored column-wise with <ROWS> rows
 */
inline double* getcolumn(const int& col,const int& ROWS,double* A);

/*
 * given a matrix <A>, conmpute <ATA> = transpose(A) * A.
 * A has dimensions ROWS * COLS. ATA has a total size of
 * COLS * (COLS + 1) / 2. Only the Upper part is computed and stored.
 * Meaning:
 * *  Two-dimensional storage of the symmetric matrix A:
 * *
 * *     a11 a12 a13 a14
 * *         a22 a23 a24
 * *             a33 a34     (aij = aji)
 * *                 a44
 * *
 * *  Packed storage of the upper triangle of A:
 * *
 * *  ATA = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]
 * *
 */
void ata(const int& COLS,const int& ROWS,double* A,double* ATA);

/*
 * compute the inverse of a positive-definet matrix <A>
 * with rows = cols = dim. Only the Upper-triangular
 * part of A is given, packed columnwise in a linear array
 * (see function ata(...)
 */
bool posdefinv(const int& dim,double* A);

/*
 * print a symmetric matrix stored column-wise
 */
void printsym(double* A,const int& dim);

/*
 * compute interquartile range
 */
double iqr(double* data,const int& start,const int* stop);
bool iqrstatistics(double* data,const int& start,const int& stop,double& q1,double& q2,double& q3,double& q4);


#endif
