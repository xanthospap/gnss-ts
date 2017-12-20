#include <algorithm>
// Eigen headers
#include "eigen3/Eigen/Core"
// see https://en.wikipedia.org/wiki/Non-negative_least_squares
// https://github.com/scipy/scipy/blob/master/scipy/optimize/nnls/nnls.f

void
h12(int mode, int lpivot, int l1, int m, double* u, int iue, double updouble* c, int ice, int icv, int ncv)
/*
C     SUBROUTINE H12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV)
C
C  CONSTRUCTION AND/OR APPLICATION OF A SINGLE
C  HOUSEHOLDER TRANSFORMATION..     Q = I + U*(U**T)/B
C
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1973 JUN 12, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C     ------------------------------------------------------------------
c                     Subroutine Arguments
c
C     MODE   = 1 OR 2   Selects Algorithm H1 to construct and apply a
c            Householder transformation, or Algorithm H2 to apply a
c            previously constructed transformation.
C     LPIVOT IS THE INDEX OF THE PIVOT ELEMENT.
C     L1,M   IF L1 .LE. M   THE TRANSFORMATION WILL BE CONSTRUCTED TO
C            ZERO ELEMENTS INDEXED FROM L1 THROUGH M.   IF L1 GT. M
C            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
C     U(),IUE,UP    On entry with MODE = 1, U() contains the pivot
c            vector.  IUE is the storage increment between elements.
c            On exit when MODE = 1, U() and UP contain quantities
c            defining the vector U of the Householder transformation.
c            on entry with MODE = 2, U() and UP should contain
c            quantities previously computed with MODE = 1.  These will
c            not be modified during the entry with MODE = 2.
C     C()    ON ENTRY with MODE = 1 or 2, C() CONTAINS A MATRIX WHICH
c            WILL BE REGARDED AS A SET OF VECTORS TO WHICH THE
c            HOUSEHOLDER TRANSFORMATION IS TO BE APPLIED.
c            ON EXIT C() CONTAINS THE SET OF TRANSFORMED VECTORS.
C     ICE    STORAGE INCREMENT BETWEEN ELEMENTS OF VECTORS IN C().
C     ICV    STORAGE INCREMENT BETWEEN VECTORS IN C().
C     NCV    NUMBER OF VECTORS IN C() TO BE TRANSFORMED. IF NCV .LE. 0
C            NO OPERATIONS WILL BE DONE ON C().
C     ------------------------------------------------------------------
*/
{
    if ((0 >= lpivot || lpivot >= l1) || l1 > m) return;
    double cl = std::abs(umat(0,lpivot));

    if (mode != 2) {
        // CONSTRUCT THE TRANSFORMATION.
        for (int j = l1; j < m; j++) {
            cl = std::max(std::abs(umat(0,j)), cl);
            if (cl <= 0e0) {
                return;
            } else {
                double clinv = 1e0 / cl;
                double sm = (umat(0,lpivot)*clinv)*(umat(0,lpivot)*clinv);
                for (
            }
}

// GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B,  COMPUTE AN
// N-VECTOR, X, THAT SOLVES THE LEAST SQUARES PROBLEM
// A * X = B  SUBJECT TO X .GE. 0 
void
nnls(const double* a, const double* b, int m, int n, double* x, int& mode)
{
    double *w, *zz;
    int    *index;
    try {
        w = new double[n];
        zz = new double[m];
        index = new int[n];
    } catch (std::bad_alloc&) {
        return;
    }
    constexpr double factor{0.01e0},
                     two{2e0},
                     zero{0e0};

    if ( m <= 0 || n <= 0 ) {
        mode = 2;
        return;
    }
    
    int iter = 0,
        itmax= 3*n;

    // INITIALIZE THE ARRAYS INDEX() AND X().
    for (int i=0; i<n; i++) {
        x[i] = zero;
        index[i] = i;
    }

    int iz2   = n,
        iz1   = 0,
        nsetp = 0,
        npp   = 0;

    //  MAIN LOOP BEGINS HERE
    //  QUIT IF ALL COEFFICIENTS ARE ALREADY IN THE SOLUTION. OR IF M COLS OF
    //  A HAVE BEEN TRIANGULARIZED.
    int j, izmax;
    double sm, wmax, asave;
    if ( !(iz1 > iz2 || nsetp > m)) {
        // COMPUTE COMPONENTS OF THE DUAL (NEGATIVE GRADIENT) VECTOR W().
        for (int iz=iz1; iz<iz2; iz++) {
            j = index[iz];
            sm = zero;
            for (int l=npp1; l<m; l++) {
                sm += amat(l,j)*b(l);
            }
            w[j] = sm;
        }
        // FIND LARGEST POSITIVE W(J).
        wmax = zero;
        for (int iz=iz1; iz<iz2; iz++) {
            j = index[iz];
            if (w[j] > wmax) {
                wmax  = w[j];
                izmax = iz;
            }
        }
        // IF WMAX .LE. 0. GO TO TERMINATION.
        // THIS INDICATES SATISFACTION OF THE KUHN-TUCKER CONDITIONS.
        if (wmax <= zero) { got to 350 } <-- TODO
        iz = izmax;
        j  = index(iz);

        // THE SIGN OF W(J) IS OK FOR J TO BE MOVED TO SET P.
        // BEGIN THE TRANSFORMATION AND CHECK NEW DIAGONAL ELEMENT TO AVOID
        // NEAR LINEAR DEPENDENCE.
        asave = amat(npp1, j);
        h12(1, npp1, npp1+1, m, amat(1,j),1,up,dummy,1,1,0)

    }
}
