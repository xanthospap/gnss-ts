#include <algorithm>
// Eigen headers
#include "eigen3/Eigen/Core"
// see https://en.wikipedia.org/wiki/Non-negative_least_squares
// https://github.com/scipy/scipy/blob/master/scipy/optimize/nnls/nnls.f
/*
void
nnls(Eigen::MatrixXd& A, Eigen::VectorXd& y, double epsilon, int m, int n)
{
    Eigen::VectorXd x = Eigen::VectorXd(n);
    for (auto i = 0; i < n; i++) x(i) = 0e0;

    std::vector<> P;
    
    std::vector<int> R;
    R.reserve(n);
    for (auto i = 0 ; i < n; i++) R.push_back(i+1);

    // w = A^T(y-Ax)
    auto w = A.transpose()*(y-A*x);

    // main loop
    int j,l;
    double d;
    while (R.size() > 0 && w.maxCoeff() > epsilon) {
        // Let j in R be the index of max(w) in w
        d = w.maxCoeff(&j);
        // Add j to P
        P.push_back(j);
        std::sort(P.begin(), P.end());
        // Remove j from R
        std::remove(R.begin(), R.end(), j);
        // Let AP be A restricted to the variables included in P
        // int k = P.size();
        Eigen::MatrixXd AP = Eigen::MatrixXd(m,P.size());
        l = 0;
        for (auto i : P) {
            AP.col(l) = A.col(i);
            ++l;
        }
        // Let s be vector of same length as x. Let sP denote the sub-vector 
        // with indexes from P, and let sR denote the sub-vector with indexes
        // from R.
        Eigen::VectorXd sP = Eigen::VectorXd(P.size());
        Eigen::VectorXd sR = Eigen::VectorXd(R.size());
        // Set sP = ((AP)ᵀ AP)−1 (AP)ᵀy
        sP = ((AP.transpose()*AP).ldlt().solve(Eigen::MatrixXd::Identity(k,k))) * (AP.transpose()*y);
        // Set sR to zero
        for (int i = 0; i < (int)R.size(); i++) sR(i) = 0e0;
        Eigen::VectorXd s = Eigen::VectorXd(n);
        int cinp=0, cinr=0
        for (int i = 0; i < n; i++) {
            if ( std::find(P.begin(), P.end(), i) != P.end() ) {
                s(i) = sP(i-cinp);
                ++cinp;
            } else {
                s(i) = sR(i-cinr);
                ++cinr;
            }
        }
        // While min(sP) ≤ 0:
        while( sP.minCoeff() <= 0e0 ) {
            // Let α = min(xi/xi - si) for i in P where si ≤ 0
            double a = std::numeric_limits<double>::max();
            for (auto i : P) {
                if (s(i) <= 0e0 && x(i)/(x(i)-s(i)) < a) a = x(i)/(x(i)-s(i));
            }
            // Set x to x + α(s - x)
            x = x + a*(s-x);
            // Move to R all indices j in P such that xj = 0
            for (auto i : P) {
                if (x(i) == 0e0) {
                    R.push_back(i);
                }
            }
        }
    }
}
*/

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
