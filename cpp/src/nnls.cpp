#include <cmath>
#include <iostream>
#include <numeric>
// Eigen headers
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/QR"


/// Construct the submatrix Ap of A, where each column of Ap is in the vector P.
///
/// For example, if P = {1,2,5,9,3}, then Ap will have the same rows as A, but
/// the 1st column of Ap will be the 2nd of A, the 2nd column of Ap will
/// be the 3rd of A, the 3rd column of Ap will be the 6th of A, etc...
/// The resulting matrix, will have size (A.rows x P.size)
auto
make_Ap(const Eigen::MatrixXd& A, const std::vector<int>& P)
noexcept
{
    // checks ...
    assert( P.size() );
    auto max_it = std::max_element(P.cbegin(), P.cend());
    assert( *max_it <= A.cols()-1 );

    int m = A.rows(),
        k = 0;
    Eigen::MatrixXd Ap = Eigen::MatrixXd(m, P.size());
    
    for (auto i : P) {
        Ap.col(k) = A.col(i);
        ++k;
    }
    
    return Ap;
}

auto
make_s(const Eigen::VectorXd& sp, const std::vector<int>& P,
    const std::vector<int>& R)
noexcept
{
    // checks ...
    std::vector<int> M;
    M.reserve(P.size() + R.size());
    std::merge(P.cbegin(), P.cend(), R.cbegin(), R.cend(),
                std::back_inserter(M));
    assert( (int)M.size() == M[M.size()-1]+1 );
    int k = 0;
    for (auto i : M) {
        assert( i == k );
        ++k;
    }

    int n = M.size();
    Eigen::VectorXd s = Eigen::VectorXd::Zero(n);
    k = 0;
    for (auto i : P) {
        s(i) = sp(k);
        ++k;
    }

    return s;
}

/// see https://www.mathworks.com/matlabcentral/fileexchange/3388-nnls-and-constrained-regression?focused=5051382&tab=function
/// and Wikipedia https://en.wikipedia.org/wiki/Non-negative_least_squares
auto
fnnls(const Eigen::MatrixXd& A, const Eigen::VectorXd& y, double tol = -1e0)
{
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    // initialize variables
    if (tol < 0e0) {
        tol = std::numeric_limits<double>::epsilon();
    }
    int m = A.rows();
    int n = A.cols();
    assert( y.rows() == m );

    std::vector<int> P,
                     R(n),
                     V(n);
    std::iota(R.begin(), R.end(), 0);
    VectorXd x = VectorXd::Zero(n);
    VectorXd w = VectorXd::Zero(n);
    w = A.transpose()*(y-A*x);

    // set up iteration criterion
    int iter  = 0,
        itmax = 30*n;

    // Main loop
    int j,k;
    while (R.size() && w.maxCoeff(&j) > tol) {
        P.push_back(j);
        R.erase(std::remove(R.begin(), R.end(), j), R.end());
        std::sort(P.begin(), P.end());
        std::sort(R.begin(), R.end());
        MatrixXd Ap = make_Ap(A, P);
        VectorXd sp = VectorXd(P.size());
        sp = Ap.colPivHouseholderQr().solve(y);
        VectorXd s  = VectorXd::Zero(n);
        s = make_s(sp, P, R);

        // Inner Loop
        double alpha {std::numeric_limits<double>::max()},
               atmp;
        while ( sp.minCoeff(&k) <= 0e0 && iter < itmax ) {
            for (auto i : P) {
                if (s(i) <= 0e0) {
                    atmp = x(i)/(x(i)-s(i));
                    if ( atmp < alpha ) {
                        alpha = atmp;
                    }
                }
            }
            x = x + alpha*(s-x);
            V.clear();
            for (auto i : P) {
                if (x(i) == 0e0) {
                    R.push_back(i);
                    V.push_back(i);
                }
            }
            for (auto i : V) {
                std::remove(P.begin(), P.end(), i);
            }
            P.erase(P.begin()+V.size(), P.end());
            std::sort(P.begin(), P.end());
            std::sort(R.begin(), R.end());
            Ap = make_Ap(A, P);
            sp = Ap.colPivHouseholderQr().solve(y);
            s = make_s(sp, P, R);
            ++iter;
        }
        x = s;
        w = A.transpose()*(y-A*x);
    }

    return x;
}
// http://digitalassets.lib.berkeley.edu/sdtr/ucb/text/394.pdf
auto
bvls(int key, const Eigen::MatrixXd& A, const Eigen::VectorXd& b, 
    const Eigen::VectorXd& bl, const Eigen::VectorXd& bu,
    double tol = -1e0)
{
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    // initialize variables
    if (tol < 0e0) {
        tol = std::numeric_limits<double>::epsilon();
    }
    int m = A.rows();
    int n = A.cols();
    assert( y.rows() == m );
    assert( bu.rows() == bl.rows() && bl.rows() == n );

    // Step 1. Initialize everything free and bound sets, initial values, etc
    int mm  = std::min(m,n);
    int mm1 = mm+1;
    int jj  = 0;
    int ifrom5 = 0;
    
    // check consistency of given bounds
    double bdiff = 0e0;
    for (int j = 0; j < n; j++) {
        bdiff = std::max(bdiff, bu(j)-bl(j));
        if ( bl(j) > bu(j) ) {
            std::cerr<<"\nInconsistent bounds in BVLS; parmeter num = "<<i;
            return 1;
        }
    }
    if (bdiff == 0e0) {
        std::cerr<<"\nNo free variables in BVLS; check input bounds.";
        return 1;
    }

    /*
     * In a fresh initialization (key = 0) bind all variables at their lower 
     * bounds. If (key != 0), use the supplied istate vector to initialize 
     * the variables. istate(n+1) contains the number of bound variables. The
     * absolute values of the first nbound=istate(n+1) entries of istate are
     * the indices of the bound variables. The sign of each entry determines 
     * whether the indicated variable is at its upper (positive) or lower 
     * (negative) bound.
     */
    int nbound, nact;
    if ( !key ) {
        nbound = n;
        nact   = 0;
        for (int j = 0; j < nbound; j++) istate(j) = -j;
    } else {
        nbound = istate(n);
    }
    nact = n-nbound;
    if (nact > mm) {
        std::cerr<<"\nToo many free variables in BVLS starting solution.";
        return 1;
    }
    for (int k = 0; k < nbound; k++) {
        int j = std::abs(istate(k));
        if (istate(k) < 0e0) x(j) = bl(j);
        if (istate(k) > 0e0) x(j) = bu(j);
    }

    /* 
     * In a warm start (key != 0) initialize the free variables to (bl+bu)/2. 
     * This is needed in case the initial qr results in free variables 
     * out-of-bounds and Steps 8-11 get executed the first time through.
     */
    for (int k = nbound; k < n; k++) {
        int kk = istate(k);
        x(kk) = (bu(kk)+bl(kk)) / 2e0;
    }

    /* Compute bnorm, the norm of the data vector b, for reference. */
    double bnorm = b.norm();

    double obj,
           worst;
    int    j,
           it;
    /* Initialization complete. Begin major loop (Loop A). */
    for (int loopA = 0; loopA < 3*n; loopA++) { // 15000
        act.col(mm1-1) = b - A*x;
        /* 
         * Step 2
         * Initialize the negative gradient vector w(*). Compute the residual 
         * vector b-a.x, the negative gradient vector w(*), and the current 
         * objective value obj = ||a.x - b||. The residual vector is stored in 
         * the mm+1'st column of act(*,*).
         */
        w   = A.transpose()*act.col(mm1-1);
        obj = act.col(mm1-1).norm();
        if (obj <= tol*bnorm || (loopA > 1 && nbound == 0)) {
            istate(n) = nbound;
            w(0) = obj;
            return;
        }
        for (int k = nbound; k < n; k++) {
            j = istate(k);
            for (int i = 0; i < m; i++) {
                act(i, mm1) += a(i,j)*x(j);
            }
        }
        if (loopA == 1 && key != 0) {
            goto 6000
        }
        worst = 0e0; // 3000
        it = 1;
        for (int i = 0; i < nbound; i++) {
            int ks = std::abs(istate(i));
            double bad = w(ks) * std::copysign(1e0, istate(i));
            if (bad < worst) {
                it    = j;
                worst = bad;
                iact  = ks;
            }
        }
        if (worst >= 0e0) {
            istate(n) = nbound;
            w(0) = obj;
            return;
        }
        if (iact == jj) {
            w(jj) = 0e0;
            goto 3000
        }
        if (istate(it) > 0e0) bound = bu(iact);
        if (istate(it) < 0e0) bound = bl(iact);
        act.col(mm1) = act.col(mm1) + bound*a.col(iact);
        ifrom5 = istate(it);
        istate(it) = istate(nbound);
        --nbound;
        ++nact;
        istate(nbound+1) = iact;
        if (mm < nact) {
            std::cerr<<"\nToo many free variables in BVLS.";
            return 1;
        }
        for (int i = 0; i < m; i++) { // 6000
            act(i, mm1+1) = act(i, mm1);
            for (int k = nbound; k < n ; k++) {
                j = istate(k);
                act(i, nact+1-k*nbound) = a(i,j);
            }
        }
        CALL qr(m, nact, act, act(1,mm1+1), zz, resq);
        if (resq < 0e0 
            || (ifrom5 > 0 && zz(nact) > bu(iact))
            || (ifrom5 < 0e0 && zz(nact) < bl(iact)) ) {
            ++nbound;
            istate(nbound) = istate(nbound) * std::copysign(1e0, x(iact)-bu(iact));
            --nact;
            act.col(mm1) = act.col(mm1) - x(iact)*a.col(iact);
            ifrom5  = 0;
            w(iact) = 0e0;
            goto 3000
        }
        if (ifrom5) jj = 0;
        ifrom = 0
        for (int k = 0; k < nact; k++) {
            j = istate(k+nbound);
            if (zz(nact+1-k) < bl(j) || zz(nact+1-k) > bu(j)) {
                goto 8000
            }
        }
        for (int k = 0; k < nact; k++) {
            j = istate(k+nbound);
            x(j) = zz(nact+1-k);
        }
        goto 15000
        double alpha = 2e0, // 8000
                alf  = alpha;
        for (int k = k1; k < nact; k++) {
            j = istate(k+nbound);
            if (zz(nact+1-k) > bu(j)) alf=(bu(j)-x(j))/(zz(nact+1-k)-x(j));
            if (zz(nact+1-k) < bl(j)) alf=(bl(j)-x(j))/(zz(nact+1-k)-x(j));
            if (alf < alpha) {
                alpha = alf;
                jj = j;
                double sj = std::copysign(1e0, zz(nact+1-k)-bl(j));
            }
        }
        for (int k= 0 ; k < nact; k++) {
            j = istate(k+nbound);
            x(j) += alpha*(zz(nact+1-k)-x(j));
        }
        nold = nbound;
        for (int k = 0; k < nact; k++) {
            j = istate(k+nold)
            if (((bu(j)-x(j)) <= 0e0) || (j == jj && sj > 0e0)) {
                x(j) = bu(j);
                istate(k+noldb) = istate(nbound+1);
                istate(nbound+1) = j;
                ++nbound;
                act.col(mm1) = act.col(mm1) - bu(j)*a.col(j);
            } else if (((x(j)-bl(j) <= 0e0) || (j == jj && sj < 0e0)) {
                x(j) = bl(j);
                istate(k+noldb) = istate(nbound+1);
                istate(nbound+1) = -j;
                ++nbound;
                act.col(mm1) = act.col(mm1) - bl(j)*a.col(j);
            }
        }
        nact = n - nbound;
        if (nact > 0) goto 6000

    } // 15000

}

int main()
{
    Eigen::MatrixXd A = Eigen::MatrixXd(4,2);
    A << 0.0372, 0.2869, 0.6861, 0.7071, 0.6233, 0.6245, 0.6344, 0.6170;
    std::cout<<"\nMatrix A:\n"<<A<<"\n";
    Eigen::VectorXd y = Eigen::VectorXd(4);
    y << 0.8587, 0.1781, 0.0747, 0.8405;
    std::cout<<"\nVector y:\n"<<y<<"\n";
    auto x = fnnls(A, y);
    std::cout<<"\nResult is:\n"<<x<<"\n";

    std::cout<<"\n";
    return 0;
}
