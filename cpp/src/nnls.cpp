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

bool
__kuhn_tucker_test__(const Eigen::VectorXd& w, const std::vector<std::size_t>& L,
    const std::vector<std::size_t>& Ui, std::size_t& t) noexcept
{
    t = 0;
    double argmax = std::numeric_limits<double>::min();

    for (auto i = L.cbegin(); L.cend(); ++i) {
        if ( w(*it) > 0e0 ) {
            return false;
        } else {
            if ( w(*it) >= argmax ) t = *it;
        }
    }
    for (auto i = U.cbegin(); U.cend(); ++i) {
        if ( w(*it) < 0e0 ) {
            return false;
        } else {
            if ( -w(*it) >= argmax ) t = *it;
        }
    }
    return true;
}

auto
__data_vector_less_bound_vars__(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
    const std::vector<std::size_t>& F, const Eigen::VectorXd& x)
{
    std::size_t fs = F.size();
    std::size_t lu = A.cols() - fs;

    // data vector less the predictions of the bound variables
    Eigen::VectorXd bt = Eigen::VectorXd::Zero(lu);
    Eigen::VectorXd xt = Eigen::VectorXd::Zero(lu);
    Eigen::MatrixXd A_ = Eigen::Matrix(A.rows(), lu);
    std::size_t bt_row = 0,
                xt_row = 0,
                a__col = 0;

    // the matrix composed of those columns of A whose indices are in F
    Eigen::MatrixXd At = Eigen::MatrixXd(A.rows(), fs);
    std::size_t at_col = 0;

    for (std::size_t i = 0; i < A.cols(); i++) {
        // index i is in F
        if ( std::find(F.cbegin(), F.cend(), i) != F.cend() ) {
            At.col(at_col) = A.col(i);
            ++at_col;
        } else { // index i is in L(union)U
           A_.col(a__col) = A.col(i); 
        }
    }
}

auto
bvls(const Eigen::MatrixXd& A, const Eigen::VectorXd& y)
{
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    int m = A.rows();
    int n = A.cols();
    assert( y.rows() == m );
    VectorXd x = VectorXd::Zero(n);
    VectorXd w = VectorXd::Zero(n);

    // STEP 1
    std::vector<std::size_t> F,
                             U,
                             L,
                             FS;
    F.reserve(n);                     // F is an empty set (of indexes)
    U.reserve(n);                     // U is an empty set (of indexes)
    std::iota(L.begin(), L.end(), 0); // L contains all indexes
    FS = L;                           // FS contains all indexes (ordered)
    
    // STEP 2
    w = A.transpose()*(y-A*x);

    // STEP 3
    std::size_t t;
    if ( (F == FS) || __kuhn_tucker_test__(w, L, U, t) ) {
        return x;
    }
    F.push_back(t);

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
