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
    std::cout<<"\n--M.size()="<<M.size()<<", "<<M[M.size()-1]+1;
    for (auto i : P) std::cout<<"\n\t(P)"<<i;
    for (auto i : R) std::cout<<"\n\t(R)"<<i;
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

void
fnnls(Eigen::MatrixXd& A, Eigen::VectorXd& y, double tol = -1e0)
{
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    // initialize variables
    if (tol < 0e0) {
        tol = 0e0/*10e0 * matlab_eps<double>(1e0) * XtX.col(0).norm() * XtX.rows()*/;

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
        std::cout<<"\n++ pushing "<<j<<" to P";
        for (auto i : R) std::cout<<"\n\t(R): "<<i;
        std::cout<<"\n++ removing "<<j<<" from R";
        R.erase(std::remove(R.begin(), R.end(), j), R.end());
        for (auto i : R) std::cout<<"\n\t(R): "<<i;
        std::sort(P.begin(), P.end());
        std::sort(R.begin(), R.end());
        // MatrixXd Ap = MatrixXd(m, P.size());
        MatrixXd Ap = make_Ap(A, P);
        VectorXd sp = VectorXd(P.size());
        sp = Ap.colPivHouseholderQr().solve(y);
        VectorXd s  = VectorXd::Zero(n);
        s = make_s(sp, P, R);

        // Inner Loop
        int m;
        double alpha {std::numeric_limits<double>::max()},
               atmp;
        while ( sp.minCoeff(&m) <= 0e0 ) {
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
                    std::cout<<"\n++ pushing "<<j<<" to R";
                    V.push_back(i);
                }
            }
            for (auto i : V) {
                std::cout<<"\n++ removing "<<i<<" from P";
                std::remove(P.begin(), P.end(), i);
            }
            P.erase(P.begin()+V.size(), P.end());
            std::sort(P.begin(), P.end());
            std::sort(R.begin(), R.end());
            Ap = make_Ap(A, P);
            sp = Ap.colPivHouseholderQr().solve(y);
            s = make_s(sp, P, R);
        }
        x = s;
        w = A.transpose()*(y-A*x);
    }

    return;
}

int main()
{
    Eigen::MatrixXd A = Eigen::MatrixXd(4,2);
    A << 0.0372, 0.2869, 0.6861, 0.7071, 0.6233, 0.6245, 0.6344, 0.6170;
    std::cout<<"\nMatrix A:\n"<<A<<"\n";
    Eigen::VectorXd y = Eigen::VectorXd(4);
    y << 0.8587, 0.1781, 0.0747, 0.8405;
    std::cout<<"\nVector y:\n"<<y<<"\n";
    fnnls(A, y);

    std::cout<<"\n";
    return 0;
}
