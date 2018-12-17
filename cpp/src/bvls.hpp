#include <set>
#include <cmath>
// Eigen headers
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/QR"
#include "eigen3/Eigen/Dense"

// http://digitalassets.lib.berkeley.edu/sdtr/ucb/text/394.pdf
bool
kuhn_tucker_test(const Eigen::VectorXd& w, const std::set<int>& L, 
    const std::set<int>& U, const std::set<int>& F)
{
    int n = w.rows();
    // check if F = {0, 1, ..., n-1}
    if ( F.size() == n ) return true;

    // check if w(j) <= 0 for all j in L
    for (auto it = L.cbegin(); it != L.cend(); ++it) {
        if ( w(*it) > 0e0 ) return false;
    }
    
    // check if w(j) >= 0 for all j in U
    for (auto it = U.cbegin(); it != U.cend(); ++it) {
        if ( w(*it) < 0e0 ) return false;
    }

    return true;
}

double
argmax(const Eigen::VectorXd& w, const std::set<int>& L,              
    const std::set<int>& U, int& idx)
{
    // 4. Find t* = argmax(s_t*w_t), t in L(union)U where st=1 if t E L and 
    // st=-1 if t E U. 
    // If t* is not unique, break ties arbitrarily.
    // If max index is in L, then the index is returned with opposite sign.
    
    double t = std::numeric_limits<double>::min();
    
    for (auto it = L.cbegin(); it != L.cend(); ++it) {
        if ( w(*it) > t ) {
            idx = *it * -1;
            t   = w(*it);
        }
    }
    for (auto it = U.cbegin(); it != U.cend(); ++it) {
        if ( -w(*it) > t ) {
            idx = *it;
            t   = w(*it);
        }
    }
    return t;
}

void
bvls(Eigen::MatrixXd& A, Eigen::VectorXd& b, Eigen::VectorXd& x_apr, int n, int m)
{
    std::set<int> L, // indexes of vars at their lower bound
                  U, // indexes of vars at their upper bound
                  F; // indexes of free vars
    
    auto x_bvls = x_apr;
    for (int i = 0; i < n; i++) {
        L.insert(i);
        x_bvls(i) = x_apr(i);
    }

    // 2. Compute w=AT(b-Ax), the negative of the gradient of the squared objective.
    auto w = A.transpose() * (b - A * x_bvls);

    // 3. Kuhn-Tucker test for convergence.
    if ( kuhn_tucker_test(w, L, U, F) ) goto 12

    // 4. Find t* = argmax(s_t*w_t), t in L(union)U where st=1 if t E L and 
    // st=-1 if t E U. 
    // If t* is not unique, break ties arbitrarily.
    int idx = 0;
    double t = argmax(w, L, U, idx);

    // 5. Move t* to the set F.
    if (idx < 0) {// idx is in L set
        F.insert(-idx);
        L.erase(-idx);
    } else {
        F.insert(idx);
        U.erase(idx);
    }
}
