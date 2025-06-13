// CustomTridiagonalSolver.cpp
#include "CustomTridiagonalSolver.hpp"
#include <vector>

std::vector<double> CustomTridiagonalSolver::solve(
    const std::vector<double>& a,
    const std::vector<double>& b,
    const std::vector<double>& c,
    const std::vector<double>& rhs
) const {
    int n = rhs.size();
    std::vector<double> x(n), c_prime(n), d_prime(n);

    c_prime[0] = c[0] / b[0];
    d_prime[0] = rhs[0] / b[0];

    for (int i = 1; i < n; ++i) {
        double m = b[i] - a[i - 1] * c_prime[i - 1];
        c_prime[i] = c[i] / m;
        d_prime[i] = (rhs[i] - a[i - 1] * d_prime[i - 1]) / m;
    }

    x[n - 1] = d_prime[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }

    return x;
}
