// CustomerTridiagonalSolver.cpp
#include "CustomerTridiagonalSolver.hpp"
#include <stdexcept>

std::vector<double> CustomerTridiagonalSolver::solve(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d) {
    int n = d.size();
    std::vector<double> cp(n), dp(n), x(n);

    if (b[0] == 0.0) throw std::runtime_error("Division by zero in Thomas algorithm");
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];

    for (int i = 1; i < n; ++i) {
        double m = b[i] - a[i - 1] * cp[i - 1];
        if (m == 0.0) throw std::runtime_error("Division by zero in Thomas algorithm");
        cp[i] = c[i] / m;
        dp[i] = (d[i] - a[i - 1] * dp[i - 1]) / m;
    }

    x[n - 1] = dp[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }

    return x;
}
