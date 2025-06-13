#ifndef CUSTOM_TRIDIAGONAL_SOLVER_HPP
#define CUSTOM_TRIDIAGONAL_SOLVER_HPP

#include <vector>

class CustomTridiagonalSolver {
public:
    std::vector<double> solve(
        const std::vector<double>& a,
        const std::vector<double>& b,
        const std::vector<double>& c,
        const std::vector<double>& d
    );
};

#endif // CUSTOM_TRIDIAGONAL_SOLVER_HPP
