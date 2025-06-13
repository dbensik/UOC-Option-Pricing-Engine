#ifndef GSL_TRIDIAGONAL_SOLVER_HPP
#define GSL_TRIDIAGONAL_SOLVER_HPP

#include <vector>

class GslTridiagonalSolver {
public:
    std::vector<double> solve(
        const std::vector<double>& a,
        const std::vector<double>& b,
        const std::vector<double>& c,
        const std::vector<double>& d
    );
};

#endif
