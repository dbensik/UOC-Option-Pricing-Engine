// TridiagonalSolverBase.hpp
#ifndef TRIDIAGONAL_SOLVER_BASE_HPP
#define TRIDIAGONAL_SOLVER_BASE_HPP

#include <vector>

class TridiagonalSolverBase {
public:
    virtual ~TridiagonalSolverBase() {}
    virtual std::vector<double> solve(
        const std::vector<double>& a,
        const std::vector<double>& b,
        const std::vector<double>& c,
        const std::vector<double>& rhs
    ) const = 0;
};

#endif // TRIDIAGONAL_SOLVER_BASE_HPP
