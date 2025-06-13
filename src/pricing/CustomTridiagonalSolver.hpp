// CustomTridiagonalSolver.hpp
#ifndef CUSTOM_TRIDIAGONAL_SOLVER_HPP
#define CUSTOM_TRIDIAGONAL_SOLVER_HPP

#include "TridiagonalSolverBase.hpp"

class CustomTridiagonalSolver : public TridiagonalSolverBase {
public:
    std::vector<double> solve(
        const std::vector<double>& a,
        const std::vector<double>& b,
        const std::vector<double>& c,
        const std::vector<double>& rhs
    ) const override;
};

#endif // CUSTOM_TRIDIAGONAL_SOLVER_HPP
