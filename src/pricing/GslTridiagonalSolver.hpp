// GslTridiagonalSolver.hpp
#ifndef GSL_TRIDIAGONAL_SOLVER_HPP
#define GSL_TRIDIAGONAL_SOLVER_HPP

#include "TridiagonalSolverBase.hpp"
#include <gsl/gsl_linalg.h>

class GslTridiagonalSolver : public TridiagonalSolverBase {
public:
    std::vector<double> solve(
        const std::vector<double>& a,
        const std::vector<double>& b,
        const std::vector<double>& c,
        const std::vector<double>& rhs
    ) const override;
};

#endif // GSL_TRIDIAGONAL_SOLVER_HPP
