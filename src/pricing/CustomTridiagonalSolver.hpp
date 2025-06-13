// CustomerTridiagonalSolver.hpp
#ifndef CUSTOMERTRIDIAGONALSOLVER_HPP
#define CUSTOMERTRIDIAGONALSOLVER_HPP

#include "TridiagonalSolverBase.hpp"
#include <vector>

class CustomerTridiagonalSolver : public TridiagonalSolverBase {
public:
    std::vector<double> solve(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d) override;
};

#endif // CUSTOMERTRIDIAGONALSOLVER_HPP
