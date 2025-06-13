
#ifndef TridiagonalSolver_hpp
#define TridiagonalSolver_hpp

#include <vector>

class TridiagonalSolver {
public:
    std::vector<double> solve(
        const std::vector<double>& a,  // lower diagonal (n-1)
        const std::vector<double>& b,  // main diagonal (n)
        const std::vector<double>& c,  // upper diagonal (n-1)
        const std::vector<double>& d   // right-hand side (n)
    );
};

#endif // TridiagonalSolver_hpp
