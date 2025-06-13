#include "core/UpOutCallOption.hpp"
#include "pricing/CustomerTridiagonalSolver.hpp"
#include "pricing/GslTridiagonalSolver.hpp"

int main() {
    CustomerTridiagonalSolver solver;
    UpOutCallOption option(100, 125, 0.3, 0.4, -0.2, 0.7, &solver);
    option.dumpPrint();
    double price = option.price();
    std::cout << "Option Price: " << price << std::endl;
    return 0;
}
