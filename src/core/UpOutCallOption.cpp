
#ifndef UpOutCallOption_cpp
#define UpOutCallOption_cpp

#include "UpOutCallOption.hpp"
#include "pricing/TridiagonalSolver.hpp"
#include "utils/MathUtils.hpp"
#include <tuple>

using std::vector;
using std::cout;
using std::endl;

UpOutCallOption::UpOutCallOption() {
    init();
}

UpOutCallOption::UpOutCallOption(double K, double B, double sigma, double nu, double theta, double Y)
    : K(K), B(B), sigma(sigma), nu(nu), theta(theta), Y(Y) {
    init();
}

UpOutCallOption::UpOutCallOption(double sigma, double nu, double theta, double Y)
    : sigma(sigma), nu(nu), theta(theta), Y(Y) {
    init();
}

UpOutCallOption::UpOutCallOption(double strike, double barrier, const vector<double>& params)
    : K(strike), B(barrier), sigma(params[0]), nu(params[1]), theta(params[2]), Y(params[3]) {
    init();
}

UpOutCallOption::UpOutCallOption(const UpOutCallOption& option2) {
    copy(option2);
}

UpOutCallOption::~UpOutCallOption() {}

UpOutCallOption& UpOutCallOption::operator=(const UpOutCallOption& option2) {
    if (this != &option2) {
        copy(option2);
    }
    return *this;
}

void UpOutCallOption::init() {
    S = 100; x = log(S);
    LogK = log(K);
    LogB = log(B);
    r = 0.0025; q = 0.015;
    T = 0.5; t = 0; tau = T - t;
    Smin = 5; LogSmin = log(Smin);
    numPrices = 100; numTimes = 100;
    deltaTau = tau / (numTimes - 1);
    deltaS = (B - Smin) / (numPrices - 1);
    deltax = (LogB - LogSmin) / (numPrices - 1);

    lambdaN = MathUtils::lambdaN(theta, sigma, nu);
    lambdaP = MathUtils::lambdaP(theta, sigma, nu);

    Bl = ((MathUtils::sigma2(sigma, deltax, Y) * deltaTau) / (2 * deltax * deltax)) - 
         (r - q + MathUtils::omega(sigma, theta, nu, deltax, Y) - 0.5 * MathUtils::sigma2(sigma, deltax, Y)) * (deltaTau / (2 * deltax));
    Bu = ((MathUtils::sigma2(sigma, deltax, Y) * deltaTau) / (2 * deltax * deltax)) + 
         (r - q + MathUtils::omega(sigma, theta, nu, deltax, Y) - 0.5 * MathUtils::sigma2(sigma, deltax, Y)) * (deltaTau / (2 * deltax));
}

void UpOutCallOption::copy(const UpOutCallOption& o2) {
    S = o2.S; x = o2.x; K = o2.K; LogK = o2.LogK;
    B = o2.B; LogB = o2.LogB; r = o2.r; q = o2.q;
    T = o2.T; t = o2.t; tau = o2.tau;
    sigma = o2.sigma; nu = o2.nu; theta = o2.theta; Y = o2.Y;
    Smin = o2.Smin; LogSmin = o2.LogSmin;
    numPrices = o2.numPrices; numTimes = o2.numTimes;
    deltaTau = o2.deltaTau; deltaS = o2.deltaS; deltax = o2.deltax;
    lambdaN = o2.lambdaN; lambdaP = o2.lambdaP;
    Bl = o2.Bl; Bu = o2.Bu;
}

void UpOutCallOption::dumpPrint() {
    cout << "\nStock price: $" << S << ", Strike: $" << K << ", Barrier: $" << B << endl;
    cout << "Rates -> r: " << r << ", q: " << q << ", Vol: " << sigma << endl;
    cout << "Params -> nu: " << nu << ", theta: " << theta << ", Y: " << Y << endl;
    cout << "Time -> T: " << T << ", t: " << t << ", tau: " << tau << endl;
    cout << "Mesh -> Prices: " << numPrices << ", Times: " << numTimes << endl;
    cout << "Steps -> deltaTau: " << deltaTau << ", deltaS: " << deltaS << ", deltax: " << deltax << endl;
    cout << "Lambdas -> N: " << lambdaN << ", P: " << lambdaP << endl;
    cout << "Coefficients -> Bl: " << Bl << ", Bu: " << Bu << endl;
}

// Updated price method with TridiagonalSolver
double UpOutCallOption::price() {
    vector<vector<double>> grid(numTimes, vector<double>(numPrices, 0.0));

    for (int j = 0; j < numPrices; ++j) {
        double logS = LogSmin + j * deltax;
        double S = exp(logS);
        grid[numTimes - 1][j] = std::max(S - K, 0.0);
    }

    for (int i = 0; i < numTimes; ++i) {
        grid[i][0] = 0.0;
        grid[i][numPrices - 1] = 0.0;
    }

    TridiagonalSolver solver;
    vector<double> a(numPrices - 2), b(numPrices - 2), c(numPrices - 2), rhs(numPrices - 2);

    for (int i = numTimes - 2; i >= 0; --i) {
        for (int j = 1; j < numPrices - 1; ++j) {
            a[j - 1] = Bl;
            b[j - 1] = 1 - 2 * Bl;
            c[j - 1] = Bu;
            rhs[j - 1] = grid[i + 1][j];
        }

        vector<double> solution = solver.solve(a, b, c, rhs);
        for (int j = 1; j < numPrices - 1; ++j) {
            grid[i][j] = solution[j - 1];
        }
    }

    int j = static_cast<int>((x - LogSmin) / deltax);
    double w = (x - (LogSmin + j * deltax)) / deltax;
    return (1 - w) * grid[0][j] + w * grid[0][j + 1];
}
