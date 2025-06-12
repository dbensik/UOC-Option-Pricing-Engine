// include/Parameters.hpp
#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

struct OptionParams {
    double S = 100.0;         // Initial stock price
    double K = 100.0;         // Strike price
    double B = 120.0;         // Barrier level
    double r = 0.0025;        // Risk-free rate
    double q = 0.015;         // Dividend yield
    double T = 0.5;           // Time to expiration (in years)
    double sigma = 0.25;      // Volatility
    double nu = 0.35;         // Jump size parameter
    double theta = -0.4;      // Asymmetry parameter
    double Y = 0.4;           // Activity parameter
};

#endif // PARAMETERS_HPP
