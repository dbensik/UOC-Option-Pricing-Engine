// MathUtils.cpp
#include "MathUtils.hpp"
#include <cmath>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_expint.h>
#include <stdexcept>

namespace MathUtils {
	double lambdaN(const double theta, const double sigma, const double nu) {
		return std::sqrt((theta * theta) / (sigma * sigma * sigma * sigma) + (2.0 / (sigma * sigma * nu))) + theta / (sigma * sigma);
	}

	double lambdaP(const double theta, const double sigma, const double nu) {
		return std::sqrt((theta * theta) / (sigma * sigma * sigma * sigma) + (2.0 / (sigma * sigma * nu))) - theta / (sigma * sigma);
	}

	// Computes the incomplete gamma function integral based on Y and x.
	// Valid for Y in [0, 1). Throws an error otherwise.
	double g1(const double Y, const double x) {
		if (x < 0.0)
			throw std::domain_error("x must be non-negative in g1");

		if (Y == 0.0)
			return std::exp(-x);
		else if (Y > 0.0 && Y < 1.0)
			return gsl_sf_gamma_inc(1.0 - Y, x);
		else
			throw std::runtime_error("Invalid value of Y in g1: expected 0 <= Y < 1");
	}

	double g2(const double Y, const double x) {
		if (x <= 0.0)
			throw std::domain_error("x must be positive in g2");

		if (Y == 0.0)
			return gsl_sf_expint_E1(x);
		else if (Y > 0.0 && Y < 1.0)
			return (std::exp(-x) * std::pow(x, -Y) - gsl_sf_gamma_inc(1.0 - Y, x)) / Y;
		else
			throw std::runtime_error("Invalid value of Y in g2: expected 0 <= Y < 1");
	}

	// Computes sigma squared using lambdaP, lambdaN, deltax, nu, and Y.
	// This function can be broken into components for clarity.
	double sigma2(const double lambdaN, const double lambdaP, const double deltax, const double nu, const double Y) {
		const double g1_Y_0 = g1(Y, 0.0);

		double powerP = std::pow(lambdaP, Y - 2.0);
		double scaleP = std::pow(lambdaP * deltax, 1.0 - Y) * std::exp(-lambdaP * deltax);
		double gammaDiffP = g1_Y_0 - g1(Y, lambdaP * deltax);
		double term1 = (1.0 / nu) * powerP * (-scaleP + (1.0 - Y) * gammaDiffP);

		double powerN = std::pow(lambdaN, Y - 2.0);
		double scaleN = std::pow(lambdaN * deltax, 1.0 - Y) * std::exp(-lambdaN * deltax);
		double gammaDiffN = g1_Y_0 - g1(Y, lambdaN * deltax);
		double term2 = (1.0 / nu) * powerN * (-scaleN + (1.0 - Y) * gammaDiffN);

		return term1 + term2;
	}

	// Computes omega using jump intensities and an auxiliary gamma-based function.
	// Expected range for Y: 0 <= Y < 1
	double omega(const double lambdaN, const double lambdaP, const double nu, const double deltax, const double Y) {
		if (Y < 0.0 || Y >= 1.0)
			throw std::runtime_error("Invalid value of Y in omega: expected 0 <= Y < 1");

		double term1 = (std::pow(lambdaP, Y) / nu) * g2(Y, lambdaP * deltax);
		double term2 = (std::pow(lambdaP - 1.0, Y) / nu) * g2(Y, (lambdaP - 1.0) * deltax);
		double term3 = (std::pow(lambdaN, Y) / nu) * g2(Y, lambdaN * deltax);
		double term4 = (std::pow(lambdaN + 1.0, Y) / nu) * g2(Y, (lambdaN + 1.0) * deltax);

		return term1 - term2 + term3 - term4;
	}

} // namespace MathUtils
