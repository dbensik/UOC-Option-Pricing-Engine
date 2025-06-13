// MathUtils.cpp
#include "MathUtils.hpp"
#include <cmath>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_expint.h>

namespace MathUtils {
	double lambdaN(double theta, double sigma, double nu) {
		return std::sqrt((theta * theta) / (sigma * sigma * sigma * sigma) + (2.0 / (sigma * sigma * nu))) + theta / (sigma * sigma);
	}

	double lambdaP(double theta, double sigma, double nu) {
		return std::sqrt((theta * theta) / (sigma * sigma * sigma * sigma) + (2.0 / (sigma * sigma * nu))) - theta / (sigma * sigma);
	}

	double g1(double Y, double x) {
		if (Y == 0)
			return std::exp(-x);
		else if (Y > 0 && Y < 1)
			return gsl_sf_gamma_inc(1 - Y, x);
		else {
			throw std::runtime_error("Invalid value of Y in g1");
		}
	}

	double g2(double Y, double x) {
		if (Y == 0)
			return gsl_sf_expint_E1(x);
		else if (Y > 0 && Y < 1)
			return (std::exp(-x) * std::pow(x, -Y) - gsl_sf_gamma_inc(1 - Y, x)) / Y;
		else {
			throw std::runtime_error("Invalid value of Y in g2");
		}
	}

	double sigma2(double lambdaN, double lambdaP, double deltax, double nu, double Y) {
		double term1 = (1 / nu) * std::pow(lambdaP, Y - 2) * (-std::pow(lambdaP * deltax, 1 - Y) * std::exp(-lambdaP * deltax) + (1 - Y) * (g1(Y, 0) - g1(Y, lambdaP * deltax)));
		double term2 = (1 / nu) * std::pow(lambdaN, Y - 2) * (-std::pow(lambdaN * deltax, 1 - Y) * std::exp(-lambdaN * deltax) + (1 - Y) * (g1(Y, 0) - g1(Y, lambdaN * deltax)));
		return term1 + term2;
	}

	double omega(double lambdaN, double lambdaP, double nu, double deltax, double Y) {
		double term1 = (std::pow(lambdaP, Y) / nu) * g2(Y, lambdaP * deltax);
		double term2 = (std::pow(lambdaP - 1, Y) / nu) * g2(Y, (lambdaP - 1) * deltax);
		double term3 = (std::pow(lambdaN, Y) / nu) * g2(Y, lambdaN * deltax);
		double term4 = (std::pow(lambdaN + 1, Y) / nu) * g2(Y, (lambdaN + 1) * deltax);
		return term1 - term2 + term3 - term4;
	}
} // namespace MathUtils
