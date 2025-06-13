// MathUtils.hpp
#ifndef MATH_UTILS_HPP
#define MATH_UTILS_HPP

namespace MathUtils {
	double lambdaN(double theta, double sigma, double nu);
	double lambdaP(double theta, double sigma, double nu);
	double sigma2(double sigma, double deltax, double Y);
	double omega(double sigma, double theta, double nu, double deltax, double Y);
}

#endif // MATH_UTILS_HPP
