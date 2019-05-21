
#ifndef uocOption_HPP
#define uocOption_HPP

//#include <boost/math/special_functions/gamma.hpp>  // for uigf
#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <time.h>
#include <gsl/gsl_multimin.h>  // not sure if this is needed
#include <vector>
#include <stdlib.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>

using namespace std;

class uocOption {
public:  //Public Variables

	double S;  // Initial stock price
	double x;
	double K;  // Strike price
	double LogK;
	double B;  // Upper barrier for the up-and-out call
	double LogB;
	double r;  // Risk-free rate
	double q;  // dividend rate on the stock
	double T;  // Expiration (years)
	double t; // "current" time
	double tau; // Time to maturity
	double sigma;  // volatility of the stock
	double nu;  //  Given   ????????
	double theta;  //  Given   ????????
	double Y;  //  Given   ????????
	double Smin;  // Need to choose this value for mesh size
	double LogSmin;
	int numPrices;
	int numTimes;
	double deltaTau;  // time step
	double deltaS;  // stock price step
	double deltax;
	double lambdaN;
	double lambdaP;
	double Bl;
	double Bu;
	double logSpot;
	int numRow, numCol;
	int arraySize;

public:// Public functions
	uocOption();							// Default call option
	
	uocOption(const double K, const double B, const double sigma, 
		const double nu, const double theta, const double Y);

	uocOption(const double sigma, const double nu, const double theta, const double Y);
	
	uocOption(const double strike, const double barrier, const vector<double>& params);
	
	uocOption(const uocOption& option2);	// Copy constructor
	
	virtual ~uocOption();					// Virtual destructor
	
	uocOption& operator = (const uocOption& option2);	// Assignmnet Operator
	
	double R(int i, int j, double **array, 
		const vector<double>& g1LambdaN, const vector<double>& g1LambdaP,
		const vector<double>& g2LambdaN, const vector<double>& g2LambdaP);

	double sigma2(double epsilon);
	
	double omega(double epsilon);
	
	double g1(double Y, double x);
	
	double g2(double Y, double x);
	
	double l();
	
	double d(int i);
	
	double u();

	void init();	// Initialise all default values
	
	void copy(const uocOption& o2);
	
	void dumpPrint();
	
	void printGrid(double **array) const;
	
	double price();

	double pricePrint();
	
	double payoff(double S) const;
	
	void writeCSV(double **array) const;
};

#endif