// UpOutCallOption.hpp
#ifndef UPOUTCALLOPTION_HPP
#define UPOUTCALLOPTION_HPP

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
#include <gsl/gsl_multimin.h>
#include <vector>
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

class UpOutCallOption {
public:
    double S, x, K, LogK, B, LogB;
    double r, q, T, t, tau, sigma, nu, theta, Y;
    double Smin, LogSmin;
    int numPrices, numTimes;
    double deltaTau, deltaS, deltax;
    double lambdaN, lambdaP, Bl, Bu;
    double logSpot;
    int numRow, numCol, arraySize;

    UpOutCallOption();
    UpOutCallOption(double K, double B, double sigma, double nu, double theta, double Y);
    UpOutCallOption(double sigma, double nu, double theta, double Y);
    UpOutCallOption(double strike, double barrier, const vector<double>& params);
    UpOutCallOption(const UpOutCallOption& option2);
    virtual ~UpOutCallOption();

    UpOutCallOption& operator = (const UpOutCallOption& option2);

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

    void init();
    void copy(const UpOutCallOption& o2);
    void dumpPrint();
    void printGrid(double **array) const;
    double price();
    double pricePrint();
    double payoff(double S) const;
    void writeCSV(double **array) const;
};

#endif
