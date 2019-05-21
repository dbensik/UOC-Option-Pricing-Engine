#ifndef optimizer_HPP
#define optimizer_HPP

#include "uocOption.hpp"

//#include <omp.h>

using namespace std;

void testTridiagonal();

double vectorMin(const vector<double>& v);

int vectorMinIndex(const vector<double>& v);

double constrainSigma(double& sigma);

double constrainNu(double& nu);

double constrainTheta(double& theta);

double constrainY(double& Y);

void testConstraints();

void gridSearch(vector<double>& myInitialTheta,
	vector<double>& myOptimalTheta,
	const vector<double>& myStrikes,
	const vector<double>& myBarriers,
	const vector< vector<double> >& marketPrices,
	vector<vector<double> >& modelPrices);

void getPriceGrid(const vector<double>& Theta,
					const vector<double>& myStrikes,
					const vector<double>& myBarriers,
					vector<vector<double> >& modelPrices);

void ssre(const vector<vector<double> >& modelPrices,
			const vector<vector<double> >& marketPrices,
			vector<double>& diffVector);

double ssre(const vector<vector<double> >& modelPrices,
	const vector<vector<double> >& marketPrices);

void ssre1(const vector<vector<double> >& modelPrices,
			const vector<vector<double> >& marketPrices,
			vector<double>& diffVector);

double ssre1(const vector<vector<double> >& modelPrices,
	const vector<vector<double> >& marketPrices);

void ssre2(const vector<vector<double> >& modelPrices,
			const vector<vector<double> >& marketPrices,
			vector<double>& diffVector);

double ssre2(const vector<vector<double> >& modelPrices,
	const vector<vector<double> >& marketPrices);

double ssd(const vector<vector<double> >& modelPrices,
	const vector<vector<double> >& marketPrices);

double round(double val, int precision);

void printTheta(const vector<double>& Theta);

void printInitialTheta(const vector<double>& myInitialTheta);

void printTestTheta(const vector<double>& myTestTheta);

void printSimplexTheta(const vector<double>& mySimplexTheta);

void printOptimalTheta(const vector<double>& myOptimalTheta);

void printConstrainedTheta(const vector<double>& myInitialTheta);

void printVector(const vector<double>& v);

void print2DVector(const vector< vector<double> >& vec2D);

double objective(vector<double>& mySimplexTheta);

double myFunc(const gsl_vector *v, void *params);

double myFunc2(const gsl_vector *v, void *params);

int nelderMeadExample();

int nelderMead1(const vector<double>& Theta);

int nelderMead2(const vector<double>& Theta);

// class nelderMead{
// public:
// 	//calculate market prices once here
// 	vector< vector<double> > marketPrices(3, vector<double>(4));
// 	getPriceGrid(myInitialTheta,myStrikes,myBarriers,marketPrices);
	
// 	double myFunc2(const gsl_vector *v, void *params, 
// 		const vector<vector<double> >& marketPrices){
// 		double *dp = (double *)params;

// 	    vector<double> mySimplexTheta(4);
// 	    vector<double> myInitialTheta(4);

// 	    mySimplexTheta[0] = gsl_vector_get(v, 0);
// 	    mySimplexTheta[1] = gsl_vector_get(v, 1);
// 	    mySimplexTheta[2] = gsl_vector_get(v, 2);
// 	    mySimplexTheta[3] = gsl_vector_get(v, 3);

// 	    myInitialTheta[0] = dp[0];
// 	    myInitialTheta[1] = dp[1];
// 	    myInitialTheta[2] = dp[2];
// 	    myInitialTheta[3] = dp[3];

// 	    vector<double> myStrikes(4);
// 	    vector<double> myBarriers(3);
// 	    myStrikes[0] = 100.0; myStrikes[1] = 105.0; myStrikes[2] = 110.0; myStrikes[3] = 120.0;
// 	    myBarriers[0] = 125.0, myBarriers[1] = 135.0, myBarriers[2] = 145.0;
	    
// 	    vector< vector<double> > marketPrices(3, vector<double>(4));
// 	    vector< vector<double> > modelPrices(3, vector<double>(4));
// 	    getPriceGrid(mySimplexTheta,myStrikes,myBarriers,modelPrices); //model
// 	    getPriceGrid(myInitialTheta,myStrikes,myBarriers,marketPrices); // market
	    
// 	    double objDiff = 0.0;
// 	    objDiff = ssre(modelPrices,marketPrices);

// 	    return objDiff;
// 	}

// 	int operator()(const vector<double>& Theta){

// 	}
// };

#endif