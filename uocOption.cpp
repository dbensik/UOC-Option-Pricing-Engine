#ifndef uocOption_cpp
#define uocOption_cpp

#include "uocOption.hpp"

/////////////////////////////////////////////////////////////////////////////////////
void uocOption::init() { // Default call option
	cout << "Inside uocOption::init()" << endl;
	//init();
}

uocOption::uocOption(){	// Initialise all default values
	//cout << "Inside uocOption::uocOption()" << endl;
	S = 1950;  // Initial stock price
	x = log(S); // Initial log price
	K = 2050;  // Strike price
	LogK = log(K);
	B = 2300;  // Upper barrier for the up-and-out call
	LogB = log(B);
	r = 0.0025;  // Risk-free rate
	q = 0.015;  // dividend rate on the stock
	T = 0.5;  // Expiration (years)
	t = 0; // "current" time
	tau = T-t; // Time to maturity
	sigma = 0.25;  // volatility of the stock
	nu = 0.35;  //  Given   ????????
	theta = -0.4;  //  Given   ????????
	Y = 0.4;  //  Given   ????????
	Smin = 1400;  // 1400
	LogSmin = log(Smin);
	numPrices = 400; // 400
	numTimes = 200; // 200
	deltaTau = (T - t)/(numTimes - 1);  // time step
	deltaS = (B - Smin)/(numPrices - 1) ;  // stock price step
	deltax = (LogB - LogSmin)/(numPrices - 1);

	lambdaN = sqrt((pow(theta,2)/pow(sigma,4))
		+(2/(sigma * sigma * nu))) + theta/(sigma * sigma);
	lambdaP = sqrt((pow(theta,2)/pow(sigma,4))
		+(2/(sigma * sigma * nu))) - theta/(sigma * sigma);

	Bl = ((uocOption::sigma2(deltax) * deltaTau) / (2 * deltax * deltax)) - 
		(r - q + uocOption::omega(deltax) - 0.5 * uocOption::sigma2(deltax))
		 * (deltaTau / (2 * deltax));
	
	Bu = ((uocOption::sigma2(deltax) * deltaTau) / (2 * deltax * deltax)) + 
		(r - q + uocOption::omega(deltax) - 0.5 * uocOption::sigma2(deltax))
		 * (deltaTau / (2 * deltax));
}

uocOption::uocOption(const double strike, const double barrier, const double sigmaParam, 
	const double nuParam, const double thetaParam, const double YParam){
    //cout << "uocOption Constructor with parameters" << endl;
    S = 100;  // Initial stock price
	x = log(S); // Initial log price
	K = strike;  // Strike price
	LogK = log(K);
	B = barrier;  // Upper barrier for the up-and-out call
	LogB = log(B);
	r = 0.0025;  // Risk-free rate
	q = 0.015;  // dividend rate on the stock
	T = 0.5;  // Expiration (years)
	t = 0; // "current" time
	tau = T-t; // Time to maturity
	sigma = sigmaParam;  // volatility of the stock
	nu = nuParam;  //  Given   ????????
	theta = thetaParam;  //  Given   ????????
	Y = YParam;  //  Given   ????????
	Smin = 5;  // Need to choose this value for mesh size
	LogSmin = log(Smin);
	numPrices = 100;
	numTimes = 100;
	deltaTau = (T - t)/(numTimes - 1);  // time step
	deltaS = (B - Smin)/(numPrices - 1) ;  // stock price step
	deltax = (LogB - LogSmin)/(numPrices - 1);

	lambdaN = sqrt((pow(theta,2)/pow(sigma,4))
		+(2/(sigma * sigma * nu))) + theta/(sigma * sigma);
	lambdaP = sqrt((pow(theta,2)/pow(sigma,4))
		+(2/(sigma * sigma * nu))) - theta/(sigma * sigma);

	Bl = ((uocOption::sigma2(deltax) * deltaTau) / (2 * deltax * deltax)) - 
		(r - q + uocOption::omega(deltax) - 0.5 * uocOption::sigma2(deltax))
		 * (deltaTau / (2 * deltax));
	
	Bu = ((uocOption::sigma2(deltax) * deltaTau) / (2 * deltax * deltax)) + 
		(r - q + uocOption::omega(deltax) - 0.5 * uocOption::sigma2(deltax))
		 * (deltaTau / (2 * deltax));  
}

uocOption::uocOption(const double sigmaParam, const double nuParam,
	const double thetaParam, const double YParam){
    //cout << "uocOption Constructor with parameters" << endl;
    S = 100;  // Initial stock price
	x = log(S); // Initial log price
	K = 100;  // Strike price
	LogK = log(K);
	B = 125;  // Upper barrier for the up-and-out call
	LogB = log(B);
	r = 0.0025;  // Risk-free rate
	q = 0.015;  // dividend rate on the stock
	T = 0.5;  // Expiration (years)
	t = 0; // "current" time
	tau = T-t; // Time to maturity
	sigma = sigmaParam;  // volatility of the stock
	nu = nuParam;  //  Given   ????????
	theta = thetaParam;  //  Given   ????????
	Y = YParam;  //  Given   ????????
	Smin = 5;  // Need to choose this value for mesh size
	LogSmin = log(Smin);
	numPrices = 100;
	numTimes = 100;
	deltaTau = (T - t)/(numTimes - 1);  // time step
	deltaS = (B - Smin)/(numPrices - 1) ;  // stock price step
	deltax = (LogB - LogSmin)/(numPrices - 1);

	lambdaN = sqrt((pow(theta,2)/pow(sigma,4))
		+(2/(sigma * sigma * nu))) + theta/(sigma * sigma);
	lambdaP = sqrt((pow(theta,2)/pow(sigma,4))
		+(2/(sigma * sigma * nu))) - theta/(sigma * sigma);

	Bl = ((uocOption::sigma2(deltax) * deltaTau) / (2 * deltax * deltax)) - 
		(r - q + uocOption::omega(deltax) - 0.5 * uocOption::sigma2(deltax))
		 * (deltaTau / (2 * deltax));
	
	Bu = ((uocOption::sigma2(deltax) * deltaTau) / (2 * deltax * deltax)) + 
		(r - q + uocOption::omega(deltax) - 0.5 * uocOption::sigma2(deltax))
		 * (deltaTau / (2 * deltax));  
}

uocOption::uocOption(const double strike, const double barrier, const vector<double>& params){
    //cout << "uocOption Constructor with parameter array" << endl;
    S = 100;  // Initial stock price
	x = log(S); // Initial log price
	K = strike;  // Strike price
	LogK = log(K);
	B = barrier;  // Upper barrier for the up-and-out call
	LogB = log(B);
	r = 0.0025;  // Risk-free rate
	q = 0.015;  // dividend rate on the stock
	T = 0.5;  // Expiration (years)
	t = 0; // "current" time
	tau = T-t; // Time to maturity
	sigma = params[0];  // volatility of the stock
	nu = params[1];  //  Given
	theta = params[2];  //  Give
	Y = params[3];  //  Given 
	Smin = 5;  // Need to choose this value for mesh size
	LogSmin = log(Smin);
	numPrices = 100;
	numTimes = 100;
	deltaTau = (T - t)/(numTimes - 1);  // time step
	deltaS = (B - Smin)/(numPrices - 1) ;  // stock price step
	deltax = (LogB - LogSmin)/(numPrices - 1);

	lambdaN = sqrt((pow(theta,2)/pow(sigma,4))
		+(2/(sigma * sigma * nu))) + theta/(sigma * sigma);
	lambdaP = sqrt((pow(theta,2)/pow(sigma,4))
		+(2/(sigma * sigma * nu))) - theta/(sigma * sigma);

	Bl = ((uocOption::sigma2(deltax) * deltaTau) / (2 * deltax * deltax)) - 
		(r - q + uocOption::omega(deltax) - 0.5 * uocOption::sigma2(deltax))
		 * (deltaTau / (2 * deltax));
	
	Bu = ((uocOption::sigma2(deltax) * deltaTau) / (2 * deltax * deltax)) + 
		(r - q + uocOption::omega(deltax) - 0.5 * uocOption::sigma2(deltax))
		 * (deltaTau / (2 * deltax));  
}

void uocOption::copy(const uocOption& o2){
	S = o2.S;  // Initial stock price
	x = o2.x;
	K = o2.K;  // Strike price
	LogK = o2.LogK;
	B = o2.B;  // Upper barrier for the up-and-out call
	LogB - o2.LogB;
	r = o2.r;  // Risk-free rate
	q = o2.q;  // dividend rate on the stock
	T = o2.T;  // Expiration (years)
	t = o2.t; // "current" time
	tau = o2.tau; // Time to maturity
	sigma = o2.sigma;  // volatility of the stock
	nu = o2.nu;  //  Given   ????????
	theta = o2.theta;  //  Given   ????????
	Y = o2.Y;  //  Given   ????????
	Smin = o2.Smin;  // Need to choose this value for mesh size
	LogSmin = o2.LogSmin;
	numPrices = o2.numPrices;
	numTimes = o2.numTimes;
	deltaTau = o2.deltaTau;  // time step
	deltaS = o2.deltaS;  // stock price step
	deltax = o2.deltax;

	lambdaN = o2.lambdaN;
	lambdaP = o2.lambdaP;
	Bl = o2.Bl;
	Bu = o2.Bu;	
}

uocOption::uocOption(const uocOption& o2){ // Copy constructor
	cout << "uocOption copy constructor" << endl;
	copy(o2);
}

uocOption::~uocOption() {
	//cout << "uocOption destructor" << endl;
	//cout << endl;// Virual Destructor
}

uocOption& uocOption::operator = (const uocOption& option2) {
	if (this == &option2) return *this;
	copy (option2);
	return *this;
}

void uocOption::dumpPrint() {
	cout << "\nStock price is $" << S << endl;
	cout << "Log Stock price is $" << x << endl;
	cout << "Strike price is $" << K << endl;
	cout << "Log Strike price is $" << LogK << endl;
	cout << "Knockout Barrier price is $" << B << endl;
	cout << "Log Knockout Barrier price is $" << LogB << endl;
	cout << "Minimum Stock Price is $" << Smin << endl;
	cout << "Minimum Log Stock Price is $" << LogSmin << endl;
	cout << "Risk-free rate is " << r*100 << "%" << endl;
	cout << "Dividend rate is " << q*100 <<"%" << endl;
	cout << "Expiration, T, (years) is " << T << endl;
	cout << "Current time, t, (years) is " << t << endl;
	cout << "Time to maturity, tau, (years) is " << tau << endl;
	cout << "sigma is " << sigma << endl;
	cout << "nu is " << nu << endl;
	cout << "theta is " << theta << endl;
	cout << "Y is " << Y << endl;
	cout << "There will be " << numTimes << " deltaTaus";
	cout << " beginning at expiration when tau = " << T-t;
	cout << " and ending at the current time when tau = " << T-0;
	cout << " with each deltaTau = " << deltaTau << endl;
	cout << "There will be " << numPrices << " price intervals";
	cout << " beginning at " << Smin;
	cout << " and ending at " << B << " with each step equal to ";
	cout << deltaS << endl;
	cout << "There will be " << numPrices << " log price intervals";
	cout << " beginning at " << LogSmin;
	cout << " and ending at " << LogB << " with each step equal to ";
	cout << deltax << endl;
	cout << "lambdaN is " << lambdaN << endl;
	cout << "lambdaP is " << lambdaP << endl;
	cout << "Bu is " << Bu << endl;//" at tau = " << tau << endl;
	cout << "Bl is " << Bl << endl;// " at tau = " << tau << endl;
	cout << "sigma2(deltax): " << sigma2(deltax) << endl;
	cout << "omega(deltax): " << omega(deltax) << endl;
	cout << "Grid is " << numPrices << " x " << numTimes << endl;
	cout << endl;
}

void uocOption::printGrid(double **array) const {
	cout << endl;
	for (int i = 0 ; i < numPrices; i++) {// loop over rows
        for (int j = 0; j < numTimes; j++) {// loop over cols
            cout << array[i][j] << " | ";
        }
        cout << endl;
    }
    cout << endl;
}

void uocOption::writeCSV(double **array) const {
	//gettimeofday(&startTime, NULL); 
	//FILE *surface;
	ofstream surface ("surface.csv"); 
    //surface = fopen("surface.csv", "w");

    if (surface.is_open())
	{
		for (int i = 0 ; i < numPrices; i++) // loop over rows
    	{
        	for (int j = 0; j < numTimes; j++) // loop over cols
        	{
            	surface << array[0+i][0+j] << ",";
            	//Morison_File << t << ";" << F << endl; 
        	}
        	surface << endl;
    	}
	}
	else cout << "Unable to open file";
	/*
	for (i1=1; i1<=nRows; i1++){
		i2 = (double)i1;
		dummy = rem(i2, 10.0);
		fprintf(surface, "%lf \n", dummy);
	}*/
    //fclose(surface);
	cout << "surface.csv saved." << endl;
	surface.close();
}

double uocOption::price(){
    int gslDim = numPrices - 2;
    gsl_vector *lower = gsl_vector_alloc (gslDim - 1); gsl_vector_set_all (lower, 0.0);
    gsl_vector *diag = gsl_vector_alloc (gslDim); gsl_vector_set_all (diag, 0.0);
    gsl_vector *upper = gsl_vector_alloc (gslDim - 1); gsl_vector_set_all (upper, 0.0);
    gsl_vector *rhs = gsl_vector_alloc (gslDim); gsl_vector_set_all (rhs, 0.0);
    gsl_vector *solution = gsl_vector_alloc (gslDim); gsl_vector_set_all (solution, 0.0);

    for (int i = 0; i < gslDim - 1; i++) {
        gsl_vector_set(lower, i, l());
        gsl_vector_set(upper, i, u());
    }

    for (int i = 0; i < gslDim; i++) {
        gsl_vector_set(diag, i, d(i+1));
    }

	/* // GSL print statements
	    for (int i = 0; i < gslDim ; i++) // gsl print statements
	    {   
	        printf ("diag_%d = %g\n", i, gsl_vector_get (diag, i));
	    }
	        for (int i = 0; i < gslDim - 1 ; i++) // gsl print statements
	    {   
	        printf ("lower_%d = %g\n", i, gsl_vector_get (lower, i));
	        printf ("upper_%d = %g\n", i, gsl_vector_get (upper, i));
	    }
	*/
    
    vector<double> w(numPrices,0.0);

	/*	Calculates GSL G Vectors
	    gsl_vector *gslG1LambdaN, *gslG1LambdaP, *gslG2LambdaN, *gslG2LambdaP;
	    
	    gslG1LambdaN = gsl_vector_alloc (numPrices - 1);
	    gslG1LambdaP = gsl_vector_alloc (numPrices - 1);
	    gslG2LambdaN = gsl_vector_alloc (numPrices - 1);
	    gslG2LambdaP = gsl_vector_alloc (numPrices - 1);

	    gsl_vector_set_zero(gslG1LambdaN);
	    gsl_vector_set_zero(gslG1LambdaP);
	    gsl_vector_set_zero(gslG2LambdaN);
	    gsl_vector_set_zero(gslG2LambdaP);
	*/

    vector<double> g1LambdaN(numPrices - 1, 0.0);
    vector<double> g1LambdaP(numPrices - 1, 0.0);
    vector<double> g2LambdaN(numPrices - 1, 0.0);
    vector<double> g2LambdaP(numPrices - 1, 0.0);
    
    for(int i = 0; i <= numPrices - 2; i++)  {// G vectors
    
        g1LambdaN[i] = g1(Y,(i+1) * deltax * lambdaN);
        //gsl_vector_set(gslG1LambdaN, i, g1(Y,(i+1) * deltax * lambdaN));
        //cout << "g1LambdaN[" << i << "] = " << g1LambdaN[i] << endl;
        //printf ("gslG1LambdaN[%d] = %g\n", i, gsl_vector_get (gslG1LambdaN, i));
        //cout << "g1(Y, " << i << " *deltax*lambdaN) = " << g1(Y, (i+1) * deltax * lambdaN) << endl;

        g1LambdaP[i] = g1(Y,(i+1) * deltax * lambdaP);
        //gsl_vector_set(gslG1LambdaP, i, g1(Y,(i+1) * deltax * lambdaP));
        //cout << "g1LambdaP[" << i << "] = " << g1LambdaP[i] << endl;
        //printf ("gslG1LambdaP[%d] = %g\n", i, gsl_vector_get (gslG1LambdaP, i));
        //cout << "g1(Y, " << i << " *deltax*lambdaP) = " << g1(Y, (i+1) * deltax * lambdaP) << endl;

        g2LambdaN[i] = g2(Y,(i+1) * deltax * lambdaN);
        //gsl_vector_set(gslG2LambdaN, i, g2(Y,(i+1) * deltax * lambdaN));
        //cout << "g2LambdaN[" << i << "] = " << g2LambdaN[i] << endl;
        //printf ("gslG2LambdaN[%d] = %g\n", i, gsl_vector_get (gslG2LambdaN, i));
        //cout << "g2(Y, " << i << " *deltax*lambdaN) = " << g2(Y, (i+1) * deltax * lambdaN) << endl;

        g2LambdaP[i] = g2(Y,(i+1) * deltax * lambdaP);
        //gsl_vector_set(gslG2LambdaP, i, g2(Y,(i+1) * deltax * lambdaP));
        //cout << "g2LambdaP[" << i << "] = " << g2LambdaP[i] << endl;
        //printf ("gslG2LambdaP[%d] = %g\n", i, gsl_vector_get (gslG2LambdaP, i));
        //cout << "g2(Y, " << i << " *deltax*lambdaP) = " << g2(Y, (i+1) * deltax * lambdaP) << endl;
    }
	
	double **grid = new double*[numPrices];
	for (int i = 0; i < numPrices; ++i){
		grid[i] = new double[numTimes];
	}
	// set array by looping over col then row
	for (int j = 0 ; j < numTimes; j++) {// loop over columns
    	for (int i = 0; i < numPrices; i++) {// loop over rows
        		grid[i][j] = 0;
        	}
    }
    
	for (int j = 0 ; j < numTimes; j++) { // loop over columns
	//for (int j = 0 ; j < 2; j++) {
    	//cout << "j = " << j << endl;
        //cout << "\ntau = " << myOption.tau - ((myOption.numTimes-1-j) * myOption.deltaTau) << endl;
    	if (j == 0) {
    		for (int i = 1; i < numPrices - 1; i++) { // loop over rows
        		w[i] = payoff(LogSmin + i * deltax);
        		grid[i][j] = w[i];
				//cout << "w[" << i  << "] = " << w[i] << endl;
        	}
        	for (int i = 1; i < numPrices - 1; i++) { // loop over rows
        		gsl_vector_set(rhs, i - 1, w[i] + (deltaTau / nu) * 
        			R(i, j, grid,g1LambdaN,g1LambdaP,g2LambdaN,g2LambdaP));
                //printf ("rhs[%d] = %g\n", i-1, gsl_vector_get(rhs, i-1));
        	}
        }
    	else {
    		gsl_linalg_solve_tridiag(diag,upper,lower,rhs,solution);
        	for (int i = 1; i < numPrices - 1 ; i++){ // loop over rows
        		w[i] = gsl_vector_get(solution, i - 1);
        		grid[i][j] = w[i];
        		//cout << "w[" << i  << "] = " << w[i] << endl;
        	}
        	for (int i = 1; i < numPrices - 1 ; i++) {// loop over rows
                gsl_vector_set(rhs,i-1,w[i] + (deltaTau / nu) *
                	R(i, j, grid, g1LambdaN, g1LambdaP, g2LambdaN, g2LambdaP));
               // printf ("gslKnown[%d] = %g\n", i-1, gsl_vector_get(gslKnown, i-1));
        	}
        }
    }
	//cout << "\ngrid Array Destructor" << endl;
	for(int i = 0; i < numPrices; ++i){
    	delete [] grid[i];
	}
	delete [] grid; 

	gsl_vector_free(lower);
	gsl_vector_free(diag);
	gsl_vector_free(upper);
	gsl_vector_free(rhs);
	gsl_vector_free(solution);

	for (int i = 0; i < numPrices; i++){
		logSpot = LogSmin + i * deltax;
		if (LogSmin + i * deltax < x && LogSmin + ((i+1) * deltax) > x){
			double y = w[i] + (w[i + 1] - w[i]) *
				((S - exp(LogSmin + (i * deltax)) ) /
				((exp(LogSmin + ((i+1) * deltax)) - exp(LogSmin + (i * deltax)))));
			// cout << "Interpolated price: $" << y << endl;
			return y;
		}
	}
}

double uocOption::pricePrint(){
    int gslDim = numPrices - 2;
    gsl_vector *lower = gsl_vector_alloc (gslDim - 1); gsl_vector_set_all (lower, 0.0);
    gsl_vector *diag = gsl_vector_alloc (gslDim); gsl_vector_set_all (diag, 0.0);
    gsl_vector *upper = gsl_vector_alloc (gslDim - 1); gsl_vector_set_all (upper, 0.0);
    gsl_vector *rhs = gsl_vector_alloc (gslDim); gsl_vector_set_all (rhs, 0.0);
    gsl_vector *solution = gsl_vector_alloc (gslDim); gsl_vector_set_all (solution, 0.0);

    for (int i = 0; i < gslDim - 1; i++) {
        gsl_vector_set(lower, i, l());
        gsl_vector_set(upper, i, u());
    }

    for (int i = 0; i < gslDim; i++) {
        gsl_vector_set(diag, i, d(i+1));
    }

	/* // GSL print statements
	    for (int i = 0; i < gslDim ; i++) // gsl print statements
	    {   
	        printf ("diag_%d = %g\n", i, gsl_vector_get (diag, i));
	    }
	        for (int i = 0; i < gslDim - 1 ; i++) // gsl print statements
	    {   
	        printf ("lower_%d = %g\n", i, gsl_vector_get (lower, i));
	        printf ("upper_%d = %g\n", i, gsl_vector_get (upper, i));
	    }
	*/
    
    vector<double> w(numPrices,0.0);

	/*	Calculates GSL G Vectors
	    gsl_vector *gslG1LambdaN, *gslG1LambdaP, *gslG2LambdaN, *gslG2LambdaP;
	    
	    gslG1LambdaN = gsl_vector_alloc (numPrices - 1);
	    gslG1LambdaP = gsl_vector_alloc (numPrices - 1);
	    gslG2LambdaN = gsl_vector_alloc (numPrices - 1);
	    gslG2LambdaP = gsl_vector_alloc (numPrices - 1);

	    gsl_vector_set_zero(gslG1LambdaN);
	    gsl_vector_set_zero(gslG1LambdaP);
	    gsl_vector_set_zero(gslG2LambdaN);
	    gsl_vector_set_zero(gslG2LambdaP);
	*/

    vector<double> g1LambdaN(numPrices - 1, 0.0);
    vector<double> g1LambdaP(numPrices - 1, 0.0);
    vector<double> g2LambdaN(numPrices - 1, 0.0);
    vector<double> g2LambdaP(numPrices - 1, 0.0);
    
    for(int i = 0; i <= numPrices - 2; i++)  {// G vectors
    
        g1LambdaN[i] = g1(Y,(i+1) * deltax * lambdaN);
        //gsl_vector_set(gslG1LambdaN, i, g1(Y,(i+1) * deltax * lambdaN));
        // cout << setprecision(10)  << "g1LambdaN[" << i << "] = " << g1LambdaN[i] << endl;
        //printf ("gslG1LambdaN[%d] = %g\n", i, gsl_vector_get (gslG1LambdaN, i));
        //cout << "g1(Y, " << i << " *deltax*lambdaN) = " << g1(Y, (i+1) * deltax * lambdaN) << endl;

        g1LambdaP[i] = g1(Y,(i+1) * deltax * lambdaP);
        //gsl_vector_set(gslG1LambdaP, i, g1(Y,(i+1) * deltax * lambdaP));
        // cout << setprecision(10) << "g1LambdaP[" << i << "] = " << g1LambdaP[i] << endl;
        //printf ("gslG1LambdaP[%d] = %g\n", i, gsl_vector_get (gslG1LambdaP, i));
        //cout << "g1(Y, " << i << " *deltax*lambdaP) = " << g1(Y, (i+1) * deltax * lambdaP) << endl;

        g2LambdaN[i] = g2(Y,(i+1) * deltax * lambdaN);
        //gsl_vector_set(gslG2LambdaN, i, g2(Y,(i+1) * deltax * lambdaN));
        // cout << setprecision(10)  << "g2LambdaN[" << i << "] = " << g2LambdaN[i] << endl;
        //printf ("gslG2LambdaN[%d] = %g\n", i, gsl_vector_get (gslG2LambdaN, i));
        //cout << "g2(Y, " << i << " *deltax*lambdaN) = " << g2(Y, (i+1) * deltax * lambdaN) << endl;

        g2LambdaP[i] = g2(Y,(i+1) * deltax * lambdaP);
        //gsl_vector_set(gslG2LambdaP, i, g2(Y,(i+1) * deltax * lambdaP));
        // cout << setprecision(10)  << "g2LambdaP[" << i << "] = " << g2LambdaP[i] << endl;
        //printf ("gslG2LambdaP[%d] = %g\n", i, gsl_vector_get (gslG2LambdaP, i));
        //cout << "g2(Y, " << i << " *deltax*lambdaP) = " << g2(Y, (i+1) * deltax * lambdaP) << endl;
    }
	
	double **grid = new double*[numPrices];
	for (int i = 0; i < numPrices; ++i){
		grid[i] = new double[numTimes];
	}
	// set array by looping over col then row
	for (int j = 0 ; j < numTimes; j++) {// loop over columns
    	for (int i = 0; i < numPrices; i++) {// loop over rows
        		grid[i][j] = 0;
        	}
    }
    
	for (int j = 0 ; j < numTimes; j++) { // loop over columns
	//for (int j = 0 ; j < 2; j++) {
    	//cout << "j = " << j << endl;
        //cout << "\ntau = " << myOption.tau - ((myOption.numTimes-1-j) * myOption.deltaTau) << endl;
    	if (j == 0) {
    		for (int i = 1; i < numPrices - 1; i++) { // loop over rows
        		w[i] = payoff(LogSmin + i * deltax);
        		grid[i][j] = w[i];
				// cout << setprecision(10) << "w[" << i  << "] = " << w[i] << endl;
        	}
        	for (int i = 1; i < numPrices - 1; i++) { // loop over rows
        		//cin.get(); 
        		gsl_vector_set(rhs, i - 1, w[i] + (deltaTau / nu) * 
        			R(i, j, grid,g1LambdaN,g1LambdaP,g2LambdaN,g2LambdaP));
        		// cout << setprecision(10) << "R: " << R(i, j, grid, g1LambdaN, g1LambdaP, g2LambdaN, g2LambdaP) << endl;
                // printf ("rhs[%d] = %g\n", i-1, gsl_vector_get(rhs, i-1));
        	}
        }
    	else {
    		gsl_linalg_solve_tridiag(diag,upper,lower,rhs,solution);
        	for (int i = 1; i < numPrices - 1 ; i++){ // loop over rows
        		w[i] = gsl_vector_get(solution, i - 1);
        		grid[i][j] = w[i];
        		// cout << setprecision(10) <<  "w[" << i  << "] = " << w[i] << endl;
        	}
        	for (int i = 1; i < numPrices - 1 ; i++) {// loop over rows
        		//cin.get();
                gsl_vector_set(rhs,i-1,w[i] + (deltaTau / nu) *
                	R(i, j, grid, g1LambdaN, g1LambdaP, g2LambdaN, g2LambdaP));
                // cout << setprecision(10) << "R: " << R(i, j, grid, g1LambdaN, g1LambdaP, g2LambdaN, g2LambdaP) << endl;
               // printf ("rhs[%d] = %g\n", i-1, gsl_vector_get(rhs, i-1));
        	}
        }
    }
    
    //writeCSV(grid);
	
	//cout << "\ngrid Array Destructor" << endl;
	for(int i = 0; i < numPrices; ++i){
    	delete [] grid[i];
	}
	delete [] grid;

	gsl_vector_free(lower);
	gsl_vector_free(diag);
	gsl_vector_free(upper);
	gsl_vector_free(rhs);
	gsl_vector_free(solution);
	
	for (int i = 0; i < numPrices; i++){
		double y = 0.0;
		logSpot = LogSmin + i * deltax;
		if (LogSmin + i * deltax < x && LogSmin + ((i+1) * deltax) > x){
			double y = w[i] + (w[i + 1] - w[i]) *
				((S - exp(LogSmin + (i * deltax)) ) /
				((exp(LogSmin + ((i+1) * deltax)) - exp(LogSmin + (i * deltax)))));
			cout << "Interpolated price: $" << y << endl;
			return y;
		}
	}
}

double uocOption::payoff(double x) const {
	
	return  max(exp(x) - K, 0.0);
}

double uocOption::R(int i, int j, double **array, 
	const vector<double>& g1LambdaN, const vector<double>& g1LambdaP,
	const vector<double>& g2LambdaN, const vector<double>& g2LambdaP) {

	double firstSum = 0.0;
	double secondSum = 0.0;
	double thirdSum = 0.0;
	double fourthSum = 0.0;

	if (i == 1){}
	else {
		for (int k = 1; k <= (i - 1); k++){	
		  	firstSum = firstSum + pow(lambdaN,Y) * (array[i - k][j] - array[i][j] - 
				(k * (array[i - k - 1][j] - array[i - k][j]))) * 
		  		(g2LambdaN[k-1] - g2LambdaN[k]);

			secondSum = secondSum + (array[i - k - 1][j] - array[i - k][j]) /
				(pow(lambdaN,1 - Y) * deltax) *
				(g1LambdaN[k-1] - g1LambdaN[k]);
		}
	}

	if (i == (numPrices - 2)){}
	else {
		for (int k = 1; k <= (numPrices - i-1 - 1); k++) {
			thirdSum = thirdSum + pow(lambdaP,Y) * (array[i + k][j] - array[i][j] - 
				(k * (array[i + k + 1][j] - array[i + k][j]))) * 
				(g2LambdaP[k-1] - g2LambdaP[k]);

			fourthSum = fourthSum + (((array[i + k + 1][j] - array[i + k][j]) / 
				(pow(lambdaP,(1 - Y)) * deltax)) *
				(g1LambdaP[k-1] - g1LambdaP[k]));
		}
	}
	double totalSum = firstSum + secondSum + thirdSum + fourthSum;
	return totalSum;
}

double uocOption::sigma2(double epsilon) {
	double result =  ((1 / nu) * pow(lambdaP,Y-2)) * (-1*(pow(lambdaP * epsilon,1-Y) *
		exp(-1 * lambdaP * epsilon)) + ((1 - Y) * (uocOption::g1(Y,0) - 
		uocOption::g1(Y,lambdaP * epsilon)))) + ((1 / nu) * pow(lambdaN,Y-2)) *
		(-1*(pow(lambdaN * epsilon,1-Y) *
		exp(-1*lambdaN * epsilon)) + ((1 - Y) * (uocOption::g1(Y,0) - 
		uocOption::g1(Y,lambdaN * epsilon)))) ;
	return result;
}

double uocOption::omega(double epsilon) {
	double result =  (((pow(lambdaP,Y))/nu) * uocOption::g2(Y,lambdaP * epsilon)) - 
		(((pow(lambdaP-1,Y))/nu) * uocOption::g2(Y,(lambdaP-1) * epsilon)) + 
		(((pow(lambdaN,Y))/nu) * uocOption::g2(Y,lambdaN * epsilon)) - 
		(((pow(lambdaN + 1,Y))/nu) * uocOption::g2(Y,(lambdaN + 1) * epsilon));
	return result;
}

double uocOption::g1(double Y, double x) {
    if (Y == 0)
        return exp(-x);
    else if (Y > 0 && Y < 1)
        return gsl_sf_gamma_inc(1-Y, x);
    else{
        printf("Error in g1(%f, %f)\n", Y, x);
        exit(1);
    }
}

double uocOption::g2(double Y, double x) {
    if (Y == 0)
        return gsl_sf_expint_E1(x);
    else if (Y > 0 && Y < 1)
        return (exp(-x) * pow(x, -Y) - gsl_sf_gamma_inc(1-Y, x)) / Y;
    else {
        printf("Error in g2(%f, %f)\n", Y, x);
        exit(2);
    }
}

double uocOption::l(){
	
	return -1 * Bl;
}  

double uocOption::d(int i){
	return 1 + r * deltaTau + Bl + Bu + deltaTau/nu * (((pow(lambdaN,Y)) * 
		(g2(Y,i * deltax * lambdaN))) + (pow(lambdaP,Y) * g2(Y,(numPrices-i) *
		deltax * lambdaP )));
}

double uocOption::u(){
	
	return -1 * Bu;
}

#endif