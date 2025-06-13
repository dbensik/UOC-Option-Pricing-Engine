
#ifndef UpOutCallOption_cpp
#define UpOutCallOption_cpp

#include "UpOutCallOption.hpp"
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

    lambdaN = sqrt((pow(theta, 2) / pow(sigma, 4)) + (2 / (sigma * sigma * nu))) + theta / (sigma * sigma);
    lambdaP = sqrt((pow(theta, 2) / pow(sigma, 4)) + (2 / (sigma * sigma * nu))) - theta / (sigma * sigma);

    Bl = ((sigma2(deltax) * deltaTau) / (2 * deltax * deltax)) - 
         (r - q + omega(deltax) - 0.5 * sigma2(deltax)) * (deltaTau / (2 * deltax));
    Bu = ((sigma2(deltax) * deltaTau) / (2 * deltax * deltax)) + 
         (r - q + omega(deltax) - 0.5 * sigma2(deltax)) * (deltaTau / (2 * deltax));
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

void UpOutCallOption::printGrid(double **array) const {
	cout << endl;
	for (int i = 0 ; i < numPrices; i++) {// loop over rows
        for (int j = 0; j < numTimes; j++) {// loop over cols
            cout << array[i][j] << " | ";
        }
        cout << endl;
    }
    cout << endl;
}

void UpOutCallOption::writeCSV(double **array) const {
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

std::tuple<gsl_vector*, gsl_vector*, gsl_vector*, gsl_vector*, gsl_vector*>
UpOutCallOption::initializeGslVectors(int gslDim) {
    gsl_vector* lower = gsl_vector_alloc(gslDim - 1);
    gsl_vector* diag = gsl_vector_alloc(gslDim);
    gsl_vector* upper = gsl_vector_alloc(gslDim - 1);
    gsl_vector* rhs = gsl_vector_alloc(gslDim);
    gsl_vector* solution = gsl_vector_alloc(gslDim);

    gsl_vector_set_all(lower, 0.0);
    gsl_vector_set_all(diag, 0.0);
    gsl_vector_set_all(upper, 0.0);
    gsl_vector_set_all(rhs, 0.0);
    gsl_vector_set_all(solution, 0.0);

    return std::make_tuple(lower, diag, upper, rhs, solution);
}

void uocOption::computeGVectors(std::vector<double>& g1LambdaN,
                                 std::vector<double>& g1LambdaP,
                                 std::vector<double>& g2LambdaN,
                                 std::vector<double>& g2LambdaP) const {
    for(int i = 0; i <= numPrices - 2; ++i) {
        g1LambdaN[i] = g1(Y, (i + 1) * deltax * lambdaN);
        g1LambdaP[i] = g1(Y, (i + 1) * deltax * lambdaP);
        g2LambdaN[i] = g2(Y, (i + 1) * deltax * lambdaN);
        g2LambdaP[i] = g2(Y, (i + 1) * deltax * lambdaP);
    }
}

double UpOutCallOption::price() {
    vector<vector<double>> grid(numTimes, vector<double>(numPrices, 0.0));

    // Set terminal condition (payoff)
    for (int j = 0; j < numPrices; ++j) {
        double logS = LogSmin + j * deltax;
        double S = exp(logS);
        grid[numTimes - 1][j] = std::max(S - K, 0.0);
    }

    // Boundary conditions
    for (int i = 0; i < numTimes; ++i) {
        grid[i][0] = 0.0;                 // S -> 0
        grid[i][numPrices - 1] = 0.0;     // S -> B (knocked out)
    }

    // Time-stepping backwards (placeholder logic)
    for (int i = numTimes - 2; i >= 0; --i) {
        for (int j = 1; j < numPrices - 1; ++j) {
            grid[i][j] = 0.5 * (grid[i + 1][j - 1] + grid[i + 1][j + 1]); // Crude average
        }
    }

    // Interpolation at x = log(S)
    int j = static_cast<int>((x - LogSmin) / deltax);
    double w = (x - (LogSmin + j * deltax)) / deltax;
    return (1 - w) * grid[0][j] + w * grid[0][j + 1];
}

double UpOutCallOption::pricePrint(){
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

double UpOutCallOption::payoff(double x) const {
	
	return  max(exp(x) - K, 0.0);
}

double UpOutCallOption::R(int i, int j, double **array, 
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

double UpOutCallOption::sigma2(double epsilon) {
	double result =  ((1 / nu) * pow(lambdaP,Y-2)) * (-1*(pow(lambdaP * epsilon,1-Y) *
		exp(-1 * lambdaP * epsilon)) + ((1 - Y) * (UpOutCallOption::g1(Y,0) - 
		UpOutCallOption::g1(Y,lambdaP * epsilon)))) + ((1 / nu) * pow(lambdaN,Y-2)) *
		(-1*(pow(lambdaN * epsilon,1-Y) *
		exp(-1*lambdaN * epsilon)) + ((1 - Y) * (UpOutCallOption::g1(Y,0) - 
		UpOutCallOption::g1(Y,lambdaN * epsilon)))) ;
	return result;
}

double UpOutCallOption::omega(double epsilon) {
	double result =  (((pow(lambdaP,Y))/nu) * UpOutCallOption::g2(Y,lambdaP * epsilon)) - 
		(((pow(lambdaP-1,Y))/nu) * UpOutCallOption::g2(Y,(lambdaP-1) * epsilon)) + 
		(((pow(lambdaN,Y))/nu) * UpOutCallOption::g2(Y,lambdaN * epsilon)) - 
		(((pow(lambdaN + 1,Y))/nu) * UpOutCallOption::g2(Y,(lambdaN + 1) * epsilon));
	return result;
}

double UpOutCallOption::g1(double Y, double x) {
    if (Y == 0)
        return exp(-x);
    else if (Y > 0 && Y < 1)
        return gsl_sf_gamma_inc(1-Y, x);
    else{
        printf("Error in g1(%f, %f)\n", Y, x);
        exit(1);
    }
}

double UpOutCallOption::g2(double Y, double x) {
    if (Y == 0)
        return gsl_sf_expint_E1(x);
    else if (Y > 0 && Y < 1)
        return (exp(-x) * pow(x, -Y) - gsl_sf_gamma_inc(1-Y, x)) / Y;
    else {
        printf("Error in g2(%f, %f)\n", Y, x);
        exit(2);
    }
}

double UpOutCallOption::l(){
	
	return -1 * Bl;
}  

double UpOutCallOption::d(int i){
	return 1 + r * deltaTau + Bl + Bu + deltaTau/nu * (((pow(lambdaN,Y)) * 
		(g2(Y,i * deltax * lambdaN))) + (pow(lambdaP,Y) * g2(Y,(numPrices-i) *
		deltax * lambdaP )));
}

double UpOutCallOption::u(){
	
	return -1 * Bu;
}

#endif
