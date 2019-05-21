#ifndef optimizer_cpp
#define optimizer_cpp

#include "optimizer.hpp"
#include "uocOption.hpp"

double vectorMin(const vector<double>& v){

  vector<double>::const_iterator it = min_element(v.begin(), v.end());
  return *it;
}

int vectorMinIndex(const vector<double>& v){

  vector<double>::const_iterator it = min_element(v.begin(), v.end());
  return it - v.begin();
}

void printVector(const vector<double>& v){
	cout << "( ";
    for (vector<double>::const_iterator it = v.begin(); it != v.end(); ++it) {
        cout  << *it  << " ";}
        cout << ")" << endl; 
}

void print2DVector(const vector< vector<double> > &vec2D){
	for (int i = 0 ; i < vec2D.size(); i++){ // loop over rows
		for (int j = 0; j < vec2D[0].size(); j++){ // loop over cols
            cout << vec2D[i][j] << " | ";
         }
         cout << endl;
    }
}

double round(double val, int precision){
    std::stringstream s;
    s << std::setprecision(precision) << std::setiosflags(std::ios_base::fixed) << val;
    s >> val;
    return val;
}

void testTridiagonal(){

	// static const double a[] = { 0.0, -1.0, -1.0, -1,0 };
    // vector<double> aVec (a, a + sizeof(a) / sizeof(a[0]) );

    // static const double b[] = { 4.0, 4.0, 4.0, 4.0 };
    // vector<double> bVec (b, b + sizeof(b) / sizeof(b[0]) );

    // static const double c[] = { -1.0, -1.0, -1,0, 0.0 };
    // vector<double> cVec (c, c + sizeof(c) / sizeof(c[0]) );

    // static const double d[] = { 5.0, 5.0, 10.0, 23.0 };
    // vector<double> dVec (d, d + sizeof(d) / sizeof(d[0]) );

    // static const double f[] = { 5.0, 5.0, 10.0, 23.0 };
    // vector<double> fVec (f, f + sizeof(f) / sizeof(f[0]) );

    // // results    { 2,  3,  5, 7  }
    // solve(aVec,bVec,cVec,dVec, fVec);
    // int gslDim = 4;
    // gsl_vector *l = gsl_vector_alloc (gslDim - 1); gsl_vector_set_all (l, -1.0);
    // gsl_vector *d = gsl_vector_alloc (gslDim); gsl_vector_set_all (d, 4.0);
    // gsl_vector *u = gsl_vector_alloc (gslDim - 1); gsl_vector_set_all (u, -1.0);
    // gsl_vector *rhs = gsl_vector_alloc (gslDim); gsl_vector_set_all (rhs, 1.0);
    // gsl_vector *solution = gsl_vector_alloc (gslDim); gsl_vector_set_zero (solution);
    
    // gsl_vector_set(rhs,0,5.0);
    // gsl_vector_set(rhs,1,5.0);
    // gsl_vector_set(rhs,2,10.0);
    // gsl_vector_set(rhs,3,23.0);
    
    // gsl_linalg_solve_tridiag(d,u,l,rhs,solution);
    // for ( int i = 0; i < gslDim; i++){
    //     cout << gsl_vector_get(solution,i) << endl;
    // }

    // gsl_vector_free(l);
    // gsl_vector_free(d);
    // gsl_vector_free(u);
    // gsl_vector_free(rhs);
    // gsl_vector_free(solution);
}

void printTheta(const vector<double>& Theta){
	cout << "Theta: " << endl;
	cout << "    sigma: " << Theta[0] << ", ";
	cout << "nu: " << Theta[1] << ", ";
	cout << "theta: " << Theta[2] << ", ";
	cout << "Y: " << Theta[3] << endl;
}

void printSimplexTheta(const vector<double>& mySimplexTheta){
	cout << "Simplex Theta: " << endl;
	cout << "    sigma: " << mySimplexTheta[0] << ", ";
	cout << "nu: " << mySimplexTheta[1] << ", ";
	cout << "theta: " << mySimplexTheta[2] << ", ";
	cout << "Y: " << mySimplexTheta[3] << endl;
}

void printInitialTheta(const vector<double>& myInitialTheta){
	cout << "Initial Theta: " << endl;
	cout << "    sigma: " << myInitialTheta[0] << ", ";
	cout << "nu: " << myInitialTheta[1] << ", ";
	cout << "theta: " << myInitialTheta[2] << ", ";
	cout << "Y: " << myInitialTheta[3] << endl;
}

void printTestTheta(const vector<double>& myTestTheta){
	cout << "Test Theta: " << endl;
	cout << "    sigma: " << myTestTheta[0] << ", ";
	cout << "nu: " << myTestTheta[1] << ", ";
	cout << "theta: " << myTestTheta[2] << ", ";
	cout << "Y: " << myTestTheta[3] << endl;
}

void printOptimalTheta(const vector<double>& myOptimalTheta){
	cout << "Optimal Theta: " << endl;
	cout << "    sigma: " << myOptimalTheta[0] << ", ";
	cout << "nu: " << myOptimalTheta[1] << ", ";
	cout << "theta: " << myOptimalTheta[2] << ", ";
	cout << "Y: " << myOptimalTheta[3] << endl;
}

void printConstrainedTheta(const vector<double>& myInitialTheta){
	cout << "Constrained Theta: " << endl;
	cout << "    sigma: " << myInitialTheta[0] << ", ";
	cout << "nu: " << myInitialTheta[1] << ", ";
	cout << "theta: " << myInitialTheta[2] << ", ";
	cout << "Y: " << myInitialTheta[3] << endl;
}

double constrainSigma(double& sigma){

	double low = 0.1;
	double high = 0.6;

	if (sigma >= low && sigma <= high){
		return sigma;
	}
	else{
		double range = high - low;
		int n = floor((sigma - low) / range);
		if (n%2 == 0){
			return sigma - n * range;
		}
		else if (n%2 != 0) {
			return high + n * range - (sigma - low);
		}
	}
}

double constrainNu(double& nu){

	double low = 0.1;
	double high = 0.8;

	if (nu >= low && nu <= high){
		return nu;
	}
	else{
		double range = high - low;
		int n = floor((nu - low) / range);
		if (n%2 == 0){
			return nu - n * range;
		}
		else if (n%2 != 0) {
			return high + n * range - (nu - low);
		}
	}
}

double constrainTheta(double& theta){

	double low = -0.5;
	double high = 0.5;

	if (theta >= low && theta <= high){
		return theta;
	}
	else{
		double range = high - low;
		int n = floor((theta - low) / range);
		if (n%2 == 0){
			return theta - n * range;
		}
		else if (n%2 != 0) {
			return high + n * range - (theta - low);
		}
	}
}

double constrainY(double& Y){

	double low = 0.1;
	double high = 0.9;

	if (Y >= low && Y <= high){
		return Y;
	}
	else{
		double range = high - low;
		int n = floor((Y - low) / range);
		if (n%2 == 0){
			return Y - n * range;
		}
		else if (n%2 != 0) {
			return high + n * range - (Y - low);
		}
	}
}

void testConstraints(){
	for (double i = -0.6; i <= 1.1; i += 0.1){
	    cout << "i: " << i << endl;
	    cout << "cSigma: " << constrainSigma(i) << endl;
	    cout << "i: " << i << endl;        
	    cout << "cNu: " << constrainNu(i) << endl;
	    cout << "i: " << i << endl;
	    cout << "cTheta: " << constrainTheta(i) << endl;
	    cout << "i: " << i << endl;
	    cout << "cY: " << constrainY(i) << endl;
	    cin.get();
    }
}

void getPriceGrid(const vector<double>& Theta,
	const vector<double>& myStrikes,
	const vector<double>& myBarriers,
	vector<vector<double> >& modelPrices){
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 4; j++) {
			uocOption myOption(myStrikes[j],myBarriers[i],Theta);
            modelPrices[i][j] = myOption.price();
		}
	}
}

void gridSearch(vector<double>& myInitialTheta,
	vector<double>& myOptimalTheta,
	const vector<double>& myStrikes,
	const vector<double>& myBarriers,
	const vector< vector<double> >& marketPrices,
	vector<vector<double> >& modelPrices){
	
	struct timeval startTime, endTime;
    long seconds, useconds;
    double mtime;
	gettimeofday(&startTime, NULL);

	double numSteps = 5.0;
	int numIter = 0;
	
	double sigmaLow = 0.0;
	double sigmaHigh = 0.6;
    double sigmaStep = (sigmaHigh - sigmaLow) / (numSteps - 1);
    
    double nuLow = 0.0;
    double nuHigh = 0.8;
    double nuStep = (nuHigh - nuLow) / (numSteps - 1);
    
    double thetaLow = -0.5;
    double thetaHigh = 0.5;
    double thetaStep = (thetaHigh - thetaLow) / (numSteps - 1);  //theta = 0 OK
    
    double yLow = 0.0; // Y = 0 OK
    double yHigh = 1.0; // Y cannot = 1
    double yStep = (yHigh - yLow) / (numSteps - 1);

	double minDiff = 100000000.0;

	// cout << "\nIteration: " << endl;
	#pragma omp parallel
	for (int i1 = 0; i1 <= 4; i1++){
		 for (int i2 = 0; i2 <= 4; i2++){
			for (int i3 = 0; i3 <= 4; i3++){
		 		for (int i4 = 0; i4 <= 4; i4++){
					numIter++;
					cout << "#" << numIter << " ";					
					myInitialTheta[0] = sigmaLow + sigmaStep * (i1);
					myInitialTheta[1] = nuLow + nuStep * (i2);
					myInitialTheta[2] = thetaLow + thetaStep * (i3);
					myInitialTheta[3] = yLow + yStep * (i4);
					// printInitialTheta(myInitialTheta);
					
					myInitialTheta[0] = constrainSigma(myInitialTheta[0]);
					myInitialTheta[1] = constrainNu(myInitialTheta[1]);
					myInitialTheta[2] = constrainTheta(myInitialTheta[2]);
					myInitialTheta[3] = constrainY(myInitialTheta[3]);
					// cout << endl; printConstrainedTheta(myInitialTheta);
					
					getPriceGrid(myInitialTheta,myStrikes,myBarriers,modelPrices);
					
					double diff = ssd(modelPrices,marketPrices);
					// cout << "minDiff: " << minDiff << endl;
					// cout << "diff: " << diff << endl;
					if (diff < minDiff){
						minDiff = diff;
						myOptimalTheta = myInitialTheta;
					}
					// printOptimalTheta(myOptimalTheta);
					// cin.get();
				}
			}
		}
	}
	cout << endl; printOptimalTheta(myOptimalTheta);

	gettimeofday(&endTime, NULL);
    seconds = endTime.tv_sec - startTime.tv_sec;
    useconds = endTime.tv_usec - startTime.tv_usec;
    mtime = ((seconds) * 1000.0 + useconds/1000.0); 
    cout << "\nTime elapsed for gridSearch was: " << mtime << " (milliseconds)" << endl;
    cout << endl;
}

void ssre(const vector<vector<double> >& modelPrices,
	const vector<vector<double> >& marketPrices,
	vector<double>& diffVector){
	double diff = 0.0;
	for (int i = 0; i < 3 ; i++) { // loop over barriers
        for (int j = 0; j < 4; j++) { // loop over strikes
        		diff += pow(abs(((modelPrices[i][j] - marketPrices[i][j]) / 
        			marketPrices[i][j])),2);
            }
        }
	for(double& d : diffVector)
  		d += diff;
}

double ssre(const vector<vector<double> >& modelPrices,
	const vector<vector<double> >& marketPrices){
	double diff = 0.0;
	for (int i = 0; i < 3 ; i++) { // loop over barriers
        for (int j = 0; j < 4; j++) { // loop over strikes
        		diff += pow(abs(((modelPrices[i][j] - marketPrices[i][j]) / 
        			marketPrices[i][j])),2);
            }
        }
    return diff;
}

void ssre1(const vector<vector<double> >& modelPrices,
	const vector<vector<double> >& marketPrices,
	vector<double>& diffVector){
	double diff = 0.0;
	for (int i = 0; i < 3 ; i++) { // loop over barriers
        for (int j = 0; j < 4; j++) { // loop over strikes
        		diff += abs((modelPrices[i][j] - marketPrices[i][j]) / 
        			marketPrices[i][j]);
            }
        }
	for(double& d : diffVector)
  		d += diff;
}

double ssre1(const vector<vector<double> >& modelPrices,
	const vector<vector<double> >& marketPrices){
	double diff = 0.0;
	for (int i = 0; i < 3 ; i++) { // loop over barriers
        for (int j = 0; j < 4; j++) { // loop over strikes
        		diff += abs((modelPrices[i][j] - marketPrices[i][j]) / 
        			marketPrices[i][j]);
            }
        }
	return diff;
}

void ssre2(const vector<vector<double> >& modelPrices,
	const vector<vector<double> >& marketPrices,
	vector<double>& diffVector){
	double diff = 0.0;
	for (int i = 0; i < 3 ; i++) { // loop over barriers
        for (int j = 0; j < 4; j++) { // loop over strikes
        		diff += abs(modelPrices[i][j] - marketPrices[i][j]);
            }
        }
	for(double& d : diffVector)
  		d += diff; // closest
}

double ssre2(const vector<vector<double> >& modelPrices,
	const vector<vector<double> >& marketPrices){
	double diff = 0.0;
	for (int i = 0; i < 3 ; i++) { // loop over barriers
        for (int j = 0; j < 4; j++) { // loop over strikes
        		diff += abs(modelPrices[i][j] - marketPrices[i][j]);
            }
        }
	return diff;
}

double ssd(const vector<vector<double> >& modelPrices,
	const vector<vector<double> >& marketPrices){
	double diff = 0.0;
	for (int i = 0; i < 3 ; i++) { // loop over barriers
        for (int j = 0; j < 4; j++) { // loop over strikes
        		diff += pow(modelPrices[i][j] - marketPrices[i][j],2);
            }
        }
	return diff;
}

double objective(vector<double>& mySimplexTheta){

	// printTheta(mySimplexTheta);
	// Constrain variables if necessary

	// double objDiff = 0.0;
	vector< vector<double> > marketPrices(3, vector<double>(4));
	vector< vector<double> > modelPrices(3, vector<double>(4));
	vector<double> myStrikes(4);
    vector<double> myBarriers(3);
    vector<double> myTestTheta(4);
	myStrikes[0] = 100.0; myStrikes[1] = 105.0; myStrikes[2] = 110.0; myStrikes[3] = 120.0;
    myBarriers[0] = 125.0, myBarriers[1] = 135.0, myBarriers[2] = 145.0;
    myTestTheta[0] = 0.32, myTestTheta[1] = 0.42, myTestTheta[2] = -0.29, myTestTheta[3] = 0.74;
	getPriceGrid(mySimplexTheta,myStrikes,myBarriers,modelPrices); //model
	getPriceGrid(myTestTheta,myStrikes,myBarriers,marketPrices); // market
	// print2DVector(marketPrices);
	// print2DVector(modelPrices);
	// objDiff = ssre(modelPrices,marketPrices);
	// cout << "Difference: " << objDiff << endl;
	// return objDiff;
	return ssre(modelPrices, marketPrices);
}

double myFunc(const gsl_vector *v, void *params){
    double x, y;
    double *dp = (double *)params;

    x = gsl_vector_get(v, 0);
    y = gsl_vector_get(v, 1);

    //cout << "x: " << x << ", y: " << y << endl; // changing
    //cout << "dp[0]: " << dp[0] << ", dp[1]: " << dp[1] << endl; // not changing

    // Paraboloid centered on (dp[0],dp[1]), with scale factors (3.0, 5.0) and minimum 10 
    return 3.0*(x-dp[0])*(x-dp[0]) + 5.0*(y-dp[1])*(y-dp[1]) + 10.0;
}

double myFunc2(const gsl_vector *v, void *params){
    //double sigmaF, nuF, thetaF, yF;
    double *dp = (double *)params;

    vector<double> mySimplexTheta(4);
    vector<double> myTestTheta(4);

    mySimplexTheta[0] = gsl_vector_get(v, 0);
    mySimplexTheta[1] = gsl_vector_get(v, 1);
    mySimplexTheta[2] = gsl_vector_get(v, 2);
    mySimplexTheta[3] = gsl_vector_get(v, 3);

    mySimplexTheta[0] = constrainSigma(mySimplexTheta[0]);
    mySimplexTheta[1] = constrainNu(mySimplexTheta[1]);
    mySimplexTheta[2] = constrainTheta(mySimplexTheta[2]);
    mySimplexTheta[3] = constrainY(mySimplexTheta[3]);

    myTestTheta[0] = dp[0];
    myTestTheta[1] = dp[1];
    myTestTheta[2] = dp[2];
    myTestTheta[3] = dp[3];

    vector<double> myStrikes(4);
    vector<double> myBarriers(3);
    myStrikes[0] = 100.0; myStrikes[1] = 105.0; myStrikes[2] = 110.0; myStrikes[3] = 120.0;
    myBarriers[0] = 125.0, myBarriers[1] = 135.0, myBarriers[2] = 145.0;
    
    vector< vector<double> > marketPrices(3, vector<double>(4));
    vector< vector<double> > modelPrices(3, vector<double>(4));
    getPriceGrid(mySimplexTheta,myStrikes,myBarriers,modelPrices); //model
    getPriceGrid(myTestTheta,myStrikes,myBarriers,marketPrices); // market
    
    double objDiff = 0.0;
    objDiff = ssd(modelPrices,marketPrices);

    return objDiff;
}

int nelderMeadExample(){

	double par[2] = {1.0, 2.0};
	//const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
    
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minex_func;
	size_t iter = 0;
	int status;
	double size;

	/* Starting point */
	x = gsl_vector_alloc (2);
	gsl_vector_set (x, 0, 5.0);
	gsl_vector_set (x, 1, 7.0);

	/* Set initial step sizes to 1 */
	ss = gsl_vector_alloc (2);
	gsl_vector_set_all (ss, 1.0);

	/* Initialize method and iterate */
	minex_func.n = 2;
	minex_func.f = myFunc;
	minex_func.params = par;
	s = gsl_multimin_fminimizer_alloc (T, 2);
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
	do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if (status)
			break;
		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 1e-4);
		if (status == GSL_SUCCESS) {
			printf ("converged to minimum at\n");
		}
		printf ("%5zu %9.3e %9.3e f() = %8.3f size = %.3f\n",
			iter,
			gsl_vector_get (s->x, 0),
			gsl_vector_get (s->x, 1),
			s->fval, size);
	}
	while (status == GSL_CONTINUE && iter < 100);

	gsl_vector_free(x);
	gsl_vector_free(ss);

	gsl_multimin_fminimizer_free (s);

	return status;
}

int nelderMead1(const vector<double>& Theta){
	struct timeval startTime, endTime;
    long seconds, useconds;
    double mtime;
	gettimeofday(&startTime, NULL);

    double par[4] = {0.32,0.42,-0.29,0.74};
	// const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;  
	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex; 

    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;
    size_t iter = 0;
    int status;
    double size;

    // /* Starting point */
    x = gsl_vector_alloc(4);
    gsl_vector_set(x,0,Theta[0]);
    gsl_vector_set(x,1,Theta[1]);
    gsl_vector_set(x,2,Theta[2]);
    gsl_vector_set(x,3,Theta[3]);

    /* Set initial step sizes to 0.1 */
    ss = gsl_vector_alloc(4);
    gsl_vector_set_all(ss, 0.1);

    /* Initialize method and iterate */
    minex_func.n = 4;
    minex_func.f = myFunc2;
    minex_func.params = par;
    s = gsl_multimin_fminimizer_alloc(T,4);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    do {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        if (status)
            break;
        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, 1e-7);
        if (status == GSL_SUCCESS) {
            printf ("\nconverged to minimum at\n");
        }
        printf("%5zu %9.3e %9.3e %9.3e %9.3e f() = %8.3e size = %.8f\n", 
            iter,
            gsl_vector_get (s->x, 0),
            gsl_vector_get (s->x, 1),
            gsl_vector_get (s->x, 2),
            gsl_vector_get (s->x, 3),
            s->fval, size);
    	//}
    }
    while (status == GSL_CONTINUE && iter < 1000);

    gsl_vector_free(x);
    gsl_vector_free(ss);

    gsl_multimin_fminimizer_free (s);

    gettimeofday(&endTime, NULL);
    seconds = endTime.tv_sec - startTime.tv_sec;
    useconds = endTime.tv_usec - startTime.tv_usec;
    mtime = ((seconds) * 1000.0 + useconds/1000.0); 
    cout << "\nTime elapsed for nelderMead1 was: " << mtime << " (milliseconds)" << endl;
    cout << endl;

    return status;
}

int nelderMead2(const vector<double>& Theta){
	struct timeval startTime, endTime;
    long seconds, useconds;
    double mtime;
	gettimeofday(&startTime, NULL);

    double par[4] = {0.32,0.42,-0.29,0.74};
	// const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex; 

    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;
    size_t iter = 0;
    int status;
    double size;

    // /* Starting point */
    x = gsl_vector_alloc(4);
    gsl_vector_set(x,0,Theta[0]);
    gsl_vector_set(x,1,Theta[1]);
    gsl_vector_set(x,2,Theta[2]);
    gsl_vector_set(x,3,Theta[3]);

    /* Set initial step sizes to 0.1 */
    ss = gsl_vector_alloc(4);
    gsl_vector_set_all(ss, 0.1);

    /* Initialize method and iterate */
    minex_func.n = 4;
    minex_func.f = myFunc2;
    minex_func.params = par;
    s = gsl_multimin_fminimizer_alloc(T,4);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    do {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        if (status)
            break;
        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, 1e-7);
        if (status == GSL_SUCCESS) {
            printf ("\nconverged to minimum at\n");
        }
        double sigma = gsl_vector_get (s->x, 0);
        double nu = gsl_vector_get (s->x, 1);
        double theta = gsl_vector_get (s->x, 2);
        double Y = gsl_vector_get (s->x, 3);
        printf("U: %5zu %9.3e %9.3e %9.3e %9.3e f() = %8.3e size = %.8f\n", // unconstrained
            iter,
            gsl_vector_get (s->x, 0),
            gsl_vector_get (s->x, 1),
            gsl_vector_get (s->x, 2),
            gsl_vector_get (s->x, 3),
            s->fval, size);
        printf("C: %5zu %9.3e %9.3e %9.3e %9.3e f() = %8.3e size = %.8f\n", // constrained
            iter,
            constrainSigma(sigma),
            constrainNu(nu),
            constrainTheta(theta),
            constrainY(Y),
            s->fval, size);
    	//}
    }
    while (status == GSL_CONTINUE && iter < 1000);

    gsl_vector_free(x);
    gsl_vector_free(ss);

    gsl_multimin_fminimizer_free (s);

    gettimeofday(&endTime, NULL);
    seconds = endTime.tv_sec - startTime.tv_sec;
    useconds = endTime.tv_usec - startTime.tv_usec;
    mtime = ((seconds) * 1000.0 + useconds/1000.0); 
    cout << "\nTime elapsed for nelderMead2 was: " << mtime << " (milliseconds)" << endl;
    cout << endl;

    return status;
}

#endif