// TestUoc.cpp

/* //  DASHBOARD
    This set of functions will be used to price and up-and-out
    call option.  By Solving the PIDE and using the functions for
    the difference equation, k(y), lambdaP, lambdaN, sigma^2, 
    omega, g1, g2.

    The difference equation will be solved for omega(x,tau) using
    a tridiagonal matrix solver from ??????

    The PIDE is solved using the following boundary and initial
    conditions:
    	Initial condition:
    	Boundary conditions:

    /// make file command
    Use pkg-config to compile and run program.cpp:
    $ g++ $(pkg-config --cflags gsl) -c TestUoc.cpp uocOption.cpp // to link to gsl.
    $ g++ TestUoc.o uocOption.o $(pkg-config --libs gsl) -o TestUoc
    $ ./TestUoc
*/

#include "uocOption.hpp"
#include "optimizer.hpp"
// #include "functions.hpp"

int main(int argc, char* argv[])
{
	// This code for measuring runtime was taken from Ali Hirsa
	struct timeval startTime, endTime;
    long seconds, useconds;
    double mtime;
	gettimeofday(&startTime, NULL);  // struct startTime passed by ref
    ///////////////////////////////////////////////////////////////////////////////
  // Corrected result from Case Study 1
    uocOption myOption;
    cout << "Corrected result from Case Study 1: " << endl;
    myOption.pricePrint();

    vector<double> myStrikes(4);
    vector<double> myBarriers(3);
    vector<double> myTestTheta(4);
    vector<double> myInitialTheta(4);
    vector<double> mySimplexTheta(4);
    vector<double> myOptimalTheta(4);
    
    myStrikes[0] = 100.0; myStrikes[1] = 105.0; myStrikes[2] = 110.0; myStrikes[3] = 120.0;
    myBarriers[0] = 125.0, myBarriers[1] = 135.0, myBarriers[2] = 145.0;
    myTestTheta[0] = 0.32, myTestTheta[1] = 0.42, myTestTheta[2] = -0.29, myTestTheta[3] = 0.74;
    myInitialTheta[0] = 0.1, myInitialTheta[1] = 0.2, myInitialTheta[2] = 0.1, myInitialTheta[3] = 0.1;

    vector< vector<double> > marketPrices(3, vector<double>(4));
    vector< vector<double> > modelPrices(3, vector<double>(4));

    getPriceGrid(myTestTheta, myStrikes, myBarriers, marketPrices);
    getPriceGrid(myInitialTheta, myStrikes, myBarriers, modelPrices);  // test purpose only
    
    cout << endl; 
    printTestTheta(myTestTheta);
    cout << "\nMarket Prices: " << endl;
    print2DVector(marketPrices); cout << endl;
    
    printInitialTheta(myInitialTheta);
    cout << "\nModel Prices: " << endl;
    print2DVector(modelPrices); cout << endl;
    
    char response;
    cout << "Do you want to run the grid search? (y/n): ";
    cin >> response;
    if (response == 'y'){ // Results from current grid search
        gridSearch(myInitialTheta, myOptimalTheta, myStrikes, myBarriers,
            marketPrices, modelPrices);
        mySimplexTheta[0] = myOptimalTheta[0];
        mySimplexTheta[1] = myOptimalTheta[1];
        mySimplexTheta[2] = myOptimalTheta[2];
        mySimplexTheta[3] = myOptimalTheta[3];
        // printSimplexTheta(mySimplexTheta);
    }
    else if (response == 'n'){ // Results from previous grid search
        mySimplexTheta[0] = 0.3;
        mySimplexTheta[1] = 0.4; 
        mySimplexTheta[2] = -0.25; 
        mySimplexTheta[3] = 0.75; 
        //printSimplexTheta(mySimplexTheta);
    }
    else {
        cout << "Enter a valid character" << endl;
    }

    cout << "\nDo you want to run Nelder Mead optimization? (y/n): ";
    cin >> response;
    if (response == 'y'){
        int num;
        cout << "\nWhich one? (1/2): ";
        cin >> num;
        char answer;
        if (num ==1){
            cout << "\nWhich Theta? Initial (i) or Optimized (o)? : ";
            cin >> answer;
            if (answer == 'i'){
                myInitialTheta[0] = 0.1;
                myInitialTheta[1] = 0.2;
                myInitialTheta[2] = 0.1;
                myInitialTheta[3] = 0.1;
                nelderMead1(myInitialTheta); // ~ mins,  iterations
            }
            else{
                nelderMead1(mySimplexTheta);
            }
        }
        else if (num == 2) {
            cout << "\nWhich Theta? Initial (i) or Optimized (o)? : ";
            cin >> answer;
            if (answer == 'i'){
                myInitialTheta[0] = 0.1;
                myInitialTheta[1] = 0.2;
                myInitialTheta[2] = 0.1;
                myInitialTheta[3] = 0.1;
                nelderMead2(myInitialTheta); // ~ mins,  iterations
            }
            else{
                nelderMead2(mySimplexTheta);
            }
        }
        else {
            cout << "\nEnter a valid number" << endl;
        }
    }
    else if (response == 'n'){}
    ////////////////////////////////////////////////////////////////
    gettimeofday(&endTime, NULL);
    seconds = endTime.tv_sec - startTime.tv_sec;
    useconds = endTime.tv_usec - startTime.tv_usec;
    mtime = ((seconds) * 1000.0 + useconds/1000.0); 
    cout << "\nTime elapsed was: " << mtime << " (milliseconds)" << endl;
    cout << endl;

	return 0;
}