# UOC Option Pricing Engine 🧠📈
This C++ project is an academic-grade implementation of an Up-and-Out Call Option pricing engine, built to solve a Partial Integro-Differential Equation (PIDE) using finite difference methods and tridiagonal solvers. In short: it’s a bunch of math that spits out the price of fancy financial derivatives while giving your CPU a mild workout.

## 💡 Purpose
The model implements:

* A finite difference scheme to numerically solve the PIDE

* Support for parameterized option characteristics (strike, barrier, etc.)

* Market-model calibration using:

* Grid search

* Nelder-Mead simplex optimization (via GSL)

* Output via CSV for 3D surface plotting

* Runtime metrics, because you deserve to know how long it took

## 🧰 Features
* PIDE solver with tridiagonal matrix linear algebra

* GSL-backed numerical routines

* Multiple constructors for flexible model configuration

* Market calibration routines

* Full debugging output (opt-in if you like pain)

* Structured modularity (finally)

* Runtime profiling (Ali Hirsa style)

## 📦 Project Structure
```bash
.
├── TestUoc.cpp         # Main driver file with I/O and runtime logic
├── uocOption.hpp/cpp   # Main option class with pricing engine and solver logic
├── optimizer.hpp/cpp   # Optimizers, grid search, Nelder-Mead, evaluation tools
├── Makefile            # Use pkg-config to compile with GSL
└── README.md           # The thing you're reading
```

## 🧪 Dependencies
* C++11 or later

* GSL (GNU Scientific Library)

* A functioning brain (optional but recommended)

🧰 Building the Project
```bash
g++ $(pkg-config --cflags gsl) -c TestUoc.cpp uocOption.cpp optimizer.cpp
g++ TestUoc.o uocOption.o optimizer.o $(pkg-config --libs gsl) -o TestUoc
./TestUoc
```

## 🎮 Usage
1. Run the executable.

2. Follow the prompts to:

      * Price a sample option.

      * Generate a pricing grid.

      * Run Nelder-Mead optimization to calibrate parameters to market data.

  3. Watch the output and CSV files rain down.

## 📈 Output
* Terminal output includes intermediate and final option prices.

* A surface.csv file is generated for 3D surface plotting.

## 🧪 Test Cases
Several hard-coded configurations simulate pricing across strike-barrier matrices and known test values. You’ll be prompted to run:

  * Grid search calibration

  * Nelder-Mead optimization (type 1 or 2)

  * Output model and market prices

## ❓ FAQ (aka stuff you may ask eventually)
**Q: Can I use this in production?**

A: Unless your production is a grad school thesis, please don't. This is a teaching engine.

**Q: Is the PIDE method better than Monte Carlo?**

A: Depends. This is deterministic and generally faster for vanilla options. MC is more flexible but slower.

**Q: Why are there so many constructors in uocOption?**

A: Because the past is a foreign country where people do things differently, and constructors were cheap.

## 🧠 Author
Project and code by @dbensik.
