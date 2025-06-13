
#include "GslTridiagonalSolver.hpp"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <stdexcept>

std::vector<double> GslTridiagonalSolver::solve(
    const std::vector<double>& a,
    const std::vector<double>& b,
    const std::vector<double>& c,
    const std::vector<double>& d
) {
    int n = b.size();

    if (a.size() != n - 1 || c.size() != n - 1 || d.size() != n) {
        throw std::invalid_argument("Invalid input sizes for GSL tridiagonal solver.");
    }

    gsl_vector* gsl_a = gsl_vector_alloc(n - 1);
    gsl_vector* gsl_b = gsl_vector_alloc(n);
    gsl_vector* gsl_c = gsl_vector_alloc(n - 1);
    gsl_vector* gsl_d = gsl_vector_alloc(n);
    gsl_vector* gsl_x = gsl_vector_alloc(n);

    for (int i = 0; i < n - 1; ++i) {
        gsl_vector_set(gsl_a, i, a[i]);
        gsl_vector_set(gsl_c, i, c[i]);
    }
    for (int i = 0; i < n; ++i) {
        gsl_vector_set(gsl_b, i, b[i]);
        gsl_vector_set(gsl_d, i, d[i]);
    }

    gsl_linalg_solve_tridiag(gsl_b, gsl_c, gsl_a, gsl_d, gsl_x);

    std::vector<double> result(n);
    for (int i = 0; i < n; ++i) {
        result[i] = gsl_vector_get(gsl_x, i);
    }

    gsl_vector_free(gsl_a);
    gsl_vector_free(gsl_b);
    gsl_vector_free(gsl_c);
    gsl_vector_free(gsl_d);
    gsl_vector_free(gsl_x);

    return result;
}
