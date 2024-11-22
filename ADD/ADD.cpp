#include <iostream>
#include <complex>
#include <memory>
#include <vector>
#include <algorithm>
#include <cmath>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

namespace nusquids {

// Utility for GSL complex numbers from std::complex
gsl_complex to_gsl(const std::complex<double>& c) {
    return gsl_complex_rect(c.real(), c.imag());
}

// Returns the PMNS matrix.
std::unique_ptr<gsl_matrix_complex, void (*)(gsl_matrix_complex*)> GetPMNS(double th12, double th13, double th23) {
    const size_t dim_SM = 3;

    // Mixing angles matrix with initialization
    gsl_matrix* th = gsl_matrix_calloc(dim_SM, dim_SM);  // calloc initializes to zero
    gsl_matrix_set(th, 0, 1, th12);
    gsl_matrix_set(th, 0, 2, th13);
    gsl_matrix_set(th, 1, 2, th23);

    // Allocate complex matrices
    gsl_matrix_complex* U = gsl_matrix_complex_alloc(dim_SM, dim_SM);
    gsl_matrix_complex_set_identity(U);

    gsl_matrix_complex* R = gsl_matrix_complex_alloc(dim_SM, dim_SM);
    gsl_matrix_complex_set_identity(R);

    gsl_matrix_complex* temp = gsl_matrix_complex_alloc(dim_SM, dim_SM);

    const gsl_complex unit = gsl_complex_rect(1.0, 0.0);
    const gsl_complex zero = gsl_complex_rect(0.0, 0.0);

    // Construct the PMNS matrix
    for (size_t j = 1; j < dim_SM; ++j) {
        for (size_t i = 0; i < j; ++i) {
            double theta = gsl_matrix_get(th, i, j);
            double c = std::cos(theta);
            double s = std::sin(theta);
            gsl_complex cp = gsl_complex_rect(s, 0);

            gsl_matrix_complex_set(R, i, i, gsl_complex_rect(c, 0));
            gsl_matrix_complex_set(R, i, j, cp);
            gsl_matrix_complex_set(R, j, i, gsl_complex_rect(-s, 0));
            gsl_matrix_complex_set(R, j, j, gsl_complex_rect(c, 0));

            // Multiply R with U
            gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, unit, R, U, zero, temp);
            gsl_matrix_complex_memcpy(U, temp);

            // Reset R
            gsl_matrix_complex_set_identity(R);
        }
    }

    // Free intermediate resources
    gsl_matrix_free(th);
    gsl_matrix_complex_free(R);
    gsl_matrix_complex_free(temp);

    return std::unique_ptr<gsl_matrix_complex, void (*)(gsl_matrix_complex*)>(U, gsl_matrix_complex_free);
}

// Optimized iniMatrices function
void iniMatrices(gsl_vector*& Lambdaq, gsl_matrix_complex*& W, double th01, double th02, double th12) {
    // Generate PMNS matrix
    auto PMNS = GetPMNS(th01, th02, th12);

    // Diagonal mass matrix
    gsl_matrix_complex* mass_diag = gsl_matrix_complex_alloc(3, 3);
    gsl_matrix_complex_set_zero(mass_diag);
    gsl_matrix_complex_set(mass_diag, 0, 0, gsl_complex_rect(1.0, 0));  // Replace 1.0 with m1
    gsl_matrix_complex_set(mass_diag, 1, 1, gsl_complex_rect(2.0, 0));  // Replace 2.0 with m2
    gsl_matrix_complex_set(mass_diag, 2, 2, gsl_complex_rect(3.0, 0));  // Replace 3.0 with m3

    gsl_matrix_complex* temp1 = gsl_matrix_complex_alloc(3, 3);
    gsl_matrix_complex* MD = gsl_matrix_complex_alloc(3, 3);

    const gsl_complex unit = gsl_complex_rect(1.0, 0.0);
    const gsl_complex zero = gsl_complex_rect(0.0, 0.0);

    // PMNS * mass_diag
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, unit, PMNS.get(), mass_diag, zero, temp1);
    // (PMNS * mass_diag) * PMNS^T
    gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, unit, temp1, PMNS.get(), zero, MD);

    // Cleanup intermediate matrices
    gsl_matrix_complex_free(mass_diag);
    gsl_matrix_complex_free(temp1);

    // Eigenvalue decomposition
    gsl_eigen_hermv_workspace* workspace = gsl_eigen_hermv_alloc(3);
    gsl_vector* eigenvalues = gsl_vector_alloc(3);
    gsl_eigen_hermv(MD, eigenvalues, W, workspace);

    // Sort eigenvalues and eigenvectors
    std::vector<std::pair<double, size_t>> eigen_pairs(3);
    for (size_t i = 0; i < 3; ++i) {
        eigen_pairs[i] = {gsl_vector_get(eigenvalues, i), i};
    }
    std::sort(eigen_pairs.begin(), eigen_pairs.end());

    gsl_matrix_complex* W_sorted = gsl_matrix_complex_alloc(3, 3);
    for (size_t i = 0; i < 3; ++i) {
        size_t idx = eigen_pairs[i].second;
        gsl_vector_set(eigenvalues, i, eigen_pairs[i].first);
        for (size_t j = 0; j < 3; ++j) {
            gsl_matrix_complex_set(W_sorted, j, i, gsl_matrix_complex_get(W, j, idx));
        }
    }
    gsl_matrix_complex_memcpy(W, W_sorted);

    gsl_vector_free(eigenvalues);
    gsl_matrix_complex_free(W_sorted);
    gsl_matrix_complex_free(MD);
    gsl_eigen_hermv_free(workspace);
}

}  // namespace nusquids
