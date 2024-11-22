#ifndef ADD_H
#define ADD_H

#include <SQuIDS/SQuIDS.h>
#include <nuSQuIDS/nuSQuIDS.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix_complex.h>
#include <memory>
#include <vector>
#include <algorithm>

namespace nusquids {

class nuSQUIDS_ADD : public nuSQUIDS {
private:
    unsigned int N_KK;
    unsigned int dim_ADD;
    double a;
    double m0;
    double m1, m2, m3;
    bool NormalOrdering;

    std::unique_ptr<gsl_vector, void(*)(gsl_vector*)> Lambdaq;
    std::unique_ptr<gsl_matrix_complex, void(*)(gsl_matrix_complex*)> W;

    // Helper method for memory allocation
    static gsl_vector* allocate_gsl_vector(size_t size) {
        return gsl_vector_alloc(size);
    }
    static gsl_matrix_complex* allocate_gsl_matrix_complex(size_t rows, size_t cols) {
        return gsl_matrix_complex_alloc(rows, cols);
    }

public:
    nuSQUIDS_ADD(marray<double, 1> E_vector, unsigned int N_KK, double a, double m0,
                 bool NormalOrdering, NeutrinoType NT = both, bool iinteraction = false,
                 std::shared_ptr<CrossSectionLibrary> ncs = nullptr)
        : nuSQUIDS(E_vector, 3 * (N_KK + 1), NT, iinteraction, ncs),
          N_KK(N_KK), dim_ADD(3 * (N_KK + 1)), a(a), m0(m0), NormalOrdering(NormalOrdering),
          Lambdaq(allocate_gsl_vector(dim_ADD - 1), gsl_vector_free),
          W(allocate_gsl_matrix_complex(dim_ADD, dim_ADD), gsl_matrix_complex_free) 
    {
        initializeMassHierarchy();
        initializeMatrices();
    }

    void iniMatrices(gsl_vector*& Lambdaq, gsl_matrix_complex*& W, double th01, double th02, double th12);
    void iniProjectors() override;
    void SetIniFlavorProyectors() override;

private:
    void initializeMassHierarchy() {
        double Deltaq_m21 = 7.65e-05; // Mass difference in eV^2
        double Deltaq_m31 = 0.00247;  // Mass difference in eV^2
        double th01 = 0.563942, th02 = 0.154085, th12 = 0.785398;

        if (NormalOrdering) {
            m1 = m0;
            m2 = std::sqrt(Deltaq_m21 + std::pow(m0, 2));
            m3 = std::sqrt(Deltaq_m31 + std::pow(m0, 2));
        } else {
            m1 = std::sqrt(Deltaq_m31 + std::pow(m0, 2));
            m2 = std::sqrt(Deltaq_m21 + Deltaq_m31 + std::pow(m0, 2));
            m3 = m0;
        }
    }

    void initializeMatrices() {
        // Initialize Lambdaq and W using the mass hierarchy and transformation parameters
        iniMatrices(Lambdaq.get(), W.get(), 0.563942, 0.154085, 0.785398);
    }
};

} // namespace nusquids

#endif // ADD_H