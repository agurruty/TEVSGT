#include "SM_Copies.h"

namespace nusquids {

  // Helper function to compute neutrino masses
  void nuSQUIDS_SM_Copies::computeMasses(double Deltaq_m21, double Deltaq_m31) {
      if (NormalOrdering) {
          m1 = m0;
          m2 = std::sqrt(Deltaq_m21 + m0 * m0);
          m3 = std::sqrt(Deltaq_m31 + m0 * m0);
      } else {
          m3 = m0;
          m1 = std::sqrt(Deltaq_m31 + m0 * m0);
          m2 = std::sqrt(Deltaq_m21 + Deltaq_m31 + m0 * m0);
      }
      mH1 = mu * m1;
      mH2 = mu * m2;
      mH3 = mu * m3;
  }

  nuSQUIDS_SM_Copies::nuSQUIDS_SM_Copies(
      const marray<double, 1>& E_vector, double N, double mu, double m0, bool NormalOrdering,
      unsigned int numneu, NeutrinoType NT, bool iinteraction, std::shared_ptr<CrossSectionLibrary> ncs)
      : nuSQUIDS(E_vector, numneu, NT, iinteraction, ncs),
        N(N), mu(mu), m0(m0), NormalOrdering(NormalOrdering) {
      
      // Mixing angles and squared mass differences
      constexpr double Deltaq_m21 = 7.65e-05; // Used `constexpr` for compile-time constants
      constexpr double Deltaq_m31 = 0.00247;
      constexpr double th01 = 0.563942, th02 = 0.154085, th12 = 0.785398;

      // Set mixing parameters
      Set_MixingAngle(0, 1, th01);
      Set_MixingAngle(0, 2, th02);
      Set_MixingAngle(1, 2, th12);
      Set_SquareMassDifference(1, Deltaq_m21);
      Set_SquareMassDifference(2, Deltaq_m31);

      computeMasses(Deltaq_m21, Deltaq_m31); // Encapsulated mass computations for clarity

      Set_SquareMassDifference(3, mH1 * mH1 - m0 * m0);
      Set_SquareMassDifference(4, mH2 * mH2 - m0 * m0);
      Set_SquareMassDifference(5, mH3 * mH3 - m0 * m0);
  }

  // Returns the PMNS matrix
  std::unique_ptr<gsl_matrix_complex, void (*)(gsl_matrix_complex*)> nuSQUIDS_SM_Copies::GetPMNS(
      double th12, double th13, double th23) {

      auto U = gsl_matrix_complex_alloc(dim, dim);
      gsl_matrix_complex_set_identity(U);

      auto addRotation = [&U](size_t i, size_t j, double theta, double delta) {
          double c = cos(theta);
          gsl_complex cp = gsl_complex_polar(sin(theta), -delta);
          gsl_complex cpc = gsl_complex_conjugate(cp);

          gsl_matrix_complex_set(U, i, i, gsl_complex_rect(c, 0));
          gsl_matrix_complex_set(U, i, j, cp);
          gsl_matrix_complex_set(U, j, i, cpc);
          gsl_matrix_complex_set(U, j, j, gsl_complex_rect(c, 0));
      };

      addRotation(0, 1, th12, 0); // Used lambda to reduce redundancy in rotation matrix construction
      addRotation(0, 2, th13, 0);
      addRotation(1, 2, th23, 0);

      return std::unique_ptr<gsl_matrix_complex, void (*)(gsl_matrix_complex*)>(U, gsl_matrix_complex_free);
  }

  // Returns the transformation matrix for SM Copies
  std::unique_ptr<gsl_matrix_complex, void (*)(gsl_matrix_complex*)> nuSQUIDS_SM_Copies::GetTransformationMatrix_SM_Copies() {
      auto PMNS = GetPMNS(); // Reused PMNS matrix directly
      auto U = gsl_matrix_complex_alloc(6, 6);
      gsl_matrix_complex_set_zero(U);

      gsl_complex cos_gsl = gsl_complex_rect(std::sqrt((N - 1.0) / N), 0); // Precomputed cos and sin values
      gsl_complex sin_gsl = gsl_complex_rect(1.0 / std::sqrt(N), 0);

      for (unsigned int n = 0; n < 2; ++n) {
          for (unsigned int i = 0; i < dim; ++i) {
              for (unsigned int j = 0; j < dim; ++j) {
                  gsl_complex element1 = gsl_complex_mul(cos_gsl, gsl_matrix_complex_get(PMNS.get(), i, j));
                  gsl_complex element2 = gsl_complex_mul(sin_gsl, gsl_matrix_complex_get(PMNS.get(), i, j));
                  element2 = gsl_complex_mul(gsl_complex_exp(gsl_complex_rect((n + 1) * M_PI, 0)), element2);

                  gsl_matrix_complex_set(U, i + 3 * n, j + 3 * n, element1);
                  gsl_matrix_complex_set(U, i + 3 - 3 * n, j + 3 * n, element2);
              }
          }
      }

      return std::unique_ptr<gsl_matrix_complex, void (*)(gsl_matrix_complex*)>(U, gsl_matrix_complex_free); // Used unique_ptr to manage U memory
  }

  void nuSQUIDS_SM_Copies::iniProjectors() {
      auto U_SM_Copies = GetTransformationMatrix_SM_Copies();

      b0_proj.resize(numneu);
      b1_proj.resize(nrhos, std::vector<squids::SU_vector>(numneu));

      for (unsigned int flv = 0; flv < numneu; ++flv) {
          b0_proj[flv] = squids::SU_vector::Projector(nsun, flv); // Replaced redundant initialization with direct method calls
      }

      for (unsigned int rho = 0; rho < nrhos; ++rho) {
          for (unsigned int flv = 0; flv < numneu; ++flv) {
              b1_proj[rho][flv] = squids::SU_vector::Projector(nsun, flv).Rotate(U_SM_Copies.get()); // Avoided separate rotate call
              AntineutrinoCPFix(rho); // Encapsulated anti-neutrino fix logic
          }
      }
  }

  void nuSQUIDS_SM_Copies::SetIniFlavorProyectors() {
      auto U_SM_Copies = GetTransformationMatrix_SM_Copies();

      for (unsigned int rho = 0; rho < nrhos; ++rho) {
          for (unsigned int flv = 0; flv < numneu; ++flv) {
              for (unsigned int e1 = 0; e1 < ne; ++e1) {
                  evol_b1_proj[rho][flv][e1] = b0_proj[flv].Rotate(U_SM_Copies.get()); // Combined initialization and rotation
                  AntineutrinoCPFix(rho);
              }
          }
      }
  }

} // namespace nusquids