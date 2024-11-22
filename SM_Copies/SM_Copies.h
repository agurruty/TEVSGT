#ifndef SM_COPIES_H
#define SM_COPIES_H

#include <SQuIDS/SQuIDS.h>
#include <nuSQuIDS/nuSQuIDS.h>
#include <memory>
#include <gsl/gsl_blas.h>
#include <complex>

namespace nusquids {

  class nuSQUIDS_SM_Copies : public nuSQUIDS {
  private:
      const int dim = 3; // Dimensions. Fixed to 3. // Used `const` for fixed value
      double N, mu, m0; // SM_Copies parameters
      double m1, m2, m3; // Light neutrino masses
      double mH1, mH2, mH3; // Heavy neutrino masses
      bool NormalOrdering; // Mass hierarchy

      void computeMasses(double Deltaq_m21, double Deltaq_m31); // Encapsulated mass computation for reusability and clarity

  public:
      nuSQUIDS_SM_Copies(
          const marray<double, 1>& E_vector,
          double N, double mu, double m0, bool NormalOrdering,
          unsigned int numneu = 3, NeutrinoType NT = both,
          bool iinteraction = false, 
          std::shared_ptr<CrossSectionLibrary> ncs = nullptr);

      std::unique_ptr<gsl_matrix_complex, void (*)(gsl_matrix_complex*)> GetPMNS(
          double th12 = 0.563942, double th13 = 0.154085, double th23 = 0.785398);

      std::unique_ptr<gsl_matrix_complex, void (*)(gsl_matrix_complex*)> GetTransformationMatrix_SM_Copies();

      void iniProjectors() override;

      void SetIniFlavorProyectors() override;
  };

} // namespace nusquids

#endif // SM_COPIES_H