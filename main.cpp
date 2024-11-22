#define USE_SM
#include <vector>
#include <iostream>
#include <fstream>
#include <nuSQuIDS/nuSQuIDS.h>

#ifdef USE_ADD
#include "ADD/ADD.h"

#elif defined(USE_SM_COPIES)
#include "SM_Copies/SM_Copies.h"

#elif defined(USE_DARKDIM)
#include "DarkDim/DarkDim.h"

#elif defined(USE_SM)
#else 
#error "Error: Model type not defined."
#endif

using namespace nusquids;

// Helper function to handle errors
void HandleError(const std::string& error_message) {
    std::cerr << "Error: " << error_message << std::endl;
    std::exit(1); // Added helper to reduce repetitive error handling
}

int main() {
    squids::Const units;

    // Model-specific variables
    #ifdef USE_ADD
    double a_min, a_max, m0_min, m0_max;
    unsigned int N_a, N_m0, N_KK = 2;
    std::cin >> a_min >> a_max >> N_a >> m0_min >> m0_max >> N_m0;
    unsigned int numneu = 3 * (N_KK + 1);
    auto a_vec = linspace_vec(a_min, a_max, N_a); // Reused linspace_vec directly to minimize redundancy
    auto m0_vec = linspace_vec(m0_min, m0_max, N_m0);

    #elif defined(USE_SM_COPIES)
    double N_min, N_max, mu_min, mu_max, m0_min, m0_max;
    unsigned int N_N, N_mu, N_m0;
    std::cin >> N_min >> N_max >> N_N >> mu_min >> mu_max >> N_mu >> m0_min >> m0_max >> N_m0;
    unsigned int numneu = 6;
    auto N_vec = linspace_vec(N_min, N_max, N_N);
    auto mu_vec = linspace_vec(mu_min, mu_max, N_mu);
    auto m0_vec = linspace_vec(m0_min, m0_max, N_m0);

    #elif defined(USE_DARKDIM)
    double a_min, a_max, m0_min, m0_max, ca1_min, ca1_max, ca2_min, ca2_max, ca3_min, ca3_max;
    unsigned int N_a, N_m0, N_ca1, N_ca2, N_ca3, N_KK = 2;
    std::cin >> a_min >> a_max >> N_a >> m0_min >> m0_max >> N_m0 >> ca1_min >> ca1_max >> N_ca1
             >> ca2_min >> ca2_max >> N_ca2 >> ca3_min >> ca3_max >> N_ca3;
    unsigned int numneu = 3 * (N_KK + 1);
    auto a_vec = linspace_vec(a_min, a_max, N_a);
    auto m0_vec = linspace_vec(m0_min, m0_max, N_m0);
    auto ca1_vec = linspace_vec(ca1_min, ca1_max, N_ca1);
    auto ca2_vec = linspace_vec(ca2_min, ca2_max, N_ca2);
    auto ca3_vec = linspace_vec(ca3_min, ca3_max, N_ca3);

    #elif defined(USE_SM)
    unsigned int numneu = 3;
    #endif

    // Common inputs for all models
    double E_min, E_max, medium_param_min, medium_param_max;
    unsigned int N_medium_param, N_energy_grid = 200, Nen;
    std::string medium, neutrino_type_str, NormalOrdering_str;
    NeutrinoType neutrino_type;
    bool NormalOrdering;

    std::cin >> E_min >> E_max >> medium >> medium_param_min >> medium_param_max >> N_medium_param
             >> neutrino_type_str >> NormalOrdering_str >> N_energy_grid >> Nen;

    // Map neutrino type
    if (neutrino_type_str == "neutrino") neutrino_type = neutrino;
    else if (neutrino_type_str == "antineutrino") neutrino_type = antineutrino;
    else HandleError("neutrino_type must be either 'neutrino' or 'antineutrino'.");

    // Map NormalOrdering
    if (NormalOrdering_str == "true") NormalOrdering = true;
    else if (NormalOrdering_str == "false") NormalOrdering = false;
    else HandleError("NormalOrdering must be either 'true' or 'false'.");

    // Initial flux ratios
    std::vector<double> initial_flux_ratios(3);
    double sum_ratios = 0.0;
    for (double& ratio : initial_flux_ratios) {
        std::cin >> ratio;
        sum_ratios += ratio;
    }
    if (sum_ratios != 1.0) HandleError("The initial flavor flux ratios must add to 1.");

    // Medium parameter vector
    auto medium_param_vec = linspace_vec(medium_param_min, medium_param_max, N_medium_param);

    // File for output
    std::ofstream file("fluxes_flavor.txt");
    if (!file.is_open()) HandleError("Failed to open output file.");

    for (const double medium_param : medium_param_vec) {
        #ifdef USE_ADD
        for (const double m0 : m0_vec) {
            for (const double a : a_vec) {
                nuSQUIDS_ADD nus(logspace(E_min * units.GeV, E_max * units.GeV, N_energy_grid), N_KK, a, m0, NormalOrdering, neutrino_type);
        #elif defined(USE_SM_COPIES)
        for (const double m0 : m0_vec) {
            for (const double mu : mu_vec) {
                for (const double N : N_vec) {
                    nuSQUIDS_SM_Copies nus(logspace(E_min * units.GeV, E_max * units.GeV, N_energy_grid), N, mu, m0, NormalOrdering, numneu, neutrino_type);
        #elif defined(USE_DARKDIM)
        for (const double ca3 : ca3_vec) {
            for (const double ca2 : ca2_vec) {
                for (const double ca1 : ca1_vec) {
                    std::array<double, 3> ca = {ca1, ca2, ca3};
                    for (const double m0 : m0_vec) {
                        for (const double a : a_vec) {
                            nuSQUIDS_DarkDim nus(logspace(E_min * units.GeV, E_max * units.GeV, N_energy_grid), N_KK, a, m0, ca, NormalOrdering, neutrino_type);
        #elif defined(USE_SM)
        nuSQUIDS nus(logspace(E_min * units.GeV, E_max * units.GeV, N_energy_grid), numneu, neutrino_type, false);
        #endif

            // Medium setup
            if (medium == "vacuum") {
                nus.Set_Body(std::make_shared<Vacuum>());
                nus.Set_Track(std::make_shared<Vacuum::Track>(medium_param * units.km));
            } else if (medium == "earth") {
                double phi = std::acos(medium_param);
                auto earth_atm = std::make_shared<EarthAtm>();
                nus.Set_Body(earth_atm);
                nus.Set_Track(std::make_shared<EarthAtm::Track>(earth_atm->MakeTrack(phi)));
            } else {
                HandleError("Invalid medium of propagation.");
            }

            // Numerical settings
            nus.Set_h_max(500.0 * units.km);
            nus.Set_GSL_step(gsl_odeiv2_step_rk4);
            nus.Set_rel_error(1.0e-5);
            nus.Set_abs_error(1.0e-5);

            // Initial state
            auto E_range = nus.GetERange();
            marray<double, 2> inistate{E_range.size(), numneu};
            double N0 = 1.0e18;

            for (size_t i = 0; i < inistate.extent(0); ++i) {
                for (size_t k = 0; k < 3; ++k) inistate[i][k] = initial_flux_ratios[k] * N0 * pow(E_range[i], -2);
                for (size_t k = 3; k < inistate.extent(1); ++k) inistate[i][k] = 0.0;
            }
            nus.Set_initial_state(inistate, flavor);
            nus.EvolveState();

            // Write results
            auto lE_vec = linspace_vec(std::log10(E_min), std::log10(E_max), Nen);
            for (const double lE : lE_vec) {
                double E = pow(10.0, lE);
                file << E;
                for (int fl = 0; fl < 3; ++fl)
                    file << " " << nus.EvalFlavor(fl, E * units.GeV) / (N0 * pow(E * units.GeV, -2));
                file << std::endl;
            }

        #ifdef USE_ADD
            }
        }
        #elif defined(USE_SM_COPIES)
                }
            }
        #elif defined(USE_DARKDIM)
                        }
                    }
                }
            }
        #endif
    }

    file.close();
    std::cout << "Done!" << std::endl;
    return 0;
}
