// Implementation of a Molecule for the hartree-fock CNDO/2 method. 

#include "Atom.h"
#include "Molecule.h"
#include "BasisFunction.h"

#include <cmath>
#include <stdexcept>
#include <cassert>
#include <iostream>
#include <armadillo>

#include <unordered_map>
#include <string>

#include <iostream>
#include <vector>
#include <algorithm> 

// 2-electron integral of two primitive Gaussians
double Molecule::I2e_pG(arma::vec &Ra, arma::vec &Rb, double sigmaA, double sigmaB) const {
    double conversion_factor = 27.211324570273; // Conversion value a.u. to eV

    double V2 = 1.0 / (sigmaA + sigmaB); 
    double Rd = arma::norm(Ra - Rb, 2);
    double T = V2 * std::pow(Rd, 2); 

    double UA = std::pow(M_PI * sigmaA, 1.5); 
    double UB = std::pow(M_PI * sigmaB, 1.5); 

    if (Rd == 0.0) {
        // If RA == RB, use equation 
        double result = UA * UB * 2.0 * std::sqrt(V2 / M_PI);
        return result * conversion_factor; // Converting a.u. to eV consistent with the semi-empirical params
    }

    // Equation 3.14
    double result = UA * UB * std::sqrt(1.0 / (Rd * Rd)) * std::erf(std::sqrt(T));

    return result * conversion_factor; // Converting a.u. to eV consistent with the semi-empirical params
}

// Function to evaluate 2-electron integral for atomic orbitals (s-type)
double Molecule::eval_gamma(const BasisFunction& basisFuncA, const BasisFunction& basisFuncB) const {
    int len = basisFuncA.alphas.size();
    assert(basisFuncB.alphas.size() == len);

    arma::vec Ra = basisFuncA.center;
    arma::vec Rb = basisFuncB.center;
    const std::vector<double>& alpha_a = basisFuncA.alphas;
    const std::vector<double>& da = basisFuncA.coeffs;
    const std::vector<double>& norm_a = basisFuncA.norms;
    const std::vector<double>& alpha_b = basisFuncB.alphas;
    const std::vector<double>& db = basisFuncB.coeffs;
    const std::vector<double>& norm_b = basisFuncB.norms;


    double gamma = 0.0;
    // Loop over k, k', l, l' as in Equation 
    for (size_t k1 = 0; k1 < len; ++k1){
        for (size_t k2 = 0; k2 < len; ++k2) {
            double sigmaA = 1.0 / (alpha_a[k1] + alpha_a[k2]); 
            for (size_t l1 = 0; l1 < len; ++l1) {
                for (size_t l2 = 0; l2 < len; ++l2) {
                    double sigmaB = 1.0 / (alpha_b[l1] + alpha_b[l2]); 
                    double I2e = I2e_pG(Ra, Rb, sigmaA, sigmaB); // Integral calculation
                    gamma += (da[k1] * norm_a[k1]) * (da[k2] * norm_a[k2]) * 
                            (db[l1] * norm_b[l1]) * (db[l2] * norm_b[l2]) * I2e; 
                }
            }
        }
    }
    return gamma;
}

// Returns the diagonal value of the fock matrix for a basis function and current electron density.
double Molecule::diagonal_value(int mu, bool isAlpha) const {
    // Set the density matrix
    const arma::mat& densityMatrix = isAlpha ? pA_matrix_ : pB_matrix_;

    // Get element atomic position
    int element = basisFunctions_[mu].element;
    int A = basisFunctions_[mu].atomic_pos;

    if (semiEmpiricalParams_.find(element) == semiEmpiricalParams_.end()) {
        throw std::invalid_argument("Element not supported in semi-empirical parameter map.");
    }

    // Get semiempirical values
    CNDO2Parameters params = semiEmpiricalParams_.at(element);

    // Determine the correct ionization/affinity energies based on orbital type
    double IA_energy;
    if (basisFunctions_[mu].orbitalType == 0 ) { // S-orbital
        IA_energy = semiEmpiricalParams_.at(element).IA_s;
    } else if (basisFunctions_[mu].orbitalType == 1) { // P-orbital
        IA_energy = semiEmpiricalParams_.at(element).IA_p;
    } else {
        throw std::runtime_error("Unknown orbital type in basis function.");
    }

    // Compute the total electron density minus the valence electrons
    double ptot_Za = pTot_(A) - atoms_[A].getValElectrons();

    // Compute either p^alpha or p^beta - 1/2
    double pA_uu_half = densityMatrix(mu, mu) - 0.5;

    // Calculate the self-repulsion term (AA)
    double gamma_AA = gamma_matrix_(A, A);
    double eRepulsion_AA_exchange = (ptot_Za - pA_uu_half) * gamma_AA;

    // Loop to calculate gamma_AB terms for interactions with other atoms
    double sum_gamma_AB = 0.0;
    for (int B = 0; B < atoms_.size(); ++B) {
        if (B != A) {
            sum_gamma_AB += (pTot_(B) - atoms_[B].getValElectrons()) * gamma_matrix_(A, B);
        }
    }

    // Equation 1.4 to calculate the final diagonal Fock matrix element
    return -IA_energy + eRepulsion_AA_exchange + sum_gamma_AB;
}

// Returns the off-diagonal value of the fock matrix for two basis functions.
double Molecule::offDiagonal_value(int mu, int nu, bool isAlpha) const{
    // Set the density matrix (alpha or beta)
    const arma::mat& densityMatrix = isAlpha ? pA_matrix_ : pB_matrix_;

    // Get AB position and element
    int elementA = basisFunctions_[mu].element;
    int A = basisFunctions_[mu].atomic_pos;
    int elementB = basisFunctions_[nu].element;
    int B = basisFunctions_[nu].atomic_pos;

    // Get semiempirical parameters
    CNDO2Parameters params_A = semiEmpiricalParams_.at(elementA);
    CNDO2Parameters params_B = semiEmpiricalParams_.at(elementB);

    // Get overlap matrix element
    double s_uv = S_overlap_matrix_(mu, nu);

    // Calculate average B
    double beta_sum = params_A.beta + params_B.beta;

    // Calculate p(mu, nu) and gamma_AB
    double p_mu_nu = densityMatrix(mu, nu);
    double gamma_AB = gamma_matrix_(A, B);

    return (0.5 * beta_sum * s_uv) - (p_mu_nu * gamma_AB);
}

// Builds the Fock matrix for the molecule.
arma::mat Molecule::build_fock_matrix(bool isAlpha) {
    arma::mat Fock_matrix;
    Fock_matrix.set_size(N_, N_);

    for (int mu = 0; mu < N_; ++mu) {
        for (int nu = 0; nu < N_; ++nu) {
        
            // Get Overlap Matrix value
            double S_uv = S_overlap_matrix_(mu, nu);

            // Use parameters for diagonal and off-diagonal calculations as needed
            double h = (mu == nu) ? diagonal_value(mu, isAlpha) : offDiagonal_value(mu, nu, isAlpha);

            // Store in the Fock matrix
            Fock_matrix(mu, nu) = h;
        }
    }
    return Fock_matrix;
}

void Molecule::build_Fmatrix_wSCF(bool verbose) {
    bool converged = false;
    int max_iterations = 1000;
    double convergence_threshold = 1e-6;

    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        // Bools for selecting denisty matrices
        bool alpha = true;
        bool beta = false;

        // 1. Build the Fock matrices F_alpha and F_beta based on current density
        FA_matrix_ = build_fock_matrix(alpha);
        FB_matrix_ = build_fock_matrix(beta);

        // 2. Solve F C = S C epsilon for alpha and beta electrons
        arma::eig_sym(epsilon_A_, cA_matrix_, FA_matrix_);
        arma::eig_sym(epsilon_B_, cB_matrix_, FB_matrix_);

        // 3. Update the alpha and beta density matrices
        arma::mat new_pA_matrix = cA_matrix_.cols(0, alpha_e_ - 1) * cA_matrix_.cols(0, alpha_e_ - 1).t();
        arma::mat new_pB_matrix = cB_matrix_.cols(0, beta_e_ - 1) * cB_matrix_.cols(0, beta_e_ - 1).t();

        // 4. Check for convergence
        double diff_alpha = arma::accu(arma::abs(new_pA_matrix - pA_matrix_));
        double diff_beta = arma::accu(arma::abs(new_pB_matrix - pB_matrix_));

        if (diff_alpha < convergence_threshold && diff_beta < convergence_threshold) {
            converged = true;
        }

        // 5. Update the density matrices for the next iteration
        pA_matrix_ = new_pA_matrix;
        pB_matrix_ = new_pB_matrix;
        calculate_pTot();

        // if verbose, print the current state
        if (verbose) printSCF(converged, iteration);

        // if converged, break the loop
        if (converged) break;
    }   

    if (converged) {
        std::cout << "SCF converged successfully." << std::endl;

    } else {
        std::cerr << "SCF did not converge within the maximum number of iterations." << std::endl;
    }
}

// Function to calculate the gamma matrix for the molecule
void Molecule::calculate_gamma() {
    int int_gamma_size = atoms_.size();
    gamma_matrix_.set_size(int_gamma_size, int_gamma_size);

    // Loop over each pair of basis functions to calculate gamma
    for (int A = 0; A < int_gamma_size; ++A) {
        for (int B = 0; B < int_gamma_size; ++B) {
            gamma_matrix_(A, B) = eval_gamma(basisFunctions_[A], basisFunctions_[B]);
        }
    }
}


// Function to build the core Hamiltonian matrix (H_core)
void Molecule::build_core_hamiltonian() {
    H_core_matrix_.set_size(N_, N_);

    for (int mu = 0; mu < N_; ++mu) {
        for (int nu = 0; nu < N_; ++nu) {
            if (mu == nu) {
                // Diagonal elements (Equation 2.6)
                int A = basisFunctions_[mu].atomic_pos;
                // Get Atomic element
                int element_A = atoms_[A].getAtomicNumber();
            // std::cout << "Element_mu";
                int Z_A = atoms_[A].getValElectrons();
                double gamma_AA = gamma_matrix_(A, A);

                // Determine the correct ionization/affinity energies based on orbital type
                double IA_energy;
                if (basisFunctions_[mu].orbitalType == 0 ) { // S-orbital
                    IA_energy = semiEmpiricalParams_.at(element_A).IA_s;
                } else if (basisFunctions_[mu].orbitalType == 1) { // P-orbital
                    IA_energy = semiEmpiricalParams_.at(element_A).IA_p;
                } else {
                    throw std::runtime_error("Unknown orbital type in basis function.");
                }
            // std::cout << " IA_energy : " << IA_energy << std::endl;

                double h_mu_mu = - IA_energy - ((Z_A - 0.5) * gamma_AA);
                double sum = 0;
                for (int B = 0; B < atoms_.size(); ++B) {
                    if (B != A) {
                        int Z_B = atoms_[B].getValElectrons();
                        sum += Z_B * gamma_matrix_(A, B);
                    }
                }
                H_core_matrix_(mu, nu) = h_mu_mu - sum;
            } else {
                // Off-diagonal elements (Equation 2.7)
            // std::cout << "(mu,nu) : (" << mu << "," << nu << ")" << std::endl;
                CNDO2Parameters params_mu = semiEmpiricalParams_.at(basisFunctions_[mu].element);
                CNDO2Parameters params_nu = semiEmpiricalParams_.at(basisFunctions_[nu].element);
            // std::cout << "Element_mu_beta = " << params_mu.beta << " & Element_nu_beta = " << params_nu.beta << std::endl;

                double beta_sum = (params_mu.beta + params_nu.beta);
                double s_uv = S_overlap_matrix_(mu, nu);
                H_core_matrix_(mu, nu) = 0.5 * beta_sum * s_uv;
            }
        }
    }
}

// Function to calculate the atomic electron populations from the basis function density matrix
void Molecule::calculate_pTot() {
    int num_atoms = atoms_.size();
    pTot_.set_size(num_atoms);
    pTot_.zeros(); // Set all elements to zero

    // Sum the densities from pTot_matrix_ for each atom
    for (int mu = 0; mu < N_; ++mu) {
        int atom_mu = basisFunctions_[mu].atomic_pos;  // Get the atom index for this basis function
        // Add the diagonal elements from the total density matrix for basis functions centered on atom_mu
        pTot_(atom_mu) += pA_matrix_(mu, mu) + pB_matrix_(mu, mu);
    }
}


// Calculates the total electronic energy of the molecule.
double Molecule::calc_totalE() {
    double conversion_factor = 27.211324570273; // a.u. to eV

    // Alpha and beta electron energies
    alpha_electron_E_ = 0.5 * arma::accu(pA_matrix_ % (H_core_matrix_ + FA_matrix_));
    beta_electron_E_ = 0.5 * arma::accu(pB_matrix_ % (H_core_matrix_ + FB_matrix_));

    // Nuclear repulsion energy
    for (int A = 0; A < atoms_.size(); ++A) {
        for (int B = A + 1; B < atoms_.size(); ++B) {
            int Z_A = atoms_[A].getValElectrons();
            int Z_B = atoms_[B].getValElectrons();
            arma::vec RA = arma::conv_to<arma::vec>::from(atoms_[A].getCenter());
            arma::vec RB = arma::conv_to<arma::vec>::from(atoms_[B].getCenter());
            double R_AB = arma::norm(RA - RB);
            nuclear_repulsion_E_ += (Z_A * Z_B) / R_AB ;
        }
    }
    // Convert to eV
    nuclear_repulsion_E_ = nuclear_repulsion_E_ * conversion_factor;

    // Total energy
    total_energy_ = alpha_electron_E_ + beta_electron_E_ + nuclear_repulsion_E_;
    return total_energy_;
}

std::string getMoleculeName(const std::string& filePath){
    if (filePath.find("N.txt") != std::string::npos) return "N";
    if (filePath.find("N2.txt") != std::string::npos) return "N2";
    if (filePath.find("O.txt") != std::string::npos) return "O";
    if (filePath.find("O2.txt") != std::string::npos) return "O2";
    return "NULL";
}

// Function to get bond Energy
void calculateBondEnergy(std::unordered_map<std::string, double>& energies) {

    // Calculate bond energy for N2
    if (energies.find("N2") != energies.end() && energies.find("N") != energies.end()) {
        double bond_energy = (2 * energies["N"]) - energies["N2"];
        std::cout << "The bond energy of N2 is " << bond_energy << " eV." << std::endl;
    }

    // Calculate bond energy for O2
    if (energies.find("O2") != energies.end() && energies.find("O") != energies.end()) {
        double bond_energy = (2 * energies["O"]) - energies["O2"];
        std::cout << "The bond energy of O2 is " << bond_energy << " eV." << std::endl;
   
    }
}


