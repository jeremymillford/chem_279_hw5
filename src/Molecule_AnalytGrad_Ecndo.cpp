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
#include <algorithm> // For std::count

// Compute the X matrix for electron overlap contributions to the gradient
void Molecule::calculate_x_matrix() {
    for (int mu = 0; mu < N_; ++mu) {
        for (int nu = 0; nu < N_; ++nu) {
            if (mu != nu){
            CNDO2Parameters params_mu = semiEmpiricalParams_.at(basisFunctions_[mu].element);
            CNDO2Parameters params_nu = semiEmpiricalParams_.at(basisFunctions_[nu].element);
            double beta_sum = (params_mu.beta + params_nu.beta);

            x_matrix_(mu, nu) = 0.5 * (pA_matrix_(mu, nu) + pB_matrix_(mu, nu)) * beta_sum;
            }
        }
    }
}

// Compute the Y matrix for electron repulsion contributions to the gradient
void Molecule::calculate_y_matrix() {
    // Loop over all atom pairs (A, B)
    for (int A = 0; A < num_atoms_; ++A) {
        for (int B = 0; B < num_atoms_; ++B) {
            if (B != A){
                double term1 = 0.0; // First summation
                double term2 = 0.0; // Second summation
                double term3 = 0.0; // Third summation
                double term4 = 0.0; // Fourth summation

                // Compute P_AA and P_BB
                double P_AA = 0.0, P_BB = 0.0;

                // Sum diagonal elements for P_AA
                for (int mu = 0; mu < N_; ++mu) {
                    if (basisFunctions_[mu].atomic_pos == A) {
                        P_AA += pA_matrix_(mu, mu) + pB_matrix_(mu, mu);
                    }
                }

                // Sum diagonal elements for P_BB
                for (int nu = 0; nu < N_; ++nu) {
                    if (basisFunctions_[nu].atomic_pos == B) {
                        P_BB += pA_matrix_(nu, nu) + pB_matrix_(nu, nu);
                    }
                }

                // Loop over all basis functions on atoms A and B
                for (int mu = 0; mu < N_; ++mu) {
                    if (basisFunctions_[mu].atomic_pos == A) {
                        double P_alpha_mu_mu = pA_matrix_(mu, mu); // Alpha density diagonal element
                        double P_beta_mu_mu = pB_matrix_(mu, mu);  // Beta density diagonal element

                        for (int nu = 0; nu < N_; ++nu) {
                            if (basisFunctions_[nu].atomic_pos == B) {
                                double P_alpha_mu_nu = pA_matrix_(mu, nu); // Alpha off-diagonal
                                double P_beta_mu_nu = pB_matrix_(mu, nu);  // Beta off-diagonal

                                // First summation term
                                term1 += -0.5 * (std::pow(P_alpha_mu_nu, 2) + std::pow(P_beta_mu_nu, 2));

                                // Second summation term (symmetric to first)
                                term2 += -0.5 * (std::pow(P_alpha_mu_nu, 2) + std::pow(P_beta_mu_nu, 2));
                            }
                        }

                        // Third summation term
                        term3 += 0.5 * (P_alpha_mu_mu + P_beta_mu_mu) *
                                    (P_BB - 2.0 * atoms_[B].getValElectrons());
                    }

                    if (basisFunctions_[mu].atomic_pos == B) {
                        double P_alpha_mu_mu = pA_matrix_(mu, mu); // Alpha density diagonal element
                        double P_beta_mu_mu = pB_matrix_(mu, mu);  // Beta density diagonal element

                        // Fourth summation term
                        term4 += 0.5 * (P_alpha_mu_mu + P_beta_mu_mu) *
                                    (P_AA - 2.0 * atoms_[A].getValElectrons());
                    }
                }
            // Combine all terms to calculate y_{AB}
            y_matrix_(A, B) = term1 + term2 + term3 + term4;
            }
        }
    }
}


double Molecule::derivative_primitive_overlap(
    double alpha_k,
    double alpha_l,
    const arma::vec& R_A,
    const arma::vec& R_B,
    const arma::vec& lA,
    const arma::vec& lB,
    int direction) const
{
    // 1D overlaps
    double Sx = Overlap_onedim(alpha_k, alpha_l, R_A(0), R_B(0), lA(0), lB(0));
    double Sy = Overlap_onedim(alpha_k, alpha_l, R_A(1), R_B(1), lA(1), lB(1));
    double Sz = Overlap_onedim(alpha_k, alpha_l, R_A(2), R_B(2), lA(2), lB(2));
    double S_prim = Sx * Sy * Sz;

    // Gaussian exponential derivative term (product rule)
    double coef = -2.0 * alpha_k * alpha_l / (alpha_k + alpha_l);
    double d_exp = coef * (R_A(direction) - R_B(direction)) * S_prim;

    // Return only the exponential derivative — angular effects are in Overlap_onedim
    return d_exp;
}




// Compute the derivative of the *contracted* overlap between basis functions μ and ν
double Molecule::derivative_contracted_overlap(int mu, int nu, int direction) const {
    double dS_mu_nu = 0.0;

    // Loop over primitives k in basisFunctions_[mu]
    size_t nk = basisFunctions_[mu].alphas.size();
    for (size_t k = 0; k < nk; ++k) {
        double alpha_mu_k = basisFunctions_[mu].alphas[k];
        double N_mu_k     = basisFunctions_[mu].norms[k];
        double d_mu_k     = basisFunctions_[mu].coeffs[k];

        // Loop over primitives l in basisFunctions_[nu]
        size_t nl = basisFunctions_[nu].alphas.size();
        for (size_t l = 0; l < nl; ++l) {
            double alpha_nu_l = basisFunctions_[nu].alphas[l];
            double N_nu_l     = basisFunctions_[nu].norms[l];
            double d_nu_l     = basisFunctions_[nu].coeffs[l];

            // Accumulate weighted primitive‐overlap derivative
            dS_mu_nu += d_mu_k * d_nu_l * N_mu_k * N_nu_l
                      * derivative_primitive_overlap(
                            alpha_mu_k, alpha_nu_l,
                            basisFunctions_[mu].center,
                            basisFunctions_[nu].center,
                            basisFunctions_[mu].ang_momenta,
                            basisFunctions_[nu].ang_momenta,
                            direction
                        );
        }
    }

    return dS_mu_nu;
}



// Function to calculate the derivative of the overlap matrix
void Molecule::calculate_overlap_derivatives_RA() {

    // Initialize the derivative matrix
    S_derivative_matrix_RA_.set_size(3, N_ * N_);
    S_derivative_matrix_RA_.fill(0.0);


    // Loop over all pairs of basis functions (mu, nu) with mu <= nu to avoid double counting
    for (int mu = 0; mu < N_; ++mu) {
        for (int nu = mu; nu < N_; ++nu) { // Symmetry: only compute for mu <= nu
            // Get the atom
            int center_mu = basisFunctions_[mu].atomic_pos;
            int center_nu = basisFunctions_[nu].atomic_pos;
            // Loop over all nuclei A
            for (int A = 0; A < num_atoms_; ++A) {
                // Check if derivative is non-zero
                if (!((center_mu == A && center_nu != A) || (center_nu == A && center_mu != A))) {
                    continue; // Skip if derivative is zero
                }

                // Flatten the (mu, nu) pair into a column index
                int index = mu * N_ + nu;

                // Calculate the derivatives in x, y, z directions
                for (int direction = 0; direction < 3; ++direction) {
                    double dS_mu_nu = derivative_contracted_overlap(mu, nu, direction);

                    // Add the derivative to the matrix (account for symmetry)
                    S_derivative_matrix_RA_(direction, index) = dS_mu_nu;

                    // If mu != nu, include the symmetric contribution
                    if (mu != nu) {
                        int sym_index = nu * N_ + mu;
                        S_derivative_matrix_RA_(direction, sym_index) = -dS_mu_nu;
                    }
                }
            }
        }
    }
}

// Compute the derivative of the two-electron integral with respect to atomic coordinates
double Molecule::I2e_pG_wrtRA(arma::vec &Ra, arma::vec &Rb, double sigmaA, double sigmaB, int coordIndex) const {
    double conversion_factor = 27.211324570273; // Conversion value a.u. to eV

    // Compute necessary terms
    double V2 = 1.0 / (sigmaA + sigmaB);         // V^2
    double V = std::sqrt(V2);                    // V = sqrt(V^2)
    arma::vec Rdiff = Ra - Rb;                   // Vector difference R_A - R_B
    double Rd = arma::norm(Rdiff, 2);            // Magnitude |R_A - R_B|
    double Rdsq = arma::dot(Rdiff, Rdiff);       // Squared magnitude |R_A - R_B|^2
    double T = V2 * Rdsq;                        // T = V^2 * |R_A - R_B|^2

    // Handle edge case: Avoid division by zero if Rd == 0
    if (Rd == 0) {
        return 0.0; // Special case: no contribution if R_A = R_B
    }

    // Compute normalization constants
    double UA = std::pow(M_PI * sigmaA, 1.5);    // U_A = (π * σ_A)^(3/2)
    double UB = std::pow(M_PI * sigmaB, 1.5);    // U_B = (π * σ_B)^(3/2)

    // Compute derivatives of overlap components
    double partial_term = Rdiff[coordIndex] / Rdsq;  // Select x, y, or z component of R_A - R_B
    double scalar_part1 = -(std::erf(std::sqrt(T)) / Rd); // Scalar part 1: -erf(sqrt(T)) / |R_A - R_B|
    double scalar_part2 = ((2 * V) / std::sqrt(M_PI)) * std::exp(-T); // Scalar part 2: (2V /sqrt(pi)) * exp(-T)
    double scalar_derivative = scalar_part1 + scalar_part2;

    // Derivative terms for the overlap
    double prefactor = (UA * UB) * partial_term; // Prefactor for the selected coordinate
    double derivative_result = prefactor * scalar_derivative;

    // Return result scaled to eV
    return derivative_result * conversion_factor;
}

// Function to evaluate 2-electron integral for atomic orbitals
double Molecule::eval_gamma_wrtRA(const BasisFunction& basisFuncA, const BasisFunction& basisFuncB, int coordIndex) const {
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


    double gamma_wrtRA = 0.0;
    // Loop over k, k', l, l' as in Equation 3.3
    for (size_t k1 = 0; k1 < len; ++k1){
        for (size_t k2 = 0; k2 < len; ++k2) {
            double sigmaA = 1.0 / (alpha_a[k1] + alpha_a[k2]); // Equation 3.10
            for (size_t l1 = 0; l1 < len; ++l1) {
                for (size_t l2 = 0; l2 < len; ++l2) {
                    double sigmaB = 1.0 / (alpha_b[l1] + alpha_b[l2]); // Equation 3.10
                    double dI2e = I2e_pG_wrtRA(Ra, Rb, sigmaA, sigmaB, coordIndex); // Integral calculation
                    gamma_wrtRA += (da[k1] * norm_a[k1]) * (da[k2] * norm_a[k2]) * 
                            (db[l1] * norm_b[l1]) * (db[l2] * norm_b[l2]) * dI2e; // Equation 3.3
                }
            }
        }
    }
    return gamma_wrtRA;
}

// Compute gamma matrix derivatives for atomic coordinates
void Molecule::calculate_gammaAB_RA() {
    // Loop over atom pairs (A, B)
    for (int A = 0; A < num_atoms_; ++A) {
        for (int B = 0; B < num_atoms_; ++B) {
            if (A != B) {
                // Loop over x, y, z coordinates
                for (int coordIndex = 0; coordIndex < 3; ++coordIndex) {
                    double gamma_value = 0.0;

                    for (int mu = 0; mu < N_; ++mu) {
                        // only s‐orbitals (orbitalType == 0) on atom A
                        if (basisFunctions_[mu].atomic_pos != A || basisFunctions_[mu].orbitalType != 0) 
                            continue;
                        for (int nu = 0; nu < N_; ++nu) {
                            // only s‐orbitals on atom B
                            if (basisFunctions_[nu].atomic_pos != B || basisFunctions_[nu].orbitalType != 0) 
                                continue;
                            // now this is the single s–s pair
                            gamma_value += eval_gamma_wrtRA(
                                basisFunctions_[mu],
                                basisFunctions_[nu],
                                coordIndex
                            );
                        }
                    }
                    

                    // Store it *with* the missing 1/2 factor
                    int columnIndex = A * num_atoms_ + B;
                    gammaAB_RA_(coordIndex, columnIndex) =  gamma_value;
                }
            }
        }
    }
}


// Compute the gradient of the nuclear repulsion term for a specific atom
arma::vec Molecule::compute_vnuc_derivative(int A) {
    double conversion_factor = 27.211324570273; // Conversion value a.u. to eV
    arma::vec vnuc_derivative(3, arma::fill::zeros);
    arma::vec RA = atoms_[A].getCenter();
    double ZA = atoms_[A].getValElectrons(); // Nuclear charge of atom A

    for (int B = 0; B < atoms_.size(); ++B) {
        if (B != A) {
            arma::vec RB = atoms_[B].getCenter();
            double ZB = atoms_[B].getValElectrons(); // Nuclear charge of atom B
            arma::vec RAB = RA - RB;
            double RAB_norm = arma::norm(RAB, 2);

            // Derivative of nuclear repulsion term
            vnuc_derivative += -ZA * ZB * RAB / std::pow(RAB_norm, 3);
        }
    }

    return vnuc_derivative * conversion_factor; // Conversion value a.u. to eV
}

// Compute the CNDO/2 energy gradient with respect to atomic coordinates
void Molecule::compute_cndo2_gradient() {
    // Initialize gradients
    arma::mat gradient(3, num_atoms_, arma::fill::zeros);
    arma::mat gradient_nuclear(3, num_atoms_, arma::fill::zeros);
    arma::mat gradient_electron(3, num_atoms_, arma::fill::zeros);

    // Compute required matrices
    calculate_x_matrix();
    calculate_y_matrix();
    calculate_overlap_derivatives_RA();
    calculate_gammaAB_RA();

    // Loop over coordinates (x, y, z)
    for (int coord = 0; coord < 3; ++coord) {
        // Iterate over all atom pairs (A, B)
        for (int A = 0; A < num_atoms_; ++A) {
            for (int B = A + 1; B < num_atoms_; ++B) { // Ensure each pair is processed only once
                double pair_contribution = 0.0;

                // Compute overlap contribution for this pair
                for (int mu = 0; mu < N_; ++mu) {
                    for (int nu = mu + 1; nu < N_; ++nu) { // mu < nu to avoid double counting
                        int columnIndex = mu * N_ + nu;
                        pair_contribution += x_matrix_(mu, nu) * S_derivative_matrix_RA_(coord, columnIndex);
                    }
                }

                // Distribute contributions symmetrically
                gradient_electron(coord, A) += pair_contribution;
                gradient_electron(coord, B) -= pair_contribution; // Opposite sign for symmetry
            }
        }
    }

    for (int coord = 0; coord < 3; ++coord) {
        // Contribution from gamma derivatives (y_AB * gamma_AB^(R_A))
        for (int A = 0; A < num_atoms_; ++A) {
            double sum_gamma = 0.0;
            for (int B = 0; B < num_atoms_; ++B) {
                if (B != A) {
                    int columnIndex = A * num_atoms_ + B;
                    sum_gamma += y_matrix_(A, B) * gammaAB_RA_(coord, columnIndex);
                }
            }
            gradient_electron(coord, A) += sum_gamma;
        }

        // Contribution from nuclear derivatives (V_nuc^(R_A))
        for (int A = 0; A < num_atoms_; ++A) {
            arma::vec vnuc = compute_vnuc_derivative(A);
            gradient_nuclear(coord, A) = vnuc[coord];
        }
    }

    // Total gradient: Add nuclear and electron contributions
    gradient = gradient_nuclear + gradient_electron;
    

    // Output
    std::cout << "S_uv_RA (3 x N*N):\n"<< S_derivative_matrix_RA_ << std::endl;
    std::cout << "gammaAB_RA:\n" << gammaAB_RA_ << std::endl;
    std::cout << "Gradient (Nuclear part):\n" << gradient_nuclear << std::endl;
    std::cout << "Gradient (Electron part - Overlap + Gamma):\n" << gradient_electron << std::endl;
    std::cout << "Gradient (Total):\n" << gradient << std::endl;

    
    print_Xuv_yAB(); // Debug print for x_uv and y_AB matrices
}