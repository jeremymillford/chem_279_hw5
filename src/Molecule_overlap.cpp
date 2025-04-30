//implementation for the overlap matrix S in Molecule class

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

// Function to calculate the exponential prefactor in Eq. 2.9
double Molecule::exponential_prefactor(double alpha, double beta, double X_A, double X_B) const {
    double exponent = (-alpha * beta * std::pow((X_A - X_B), 2)) / (alpha + beta);

    double root_term = std::sqrt(M_PI / (alpha + beta));

    return std::exp(exponent) * root_term;
}

// Function to calculate factorial (n!)
int Molecule::factorial(int n) const {
    if (n < 0) {
        throw std::invalid_argument("Factorial is not defined for negative numbers.");
    }
    int fact = 1;
    for (int i = 2; i <= n; ++i) {
        fact *= i;
    }
    return fact;
}

// Function to calculate double factorial (n!!)
int Molecule::double_factorial(int n) const {
    if (n < -1) {
        throw std::invalid_argument("Double factorial is not defined for negative numbers.");
    }
    if (n == 0 || n == -1) return 1;  // Handle 0!! = 1 and -1!! = 1 by definition
    int fact = 1;
    for (int i = n; i > 0; i -= 2) {
        fact *= i;
    }
    return fact;
}

// Function to calculate binomial coefficient
double Molecule::calculate_binomial(int l, int i) const {
    if(l < i || i < 0){
        throw std::invalid_argument("Comination: the number of elements should be bigger than selection numbers AND two numbers should be positive\n");
        return EXIT_FAILURE;
    }

    return factorial(l) / (factorial(i) * factorial(l - i));
}

// Function to calculate the center of the product Gaussian
double Molecule::calculate_RP(double alpha, double beta, double X_A, double X_B) const {
    return ((alpha * X_A) + (beta * X_B)) / (alpha + beta);
}

// Helper function to calculate the power terms with edge case handling
double Molecule::calculate_power_term(double Rp, double center, int l, int i) const {
    return std::pow(Rp - center, l - i);
}

double Molecule::calculate_term(double alpha, double beta, double X_A, double X_B,
    int lA, int lB, int i, int j) const {
double Rp = calculate_RP(alpha, beta, X_A, X_B);
double pow_XA = calculate_power_term(Rp, X_A, lA, i);
double pow_XB = calculate_power_term(Rp, X_B, lB, j);

//  correct CNDO2 normalization
double denom = std::pow(2.0 * (alpha + beta), i + j);



return (pow_XA * pow_XB) / denom;
}


double Molecule::Overlap_onedim(double alpha, double beta, double X_A, double X_B, int lA, int lB) const
{
    double prefactor = exponential_prefactor(alpha, beta, X_A, X_B);

    double result = 0.0;

    for (int i = 0; i <= lA; ++i) {
        for (int j = 0; j <= lB; ++j) {
            //  DO NOT filter i + j — Cartesian GTOs include all terms
            int binomial_A = calculate_binomial(lA, i);
            int binomial_B = calculate_binomial(lB, j);

            double term = calculate_term(alpha, beta, X_A, X_B, lA, lB, i, j);

            result += binomial_A * binomial_B * term;
        }
    }

    return result * prefactor;
}

// Calculates the primitive overlap between two Gaussian primitives.
double Molecule::calculate_primitive_overlap(double alpha_k, double alpha_l, const arma::vec& R_A, const arma::vec& R_B, const arma::vec lA, const arma::vec lB) const {
    // Calculate overlap in x, y, and z dimensions
    double S_x = Overlap_onedim( alpha_k, alpha_l, R_A(0), R_B(0), lA(0), lB(0));  // x-direction
    double S_y = Overlap_onedim( alpha_k, alpha_l, R_A(1), R_B(1), lA(1), lB(1));  // y-direction
    double S_z = Overlap_onedim( alpha_k, alpha_l, R_A(2), R_B(2), lA(2), lB(2));  // z-direction

    // Return total overlap as the product of the one-dimensional overlaps
    return S_x * S_y * S_z;
}

// Calculates the contracted overlap between two basis functions.
double Molecule::calculate_contracted_overlap(int mu, int nu) const {
    double S_mu_nu = 0.0;
    
    // Loop over primitives in basis function mu
    for (int k = 0; k < 3; ++k) {
        double d_mu_k = basisFunctions_[mu].coeffs[k];
        double N_mu_k = basisFunctions_[mu].norms[k];
        double alpha_mu_k = basisFunctions_[mu].alphas[k];

        // Loop over primitives in basis function nu
        for (int l = 0; l < 3; ++l) {

            double d_nu_l = basisFunctions_[nu].coeffs[l];
            double N_nu_l = basisFunctions_[nu].norms[l];
            double alpha_nu_l = basisFunctions_[nu].alphas[l];

            // Get the centers of the primitives
            arma::vec R_A = basisFunctions_[mu].center;  // Coordinates of center of primitive k
            arma::vec R_B = basisFunctions_[nu].center;  // Coordinates of center of primitive l

            // Get the angular momentum of the primitives
            arma::vec lA = basisFunctions_[mu].getL();
            arma::vec lB = basisFunctions_[nu].getL();

            // Calculate the primitive overlap integral S_kl
            double S_kl = calculate_primitive_overlap(alpha_mu_k, alpha_nu_l, R_A, R_B, lA, lB);


            // Add the contribution to the total overlap
            S_mu_nu += d_mu_k * d_nu_l * N_mu_k * N_nu_l * S_kl;
        }
    }

    return S_mu_nu;  // Return the total overlap integral for the contracted basis functions
}

// Function to calculate overlay matrix
void Molecule::calculate_overlap_matrix() {
    // Loop over all pairs of basis functions (mu, nu)
    for (int mu = 0; mu < N_; ++mu) {
        for (int nu = 0; nu < N_; ++nu) {
            // Calculate the contracted overlap integral between basis function mu and nu
            double S_mu_nu = calculate_contracted_overlap(mu, nu);

            // Store the result in the overlap matrix
            S_overlap_matrix_(mu, nu) = S_mu_nu;
        }
    }
}