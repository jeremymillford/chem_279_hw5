  //Implementation of a Molecule initialization functions. 



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


// Constructor and intializing Functions
Molecule::Molecule(const std::vector<Atom>& Atoms)
    : atoms_(Atoms){
    // Initalize process for MO
    Molecule_init();
}

// Function to generate electron configuration considering only valence electrons
std::vector<std::vector<int>> generateValenceElectronConfiguration(int valenceElectrons) {
    // Electron capacity for valence orbitals (e.g., 2s, 2p for second row elements)
    std::vector<std::vector<int>> orbitalCapacities = {{2}, {2, 2, 2}}; // 2s and 2p
    std::vector<std::vector<int>> electronConfiguration = {{0}, {0, 0, 0}};

    int remainingElectrons = valenceElectrons;

    for (int i = 0; i < orbitalCapacities.size() && remainingElectrons > 0; ++i) {
        for (int j = 0; j < orbitalCapacities[i].size() && remainingElectrons > 0; ++j) {
            int electronsToFill = std::min(remainingElectrons, orbitalCapacities[i][j]);
            electronConfiguration[i][j] = electronsToFill;
            remainingElectrons -= electronsToFill;
        }
    }

    return electronConfiguration;
}

// Function to count unpaired electrons in a given electron configuration
int countUnpairedElectrons(const std::vector<std::vector<int>>& electronConfiguration) {
    int unpairedElectrons = 0;

    // Loop through each orbital in the configuration
    for (const auto& orbital : electronConfiguration) {
        for (int electrons : orbital) {
            // Check if the orbital has unpaired electrons
            if (electrons % 2 != 0) {
                unpairedElectrons++;
            }
        }
    }
    return unpairedElectrons;
}

// Function to determine the multiplicity based on the number of unpaired electrons
int calculateMultiplicity(int unpairedElectrons) {
    return unpairedElectrons == 0 ? 1 : (2 * unpairedElectrons + 1);
}

// Function to determine the number of alpha and beta electrons based on atoms of the molecule
void Molecule::calculate_alpha_beta_electrons() {
    int valenceElectrons = 0;

    for (const auto& atom : atoms_) {
        int atomElectrons = atom.getAtomicNumber();
        
        // Only consider valence electrons (ignoring 1s core electrons)
        if (atomElectrons > 2) {
            valenceElectrons += (atomElectrons - 2); // Subtract 2 for the 1s electrons
        } else {
            valenceElectrons += atomElectrons; // For hydrogen, take all electrons
        }
    }

    std::vector<std::vector<int>> electronConfiguration = generateValenceElectronConfiguration(valenceElectrons);
    int unpairedElectrons = countUnpairedElectrons(electronConfiguration);
    int multiplicity = calculateMultiplicity(unpairedElectrons);

    // Calculate the number of alpha and beta electrons based on total electron count and unpaired electrons
    alpha_e_ = (valenceElectrons + unpairedElectrons) / 2;
    beta_e_ = (valenceElectrons - unpairedElectrons) / 2;

    if (alpha_e_ < 0 || beta_e_ < 0) {
        throw std::runtime_error("Calculated alpha or beta electrons are negative, check your input and calculations.");
    }

    std::cout << "Valence Electrons: " << valenceElectrons << std::endl;
    std::cout << "Unpaired Electrons: " << unpairedElectrons << std::endl;
    std::cout << "Multiplicity: " << multiplicity << std::endl;
    std::cout << "Alpha Electrons (p): " << alpha_e_ << std::endl;
    std::cout << "Beta Electrons (q): " << beta_e_ << std::endl;
}

// Builds the basis functions for a given atom.
void Molecule::buildBasisFunctions() {
    
    int pos = 0; // Variable to record the position in the internal atoms list
    for (auto atom : atoms_){
        int E = atom.getAtomicNumber();
        std::vector<std::vector<double>>  coeffs;
        arma::vec center = { atom.getX(), atom.getY(), atom.getZ() };
        arma::vec L_s = { 0, 0, 0 };    // s: l=m=n=0
        arma::vec L_px = { 1, 0, 0 };   // px: l=1, m=n=0
        arma::vec L_py = { 0, 1, 0 };   // py: l=0, m=1, n=0
        arma::vec L_pz = { 0, 0, 1 };   // pz: l=m=0, n=1

        if (E == 1) {
            coeffs = {
                { 0.15432897, 0.53532814, 0.44463454 } // 1s orbital
            };
        }
        else if (E == 6 || E == 7 || E == 8 || E == 9){
            coeffs = {
                { -0.09996723, 0.39951283, 0.70011547 }, // 2s orbital
                { 0.15591627, 0.60768372, 0.39195739 }   // 2p orbitals
            };
        }
        else {
            throw std::invalid_argument("Error: Basis functions coeffs are not defined for atomic number " + std::to_string(E));
        }

        if (E == 1) {
            // Basis functions for H (1s orbital)
            std::vector<double> alphas = { 3.42525091, 0.62391373, 0.16885540 };
            basisFunctions_.emplace_back(pos, E, center, L_s, alphas, coeffs[0]);
        } else if (E == 6) {
            // Basis functions for C (2s and 2p orbitals)
            std::vector<double> alphas = { 2.94124940, 0.68348310, 0.22228990 };

            basisFunctions_.emplace_back(pos, E, center, L_s, alphas, coeffs[0]);
            basisFunctions_.emplace_back(pos, E, center, L_px, alphas, coeffs[1]);
            basisFunctions_.emplace_back(pos, E, center, L_py, alphas, coeffs[1]);
            basisFunctions_.emplace_back(pos, E, center, L_pz, alphas, coeffs[1]);
        } else if (E == 7) {
            // Basis functions for N (2s and 2p orbitals)
            std::vector<double> alphas = { 3.78045590, 0.87849660, 0.28571440 };

            basisFunctions_.emplace_back(pos, E, center, L_s, alphas, coeffs[0]);
            basisFunctions_.emplace_back(pos, E, center, L_px, alphas, coeffs[1]);
            basisFunctions_.emplace_back(pos, E, center, L_py, alphas, coeffs[1]);
            basisFunctions_.emplace_back(pos, E, center, L_pz, alphas, coeffs[1]);
        } else if (E == 8) {
            // Basis functions for O (2s and 2p orbitals)
            std::vector<double> alphas = { 5.03315130, 1.16959610, 0.38038900 };

            basisFunctions_.emplace_back(pos, E, center, L_s, alphas, coeffs[0]);
            basisFunctions_.emplace_back(pos, E, center, L_px, alphas, coeffs[1]);
            basisFunctions_.emplace_back(pos, E, center, L_py, alphas, coeffs[1]);
            basisFunctions_.emplace_back(pos, E, center, L_pz, alphas, coeffs[1]);
        } else if (E == 9) {
            // Basis functions for F (2s and 2p orbitals)
            std::vector<double> alphas = { 6.46480320, 1.50228120, 0.48858850 };

            basisFunctions_.emplace_back(pos, E, center, L_s, alphas, coeffs[0]);
            basisFunctions_.emplace_back(pos, E, center, L_px, alphas, coeffs[1]);
            basisFunctions_.emplace_back(pos, E, center, L_py, alphas, coeffs[1]);
            basisFunctions_.emplace_back(pos, E, center, L_pz, alphas, coeffs[1]);
        } else {
            // Error handling for unsupported atoms
            throw std::invalid_argument("Error: Basis functions for atomic number " + std::to_string(E) + " are not defined.");
        }
        // Increment the atomic pos
        pos++;
    }
}

// Initializes molecule properties, including valence electrons and basis functions.
void Molecule::Molecule_init(){
    // Accummulation number of basis functions and valence electrons for the MO from each atom
    for (auto& atom : atoms_) {
        N_ += atom.getN();
        val_electrons_ += atom.getValElectrons();
    }

    // Build Basis Functions that 
    buildBasisFunctions();

    // Get the p and q numbers for the valence electrons and charge of MO
    calculate_alpha_beta_electrons();

    // Set S_overlap_matrix_ to be a N x N matrix and initialize it to zeros
    S_overlap_matrix_.set_size(N_, N_);
    S_overlap_matrix_.zeros();

    // Evaluate the S_overlap_matrix
    calculate_overlap_matrix();

    // Initialize alpha, beta, and total density matrices to zero (initial guess)
    pA_matrix_.set_size(N_, N_);
    pB_matrix_.set_size(N_, N_);
    pTot_.set_size(atoms_.size());
    cA_matrix_.set_size(N_, N_);
    cB_matrix_.set_size(N_, N_);
    pA_matrix_.zeros();
    pB_matrix_.zeros();
    pTot_.zeros();
    cA_matrix_.zeros();
    cB_matrix_.zeros();

    // // Initialize alpha and beta fock matrices to zero
    FA_matrix_.set_size(N_, N_);
    FB_matrix_.set_size(N_, N_);
    FA_matrix_.zeros();
    FB_matrix_.zeros();

    // Calculate Gamma and build core hamiltonian
    calculate_gamma();
    build_core_hamiltonian();

    // Initialize Energy Derivative matrices
    x_matrix_.set_size(N_, N_);
    x_matrix_.zeros();
    y_matrix_.set_size(atoms_.size(), atoms_.size());
    y_matrix_.zeros();
    S_derivative_matrix_.set_size(3, N_ * N_);
    S_derivative_matrix_.zeros();
    gammaAB_RA_.set_size(3, atoms_.size() * atoms_.size());
    gammaAB_RA_.zeros();
    gradient_.set_size(3, atoms_.size());
    gradient_.zeros();
}

// Adds a new atom to the molecule.
void Molecule::addAtom(const Atom& atom) {
    atoms_.push_back(atom);

    // Reset the private members
    N_ = 0;
    val_electrons_ = 0;
    numCarbons_ = 0;
    numHydrogens_ = 0;

    // Clear the basisFunctions_ vector to remove the previous basis functions
    basisFunctions_.clear();

    // Reinitialize the molecule
    Molecule_init();
}