 // implementation of a Molecule utility functions

#include "Atom.h"
#include "Molecule.h"
#include "BasisFunction.h"

#include <cmath>
#include <stdexcept>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <armadillo>

#include <unordered_map>
#include <string>

#include <iostream>
#include <vector>
#include <algorithm> // For std::count

// Get Functions
int Molecule::getN() const {return N_;}
int Molecule::getValElectrons() const {return val_electrons_;}
std::vector<Atom> Molecule::getAtoms() const {return atoms_;}

// Print Functions
std::ostream& operator<<(std::ostream& os, const Molecule& molecule) {
    os << "Basis Functions (N): " << molecule.getN() 
    <<  "\nValence Elections (2n): " << molecule.getValElectrons() 
    <<  "\nAtoms:\n";
    for (auto& atom : molecule.getAtoms())
    os << atom << "\n";
    return os;
}

void Molecule::printBasisFunctions() const {
    int iter = 0;
    for (auto& func : basisFunctions_){
        std::cout << "BasisFunctions " << iter << ":\n" << func << std::endl;
        iter++;
    }
}

void Molecule::printS() const {
    std::cout << "Overlap Matrix:" << std::endl;
    std::cout << S_overlap_matrix_ << std::endl;
}

void Molecule::printH() const {
    std::cout << "H_core Matrix:" << std::endl;
    std::cout << H_core_matrix_ << std::endl;
}

void Molecule::printGamma() const {
    std::cout << "Gamma Matrix:" << std::endl;
    std::cout << gamma_matrix_ << std::endl;
}

void Molecule::printEnergies(const std::string& file) const {
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Molecule specificied in " << file << std::endl;
    std::cout<< "Nuclear repulsion energy of " << nuclear_repulsion_E_ << " eV" << std::endl;
    std::cout << "Electron energy of " << alpha_electron_E_ + beta_electron_E_ << " eV" << std::endl;
    std::cout << "Total energy " << total_energy_ << " eV" << std::endl;
}

void Molecule::printSCF(bool converged, int it) const {
    std::cout << "Iteration: " << it << std::endl;
    std::cout << "p: " << alpha_e_ << ", q: " << beta_e_ << std::endl;

    std::cout << "F_alpha Matrix:\n" << FA_matrix_ << std::endl;
    std::cout << "F_beta Matrix:\n" << FB_matrix_ << std::endl;

    std::cout << "C_alpha Matrix:\n" << cA_matrix_ << std::endl;
    std::cout << "C_beta Matrix:\n" << cB_matrix_ << std::endl;

    std::cout << "P_alpha Matrix:\n" << pA_matrix_ << std::endl;
    std::cout << "P_beta Matrix:\n" << pB_matrix_ << std::endl;

    std::cout << "P_total:\n" << pTot_ << std::endl;

    if (converged) {
        std::cout << "Eigenvalues (Alpha):\n" << epsilon_A_ << std::endl;
        std::cout << "Eigenvalues (Beta):\n" << epsilon_B_ << std::endl;
    }
}
void Molecule::print_Xuv_yAB() const {
    // Print x matrix
    std::cout << "x_matrix_:\n" << x_matrix_ << std::endl;

    // Print x matrix
    std::cout << "y_matrix_:\n" << y_matrix_ << std::endl;

    std::cout << "----------------------------------------------------------" << std::endl;
}