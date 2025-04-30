

#include "Atom.h"
#include "BasisFunction.h"
#include <armadillo>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <set>
#include <map>

#ifndef MOLECULE_H
#define MOLECULE_H



// Define a struct to hold the semi-empirical parameters for each element
struct CNDO2Parameters {
    double IA_s;  ///< Semi-empirical parameter for s-orbital (1/2 (Is + As)).
    double IA_p;  ///< Semi-empirical parameter for p-orbital (1/2 (Ip + Ap)).
    double beta;  ///< Semi-empirical parameter for bonding interaction (-β).
};




class Molecule {
public:

    Molecule(const std::vector<Atom>& Atoms);


    friend std::ostream& operator<<(std::ostream& os, const Molecule& molecule);


    void addAtom(const Atom& atom);

    
    std::vector<Atom> getAtoms() const;

  
    void printBasisFunctions() const;

    
    void printS() const;

    
    void printGamma() const;
    

    void printH() const;

  
    double calculate_contracted_overlap(int mu, int nu) const;

   
    double calculate_primitive_overlap(double alpha_k, double alpha_l, const arma::vec& R_A, const arma::vec& R_B, const arma::vec lA, const arma::vec lB) const;

  
    void calculate_overlap_matrix();


    double calc_totalE();

   
    void build_Fmatrix_wSCF(bool verbose = false);

  
    void printEnergies(const std::string& file) const;


    void compute_cndo2_gradient();

    double derivative_primitive_overlap(double alpha, double beta, const arma::vec& RA, const arma::vec& RB, int dim) const;

private:
   
    void buildBasisFunctions();
  
    void Molecule_init();            

    // Class members for storing atoms, basis functions, matrices, and properties
    std::vector<Atom> atoms_;                   ///< Vector containing the atoms of the molecule.
    std::vector<BasisFunction> basisFunctions_; ///< List of basis functions in the molecule.

    arma::mat S_overlap_matrix_;                ///< Overlap matrix for basis functions (S_uv).
    arma::mat FA_matrix_;                       ///< Alpha Fock matrix.
    arma::mat FB_matrix_;                       ///< Beta Fock matrix.
    arma::mat pA_matrix_;                       ///< Alpha density matrix.
    arma::mat pB_matrix_;                       ///< Beta density matrix.
    arma::mat H_core_matrix_;                   ///< Core Hamiltonian matrix.

    arma::mat pTot_mu_nu_;                      ///< Total density matrix for basis functions.
    arma::vec pTot_;                            ///< Total density vector for each atom.
    arma::mat VA_matrix_;                       ///< Alpha eigenvectors matrix from SCF solution.
    arma::mat VB_matrix_;                       ///< Beta eigenvectors matrix from SCF solution.
    arma::vec epsilon_A_;                       ///< Alpha eigenvalues vector (orbital energies).
    arma::vec epsilon_B_;                       ///< Beta eigenvalues vector (orbital energies).
    arma::mat cA_matrix_;                       ///< Alpha molecular orbital (MO) coefficients matrix.
    arma::mat cB_matrix_;                       ///< Beta molecular orbital (MO) coefficients matrix.
    arma::mat gamma_matrix_;                    ///< Gamma matrix for electron repulsion integrals.

    arma::mat x_matrix_;                        ///< X matrix for CNDO/2 gradient computation (electron overlap terms).
    arma::mat y_matrix_;                        ///< Y matrix for CNDO/2 gradient computation (electron repulsion terms).

    arma::mat S_derivative_matrix_;             ///< Derivative of the overlap matrix with respect to atomic coordinates.
    arma::mat gradient_;                        ///< Total energy gradient matrix (nuclear + electronic contributions).
    arma::mat gammaAB_RA_;                      ///< Derivatives of the gamma matrix with respect to atomic coordinates.
    arma::mat S_derivative_matrix_RA_;          ///< Overlap matrix derivatives with respect to atomic coordinates.


    int N_ = 0;                                 ///< Number of basis functions in the molecule.
    int num_atoms_ = atoms_.size();
    int val_electrons_ = 0;                     ///< Number of valence electrons in the molecule.
    int numCarbons_ = 0;                        ///< Number of carbon atoms.
    int numHydrogens_ = 0;                      ///< Number of hydrogen atoms.
    int alpha_e_, beta_e_ = 0;                  ///< Alpha and beta electron counts.
    double alpha_electron_E_, beta_electron_E_ = 0; ///< Energy associated with alpha and beta electrons.
    int charge_ = 0;                            ///< Total charge of the molecule.
    double nuclear_repulsion_E_ = 0;            ///< Nuclear repulsion energy.
    double total_energy_ = 0;                   ///< Total electronic energy.
    std::unordered_map<int, CNDO2Parameters> semiEmpiricalParams_ = { // Map with semi-empirical parameters for H, C, N, O, F
        {1, {7.176, 0.0, -9}},    // H: Ionization energy and beta (only s-orbital)
        {6, {14.051, 5.572, -21}}, // C: Ionization and p-orbital values
        {7, {19.316, 7.275, -25}}, // N: Ionization and p-orbital values
        {8, {25.390, 9.111, -31}}, // O: Ionization and p-orbital values
        {9, {32.272, 11.080, -39}} // F: Ionization and p-orbital values
    };

    // Helper Functions for build and calculate functions
    
    double exponential_prefactor(double alpha, double beta, double X_A, double X_B) const;
   
    int factorial(int n) const;
    
    int double_factorial(int n) const;
   
    double calculate_binomial(int l, int i) const;
    
    double calculate_RP(double alpha, double beta, double X_A, double X_B) const;
  
    double calculate_power_term(double Rp, double center, int l, int i) const;
   
    double calculate_term(double alpha, double beta, double X_A, double X_B, int lA, int lB, int i, int j) const;
   
    double Overlap_onedim(double alpha, double beta, double X_A, double X_B, int lA, int lB) const;
    
    void calculate_gamma();
   
    void build_core_hamiltonian();
   
    int getN() const;
  
    int getValElectrons() const;

    // Hartree-Fock matrix functions (hartreeFock)
   
    double diagonal_value(int mu, bool isAlpha) const;
  
    double offDiagonal_value(int mu, int nu, bool isAlpha) const;
   
    arma::mat build_fock_matrix(bool isAlpha);
  
    void calculate_alpha_beta_electrons();
    
    double I2e_pG(arma::vec &Ra, arma::vec &Rb, double sigmaA, double sigmaB) const;
   
    double eval_gamma(const BasisFunction& basisFuncA, const BasisFunction& basisFuncB) const;
   
    void calculate_pTot();
   
    void printSCF(bool converged, int it) const;


    void calculate_x_matrix();
  
    void calculate_y_matrix();
   
    double I2e_pG_wrtRA(arma::vec &Ra, arma::vec &Rb, double sigmaA, double sigmaB, int coordIndex) const;
    
    double eval_gamma_wrtRA(const BasisFunction& basisFuncA, const BasisFunction& basisFuncB, int coordIndex) const;    
  
    arma::vec compute_vnuc_derivative(int A);
    
    void print_Xuv_yAB() const;
   
    void calculate_gammaAB_RA() ;
   
    void calculate_overlap_derivatives_RA();
    
    double calculate_primitive_overlap(double alpha_k, double alpha_l, const arma::vec& R_A, const arma::vec& R_B) const ;
  
    double derivative_contracted_overlap(int mu, int nu, int direction ) const;
   
    double derivative_primitive_overlap(double alpha, double beta, const arma::vec& R_A, const arma::vec& R_B, 
                                              const arma::vec& lA, const arma::vec& lB, int direction) const;
};

void calculateBondEnergy(std::unordered_map<std::string, double>& energies);

std::string getMoleculeName(const std::string& filePath);

#endif