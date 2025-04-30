//Class declarations for BasisFunction class

#include "Atom.h"
#include <armadillo>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <set>
#include <map>

#ifndef BASISFUNC_H
#define BASISFUNC_H


class BasisFunction {
public:
    int element;                               
    int val_electrons;                          
    int orbitalType;                            
    int atomic_pos;                             
    arma::vec center;                         
    arma::vec ang_momenta;                     
    std::vector<double> alphas;                
    std::vector<double> coeffs;                 
    std::vector<double> norms;                  


    BasisFunction(const int pos, const int E, const arma::vec& R, arma::vec L, 
                  const std::vector<double>& a, const std::vector<double>& d)
        : atomic_pos(pos), element(E), center(R), ang_momenta(L), alphas(a), coeffs(d) {
        // Compute normalization constants
        for (size_t i = 0; i < alphas.size(); ++i) {
            double norm = computeNormalization(alphas[i], ang_momenta);
            norms.push_back(norm);
        }

        // Set orbital type based on accumulated angular momentum quantum numbers
        orbitalType = arma::accu(ang_momenta);

        // Determine valence electrons based on atomic number
        if (E == 6) { // Carbon
            val_electrons = 4;
        } else if (E == 1) { // Hydrogen
            val_electrons = 1;
        } else if (E == 7) { // Nitrogen
            val_electrons = 5;
        } else if (E == 8) { // Oxygen
            val_electrons = 6;
        } else if (E == 9) { // Fluorine
            val_electrons = 7;
        } else {
            throw std::invalid_argument("Unsupported element with atomic number: " + std::to_string(element));
        }
    }


    arma::vec getL() const { return ang_momenta; }


    friend std::ostream& operator<<(std::ostream& os, const BasisFunction& bf) {
        os << "Basis Function Center: (" << bf.center(0) << ", " << bf.center(1) << ", " << bf.center(2) << ")\n";
        os << "Angular Momentum (l, m, n): (" << bf.ang_momenta(0) << ", " << bf.ang_momenta(1) << ", " << bf.ang_momenta(2) << ")\n";
        os << "Exponents (alphas): ";
        for (const auto& alpha : bf.alphas) os << alpha << " ";
        os << "\n";
        os << "Contraction Coefficients (coeffs): ";
        for (const auto& coeff : bf.coeffs) os << coeff << " ";
        os << "\n";
        os << "Normalization Constants (norms): ";
        for (const auto& norm : bf.norms) os << norm << " ";
        os << "\n";
        return os;
    }
private:

    double computeNormalization(double alpha, arma::vec ang_momenta) {
        return pow(2 * alpha / M_PI, 3.0 / 4.0) *
               pow(4 * alpha, arma::accu(ang_momenta) / 2.0) / 
               (factorial(2*ang_momenta(0) - 1) * factorial(2*ang_momenta(1) - 1) * factorial(2*ang_momenta(2) - 1));
    }


    int factorial(int n) {
        return (n <= 1) ? 1 : n * factorial(n - 2);
    }
};

#endif