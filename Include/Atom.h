//Class declarations for Atom class

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <set>
#include <map>
#include <armadillo>

#ifndef ATOM_H
#define ATOM_H

class Atom {
public:

    Atom(int atomicNumber, double x, double y, double z);
    friend std::ostream& operator<<(std::ostream& os, const Atom& atom);
    double distanceTo(const Atom& other) const;
    double get_dx(const Atom& other) const;
    double get_dy(const Atom& other) const;
    double get_dz(const Atom& other) const;
    int getN() const;
    int getValElectrons() const;   
    int getAtomicNumber() const;
    std::vector<double> getCenter() const;
    double getX() const;
    double getY() const;
    double getZ() const;
    void setX(double newX);
    void setY(double newY);
    void setZ(double newZ);
    arma::vec position; 
    // Atom(int atomic_number_, const arma::vec& position_);
 

private:
    int ValElectrons_; ///< Number of valence electrons in atom.
    int N_;            ///< Number of functions for atom.
    int atomicNumber_; ///< The atomic number of the atom.
    double x_;         ///< The x-coordinate of the atom in space.
    double y_;         ///< The y-coordinate of the atom in space.
    double z_;         ///< The z-coordinate of the atom in space.
};

#endif
