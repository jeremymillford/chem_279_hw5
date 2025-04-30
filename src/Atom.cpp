//Atom Class Implementation


#include "Atom.h"
#include <cmath>
#include <iostream>

Atom::Atom(int atomicNumber, double x, double y, double z)
    : atomicNumber_(atomicNumber), x_(x), y_(y), z_(z) {
    if (atomicNumber_ == 1) { // Hydrogen contributes 1 basis function and 1 electron
        N_ = 1;
        ValElectrons_ = 1;
    }
    else if (atomicNumber_ == 6) { // Carbon contributes 4 orbital basis functions and 4 electrons
        N_ = 4;
        ValElectrons_ = 4; 
    } else if (atomicNumber_ == 7) { // Nitrogen contributes 4 orbital basis functions  and 5 electrons
        N_ = 4;
        ValElectrons_ = 5;
    } else if (atomicNumber_ == 8) { // Oxygen contributes 4 orbital basis functions  and 6 electrons
        N_ = 4;
        ValElectrons_ = 6;
    } else if (atomicNumber_ == 9) { // Fluorine contributes 4 orbital basis functions  and 7 electrons
        N_ = 4;
        ValElectrons_ = 7;
    } else {
        throw std::invalid_argument("Unsupported element with atomic number: " + std::to_string(atomicNumber_));
    }
}


std::ostream& operator<<(std::ostream& os, const Atom& atom) {
    os << "[" << atom.atomicNumber_ << "](" << atom.x_ << "," << atom.y_ << "," << atom.z_ << ")";
    return os;
}

double Atom::distanceTo(const Atom& other) const {
    double dx = x_ - other.getX();
    double dy = y_ - other.getY();
    double dz = z_ - other.getZ();
    return std::sqrt( (dx * dx) + (dy * dy) + (dz * dz));
}

double Atom::get_dx(const Atom& other) const {return x_ - other.getX();}
double Atom::get_dy(const Atom& other) const {return y_ - other.getY();}
double Atom::get_dz(const Atom& other) const {return z_ - other.getZ();}

int Atom::getValElectrons() const { return ValElectrons_; }
int Atom::getN() const { return N_; }
int Atom::getAtomicNumber() const { return atomicNumber_; }
double Atom::getX() const { return x_; }
double Atom::getY() const { return y_; }
double Atom::getZ() const { return z_; }
std::vector<double> Atom::getCenter() const { return std::vector<double> {x_, y_, z_};}

void Atom::setX(double newX) { x_ = newX; }
void Atom::setY(double newY) { y_ = newY; }
void Atom::setZ(double newZ) { z_ = newZ; }