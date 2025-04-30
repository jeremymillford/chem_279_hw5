//Reader Class implmentation
#include "Reader.h"
#include "Atom.h"
#include "AtomReader.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <stdexcept>

// Derived class that parses the standard format
Atom AtomReader::parseLine(const std::string& line) const {
    std::stringstream ss(line);
    int E;
    double x, y, z;

    // Try to extract exactly 4 values
    if (!(ss >> E >> x >> y >> z)) {
        throw std::invalid_argument("Error: Incorrect format, could not extract 4 values.");
    }

    // Ensure no extra values exist on the line
    std::string extra;
    if (ss >> extra) {
        throw std::invalid_argument("Error: More than 4 values provided on a single line.");
    }

    // Return a Atom object based on the parsed values
    return Atom(E, x, y, z);
}