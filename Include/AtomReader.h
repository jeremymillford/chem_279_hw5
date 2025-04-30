//Class declarations for AtomReader functions

#ifndef ATOMREADER_H
#define ATOMREADER_H

#include "Reader.h"
#include "Atom.h"
#include <string>
#include <vector>


class AtomReader : public Reader<Atom> {
protected:

    Atom parseLine(const std::string& line) const override;
};

#endif 