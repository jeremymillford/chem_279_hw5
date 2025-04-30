

#include "AtomReader.h"
#include "Reader.h"
#include "Molecule.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>  
#include <filesystem>
#include <armadillo>

int main(int argc, char* argv[]) {
    // Create an output file stream to redirect std::cout
    std::ofstream outFile("../student_output/hw5_2_test_log.txt");  // File to store all cout output
    if (!outFile) {
        std::cerr << "Error: Unable to open output file for logging." << std::endl;
        return 1;
    }

    // Redirect std::cout to the output file
    std::streambuf* coutBuffer = std::cout.rdbuf();  // Save old buffer for restoration later
    std::cout.rdbuf(outFile.rdbuf());  // Redirect cout to outFile

    // Print out to terminal the start of the executable
    std::cout << "TESTING" << std::endl;

    // Set the directory path to the input files and display path that is being iterated in
    std::string inputDirPath = "../atoms";
    std::cout << "Input Directory path is \""<< inputDirPath <<"\"" << std::endl;

    // Create a vector to store the individual input file paths
    std::vector<std::string> filePaths;

    // Iterate over the input directory and get the file paths and print the files
    std::cout << "Input Files:" << std::endl;
    for (const auto& entry : std::filesystem::directory_iterator(inputDirPath)) {
        if (entry.is_regular_file()) {  // Only process regular files
            std::filesystem::path filename = entry.path().filename(); // Convert path to filename to print
            filePaths.push_back(entry.path().string());  // Convert path to string and store it
            // Echo filenames found
            std::cout << filename << std::endl;
        }   
    }

    // Create reader object for expected format.
    AtomReader atomreader;

    // Map to store energies with molecule names as keys
    std::unordered_map<std::string, double> energies;

    // Loop over each file and pass it to the Reader function
    for (const auto& filePath : filePaths) {
        
        // Notify Processing has begun
        std::cout << "_______________________________________________________________" << std::endl;
        std::cout << "Processing file: " << filePath << std::endl;
        
        // Try-catch block in case the reader finds an atom not allowed
        try {
            // Call the AtomReader function for each file
            std::vector<Atom> atoms = atomreader.readFile(filePath);

            // Check if any Gaussians were inputted
            if (atoms.empty()) { 
                std::cout << "No atoms in system" << std::endl;
                continue;
            } 

            Molecule molecule(atoms);
        
            // Print result
            std::cout.precision(5);
            std::cout << "_______________________________________________________________" << std::endl;
            std::cout << "Molecule: \n" << molecule << std::endl;

            // // Print the Gamma, Overlap, and H_core matrices for the MO
            // molecule.printGamma();
            // molecule.printS();
            // molecule.printH();

            // Build Fock matrix with Unrestricted Hartree-Fock CNDO/2 method with verbose as true
            molecule.build_Fmatrix_wSCF(false);

            // Calculate the total energy of the molecule. 
            double totalE = molecule.calc_totalE();

            // Print calculate Energies
            molecule.printEnergies(filePath);

            // Compute gradient
            molecule.compute_cndo2_gradient();

        //     // Print Gradient
        //    std::cout << "Gradient: \n" << gradient << std::endl;


            // // Store energy using molecule name as the key
            // std::string molecule_name = getMoleculeName(filePath);
            // energies[molecule_name] = totalE;


        } catch (const std::exception& e) {
            // Catch any exception and print the error message, but continue processing next file
            std::cerr << "Error reading file " << filePath << ": " << e.what() << std::endl;
            std::cerr << "Skipping file and continuing to the next...\n";
        }
    }
    // // Get Bond energy of N2 and O2
    // calculateBondEnergy(energies);

    // Notify that all files have been processed
    std::cout << "_______________________________________________________________" << std::endl;
    std::cout << "\nAll files processed.\n";

    // Restore std::cout to its original buffer (the terminal)
    std::cout.rdbuf(coutBuffer);

    // Close the output file
    outFile.close();
    std::cout << "Output written to student_output \n" << outFile.rdbuf() << std::endl;
    return 0;
}
