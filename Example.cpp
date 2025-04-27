#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include "Code/newPersCoH.hh"  // Ensure this includes the updated classes/functions
#include "Code/Neighborhood.hh"  // Include this if you need direct access to NeighborhoodGraph
#include "Code/Reader.hh"

// For parallelization, ensure this header is available and correctly set up
#include "Code/par_for.hh"

int main() {
    // Input and output paths
    std::string path_in = "Data/";
    std::string path_out = "Data/";

    // Input data file
    std::string data_filename = "CrossData.txt";

    // Parse the input data using the updated function
    std::vector<std::vector<long double>> data = parse2DCsvFile<long double>(path_in + data_filename);

    // Maximum dimension for homology
    int n = 1;

    // Local radius for local homology computation
    long double rad = 0.1;

    // File prefix for output files
    std::string fileprefix = "LocCoHCross_output";

    // Process each point in the dataset
    for (size_t i = 0; i < data.size(); i++) {
        std::vector<long double> z = data[i];

        // Create an instance of the SimplicialComplex class
        SimplicialComplex simplicialComplex(n);

        // Compute local copairings using updated method
        std::vector<std::vector<std::pair<long double, long double>>> LCP = simplicialComplex.computeLocalCopairings(data, n, z, rad);

        // Output results to a file
        std::ofstream myfile;
        std::string filename = fileprefix + std::to_string(i) + ".txt";
        myfile.open(path_out + filename);
        if (!myfile.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            continue;
        }

        myfile << std::fixed << std::setprecision(6);
        for (size_t j = 0; j < LCP.size(); j++) {
            if (!LCP[j].empty()) {
                std::cout << "Total amount of local " << j << "-pairs: " << LCP[j].size() << std::endl;
                for (const auto& pair : LCP[j]) {
                    myfile << pair.first << "," << pair.second << "," << j << std::endl;
                }
            }
        }
        myfile.close();
    }

    // Uncomment the following for parallel processing
    // pl::async_par_for(0, data.size(), [&](unsigned i) { ... }, false);

    return 0;
}
