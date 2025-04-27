#include <iostream>
#include <vector>
#include "newRipsComplex.hh"

int main() {
    // Define a set of points in 2D space
    std::vector<std::vector<long double>> points = {
        {0.0, 0.0},
        {1.0, 1.0},
        {2.0, 2.0},
        {3.0, 3.0}
    };

    // Define epsilon and maximum simplex dimension (k)
    long double epsilon = 2.5;
    int k = 2;

    // Create a Complex object
    Complex complex(k + 2);

    // Test IncrementalVR
    std::cout << "Testing IncrementalVR..." << std::endl;
    complex.computeIncrementalVR(points, epsilon, k);
    std::cout << "IncrementalVR computation completed." << std::endl;

    // Test RelIncrementalVR
    std::cout << "Testing RelIncrementalVR..." << std::endl;
    complex.computeRelIncrementalVR(points, epsilon, k);
    std::cout << "RelIncrementalVR computation completed." << std::endl;

    // Test BndIncrementalVR
    std::cout << "Testing BndIncrementalVR..." << std::endl;
    complex.computeBndIncrementalVR(points, epsilon, k);
    std::cout << "BndIncrementalVR computation completed." << std::endl;

    // Output results (for demonstration purposes)
    std::cout << "Simplices:" << std::endl;
    for (size_t i = 0; i < complex.simplices.size(); ++i) {
        std::cout << "Dimension " << i << ":" << std::endl;
        for (const auto& simplex : complex.simplices[i]) {
            std::cout << "  Diameter: " << simplex.first << ", Vertices: ";
            for (const auto& vertex : simplex.second) {
                std::cout << vertex << " ";
            }
            std::cout << std::endl;
        }
    }

    return 0;
}