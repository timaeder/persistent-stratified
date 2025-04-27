#include <iostream>
#include <vector>
#include "newPersCoH.hh"

int main() {
    // Define a small set of points in 2D space
    std::vector<std::vector<long double>> points = {
        {0.0, 0.0},
        {1.0, 1.0},
        {2.0, 2.0},
        {3.0, 3.0}
    };

    // Define epsilon and maximum simplex dimension (k)
    long double epsilon = 5.0; // Increased epsilon
    int k = 2;

    // Define a query point for local copairings
    std::vector<long double> queryPoint = {1.5, 1.5};
    long double radius = 3.0; // Increased radius

    // Create a Complex object
    Complex complex(k + 2);

    // Compute neighborhoods
    Neighborhood neighborhoodCalculator(points, epsilon);
    auto neighborhoods = neighborhoodCalculator.computeNeighborhood();
    complex.setNeighborhood(neighborhoods);

    // Debug: Print neighborhoods
    std::cout << "Neighborhoods:" << std::endl;
    for (const auto& neighbors : neighborhoods.front()) {
        for (const auto& neighbor : neighbors) {
            std::cout << neighbor << " ";
        }
        std::cout << std::endl;
    }

    // Populate IDs map by computing IncrementalVR
    complex.computeIncrementalVR(points, epsilon, k);

    // Debug: Print simplices
    for (size_t d = 0; d < complex.simplices.size(); ++d) {
        std::cout << "Simplices (dimension " << d << "):" << std::endl;
        for (const auto& simplex : complex.simplices[d]) {
            std::cout << "  Diameter: " << simplex.first << ", Simplex: ";
            for (const auto& vertex : simplex.second) {
                std::cout << vertex << " ";
            }
            std::cout << std::endl;
        }
    }

    // Test SimpCoBound
    std::vector<int> simplex = {0, 1};
    auto cobound = SimpCoBound(simplex, neighborhoods.front(), complex.IDs, k);
    std::cout << "SimpCoBound result: ";
    for (const auto& id : cobound) {
        std::cout << id << " ";
    }
    std::cout << std::endl;

    // Test RelSimpCoBound
    auto relCobound = RelSimpCoBound(points, simplex, neighborhoods.front(), complex.IDs, k, epsilon);
    std::cout << "RelSimpCoBound result: ";
    for (const auto& id : relCobound) {
        std::cout << id << " ";
    }
    std::cout << std::endl;

    // Test AddBndry
    std::vector<int> cbID1 = {1, 2, 3};
    std::vector<int> cbID2 = {2, 3, 4};
    AddBndry(cbID1, cbID2);
    std::cout << "AddBndry result: ";
    for (const auto& id : cbID1) {
        std::cout << id << " ";
    }
    std::cout << std::endl;

    // Test Copairings
    auto copairings = Copairings(points, complex, epsilon, k);
    std::cout << "Copairings result: " << std::endl;
    for (size_t d = 0; d < copairings.size(); ++d) {
        std::cout << "Dimension " << d << ":" << std::endl;
        for (const auto& pair : copairings[d]) {
            std::cout << "  (" << pair.first << ", " << pair.second << ")" << std::endl;
        }
    }

    // Test CLocGCopairings
    auto localCopairings = CLocGCopairings(points, k, queryPoint, radius);
    std::cout << "CLocGCopairings result: " << std::endl;
    for (size_t d = 0; d < localCopairings.size(); ++d) {
        std::cout << "Dimension " << d << ":" << std::endl;
        for (const auto& pair : localCopairings[d]) {
            std::cout << "  (" << pair.first << ", " << pair.second << ")" << std::endl;
        }
    }

    // Test BndLocGCopairings
    auto boundaryCopairings = BndLocGCopairings(points, k, queryPoint, radius);
    std::cout << "BndLocGCopairings result: " << std::endl;
    for (size_t d = 0; d < boundaryCopairings.size(); ++d) {
        std::cout << "Dimension " << d << ":" << std::endl;
        for (const auto& pair : boundaryCopairings[d]) {
            std::cout << "  (" << pair.first << ", " << pair.second << ")" << std::endl;
        }
    }

    return 0;
}