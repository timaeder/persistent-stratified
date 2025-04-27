#include "newNeighborhood.hh"

int main() {
    std::vector<std::vector<long double>> points = {{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}};
    long double epsilon = 1.5;

    Neighborhood neighborhood(points, epsilon);

    auto standardNeighborhood = neighborhood.computeNeighborhood();
    auto cvrNeighborhood = neighborhood.computeCVRNeighborhood();
    auto bndNeighborhood = neighborhood.computeBndNeighborhood();

    // Use the computed neighborhoods as needed
}