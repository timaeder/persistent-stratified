#include <string>
#include <vector>
#include <sstream> // istringstream
#include <iostream> // cout
#include <fstream> // ifstream
#include <stdexcept> // invalid_argument

/**
 * Reads a CSV file into a table, exported as a vector of vector of a specified numeric type.
 * @tparam T Numeric type (e.g., float, double, long double).
 * @param inputFileName Input file name (full path).
 * @return Data as a vector of vector of the specified numeric type.
 */
template <typename T>
std::vector<std::vector<T>> parse2DCsvFile(const std::string& inputFileName) {
    std::vector<std::vector<T>> data;
    std::ifstream inputFile(inputFileName);
    if (!inputFile.is_open()) {
        throw std::invalid_argument("Could not open file: " + inputFileName);
    }

    std::string line;
    int lineNumber = 0;

    while (std::getline(inputFile, line)) {
        lineNumber++;
        if (line.empty() || line[0] == '#') {
            continue; // Skip empty lines or comments
        }

        std::istringstream lineStream(line);
        std::vector<T> record;
        std::string cell;

        while (std::getline(lineStream, cell, ',')) {
            try {
                record.push_back(static_cast<T>(std::stold(cell)));
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid number found in file " << inputFileName
                          << " at line " << lineNumber << ": " << cell << std::endl;
                throw;
            } catch (const std::out_of_range& e) {
                std::cerr << "Number out of range in file " << inputFileName
                          << " at line " << lineNumber << ": " << cell << std::endl;
                throw;
            }
        }

        data.push_back(record);
    }

    if (inputFile.bad()) {
        throw std::ios_base::failure("Error reading file: " + inputFileName);
    }

    return data;
}