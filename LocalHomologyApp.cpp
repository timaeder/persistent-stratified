#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <limits>
#include <sstream>
#include <filesystem>
#include "Code/PersCoH.hh"
#include "Code/par_for.hh"
#include "Code/Reader.hh"

int main()
{
    // User input variables
    std::string data_filepath;
    std::string output_dir;
    long double rad;
    int n;

    // Prompt user for the input file path
    std::cout << "Enter the path to the data file: ";
    std::getline(std::cin, data_filepath);

    // Check if the file exists
    if (!std::filesystem::exists(data_filepath))
    {
        std::cerr << "Error: File not found at " << data_filepath << std::endl;
        return 1;
    }

    // Prompt user for the maximum radius
    std::cout << "Enter the maximum radius (double type, e.g., 0.1): ";
    std::string rad_input;
    std::getline(std::cin, rad_input);
    std::istringstream rad_stream(rad_input);
    if (!(rad_stream >> rad) || rad <= 0)
    {
        std::cerr << "Error: Invalid radius. Please enter a positive number." << std::endl;
        return 1;
    }

    // Prompt user for the maximum dimension
    std::cout << "Enter the maximum dimension (integer type): ";
    std::string dim_input;
    std::getline(std::cin, dim_input);
    std::istringstream dim_stream(dim_input);
    if (!(dim_stream >> n) || n < 0)
    {
        std::cerr << "Error: Invalid dimension. Please enter a non-negative integer." << std::endl;
        return 1;
    }

    // Prompt user for the output directory
    std::cout << "Enter the output directory (press Enter to use the current directory): ";
    std::getline(std::cin, output_dir);

    // Use the current directory if no output directory is specified
    if (output_dir.empty())
    {
        output_dir = std::filesystem::current_path().string();
    }

    // Check if the output directory exists
    if (!std::filesystem::exists(output_dir))
    {
        std::cerr << "Error: Output directory does not exist: " << output_dir << std::endl;
        return 1;
    }

    // Parse the input data file
    std::vector<std::vector<long double>> data;
    try
    {
        data = parse2DCsvFile<long double>(data_filepath);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: Failed to parse the data file. " << e.what() << std::endl;
        return 1;
    }

    // Check if the data is empty
    if (data.empty())
    {
        std::cerr << "Error: The data file is empty or invalid." << std::endl;
        return 1;
    }

    // File prefix for output files
    std::string fileprefix = "LocCoH_output";

    // Process each point in the dataset
    for (size_t i = 0; i < data.size(); i++)
    {
        std::vector<long double> z = data[i];

        // Compute local copairings
        std::vector<std::vector<std::pair<long double, long double>>> LCP = CLocGCopairings(data, n, z, rad);

        // Output results to a file
        std::ofstream myfile;
        std::string filename = output_dir + "/" + fileprefix + std::to_string(i) + ".txt";
        myfile.open(filename);
        if (!myfile.is_open())
        {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            continue;
        }

        myfile << std::fixed << std::setprecision(6);
        for (size_t j = 0; j < LCP.size(); j++)
        {
            if (!LCP[j].empty())
            {
                std::cout << "Total amount of local " << j << "-pairs: " << LCP[j].size() << std::endl;
                for (const auto& pair : LCP[j])
                {
                    myfile << pair.first << "," << pair.second << "," << j << std::endl;
                }
            }
        }
        myfile.close();
    }

    std::cout << "Processing complete. Output files have been saved in: " << output_dir << std::endl;

    return 0;
}