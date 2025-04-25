#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include "Code/PersCoH.hh"
// For parallelization
#include "Code/par_for.hh"
#include "Code/Reader.hh"

int main()
{
    // Input and output paths
    string path_in = "Data/";
    string path_out = "Data/";

    // Input data file
    string data_filename = "CrossData.txt";

    // Parse the input data
    vector<vector<long double>> data = parse2DCsvFile(path_in + data_filename);

    // Maximum dimension for homology
    int n = 1;

    // Local radius for local homology computation
    long double rad = 0.1;

    // File prefix for output files
    string fileprefix = "LocCoHCross_output";

    // Process each point in the dataset
    // pl::async_par_for(0, data.size(), [&](unsigned i)
    // {
    for (int i = 0; i < data.size(); i++)
    {
        vector<long double> z = data[i];

        // Compute local copairings
        vector<vector<pair<long double, long double>>> LCP = CLocGCopairings(data, n, z, rad);

        // Output results to a file
        ofstream myfile;
        string filename = fileprefix + to_string(i) + ".txt";
        myfile.open(path_out + filename);
        if (!myfile.is_open())
        {
            cerr << "Error: Could not open file " << filename << endl;
            continue;
        }

        myfile << fixed << setprecision(6);
        for (int j = 0; j < LCP.size(); j++)
        {
            if (!LCP[j].empty())
            {
                cout << "Total amount of local " << j << "-pairs: " << LCP[j].size() << endl;
                for (const auto& pair : LCP[j])
                {
                    myfile << pair.first << "," << pair.second << "," << j << endl;
                }
            }
        }
        myfile.close();
    }

    // Uncomment the following for parallel processing
    // pl::async_par_for(0, data.size(), [&](unsigned i) { ... }, false);

    return 0;
}
