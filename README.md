# persistent-stratified
This project provides a collection of functions and tools to approximate stratifications and determine **persistent stratified homotopy types** from point cloud data. These methods are based on the theoretical framework introduced by Mäder and Waas in their paper *"From Samples to Persistent Stratified Homotopy Types"*.

# Local Homology Computation Project

This project computes **local homology** using **Vietoris-Rips complexes** for a given set of data points. The program processes a dataset, computes **local cohomology pairings** for each point, and outputs the results to user-specified files. The computations are parameterized by a maximum radius and dimension, allowing for flexible analysis of the data's topological structure.

## Purpose of the Code

The primary goal of this project is to implement mathematical concepts from the aforementioned theoretical work and make them accessible for practical applications. The code is designed to analyze the local topological structure of data, enabling the investigation of stratified spaces that arise in real-world datasets.

This implementation focuses on constructing **Vietoris-Rips complexes** and computing **local homology** to approximate stratifications. The results can be used to study the geometry and topology of data, particularly in the context of stratification learning.

For a complete understanding of the mathematical background, we recommend referring to the paper.

## How to Use the Code

This project includes an example application based on **local homology**, which demonstrates how to use the provided functionality. The example application is implemented in `LocalHomologyApp.cpp`. It processes a dataset of points, computes local cohomology pairings, and outputs the results for further analysis.

### Example Workflow

#### 1. Compile the Program
Use a C++17-compatible compiler to compile the example application:
```bash
g++ -O3 LocalHomologyApp.cpp -o LocalHomologyApp
```

#### 2. Run the Program
Run the compiled program using the following command:
```bash
./LocalHomologyApp
```

The program will prompt you for the following inputs:

##### Input File Path
Enter the path to the data file containing the dataset.
The file should be in CSV format, where each row represents a point in Euclidean space.

**Example:**
```bash
Enter the path to the data file: Data\CrossData.txt
```

##### Maximum Radius
Enter the maximum radius (a positive floating-point number) for constructing the Vietoris-Rips complex.

**Example:**
```bash
Enter the maximum radius (double type, e.g., 0.1): 0.1
```

##### Maximum Dimension
Enter the maximum dimension (a non-negative integer) for simplices in the complex.

**Example:**
```bash
Enter the maximum dimension (integer type): 1
```

##### Output Directory
Enter the directory where the output files should be saved.
If you press Enter without specifying a directory, the program will use the current directory.

**Example:**
```bash
Enter the output directory (press Enter to use the current directory): Output/
```

#### 3. Output
The program generates one output file per point in the dataset.
Each file contains the local cohomology pairings for the corresponding point.
The files are named LocCoH_output<i>.txt, where <i> is the index of the point in the dataset.

**Example:**
```bash
Processing complete. Output files have been saved in: Output/
```

#### 4. Error Handling  
The program includes robust error handling:

- Checks if the input file exists.  
- Validates the radius and dimension inputs.  
- Ensures the output directory exists or defaults to the current directory.  
- Handles parsing errors for the input data file.  

### Input File Format  
The input file should be a CSV file where:  

- Each row represents a point in Euclidean space.  
- Each column corresponds to a coordinate of the point.  

**Example:**
```
0.0,1.0
1.0,0.0
0.5,0.5
```

### Output File Format  
Each output file contains the local cohomology pairings for a specific point. The format is:
```
0.000000,0.100000,0
0.050000,0.150000,1
```

### Visualization  
To visualize the results, a Python script (`LocCoH_plot.py`) is included. This script demonstrates how to process and visualize the output of the example application.  

### Applications  
This project is particularly suited for stratification learning, where the goal is to analyze the local topological structure of data to identify stratified regions or features. It can also serve as a starting point for researchers and practitioners interested in applying topological methods to investigate stratified spaces in real-world data.

### Further Information
For a detailed explanation of the mathematical concepts and methods implemented in this project, please refer to the paper "From Samples to Persistent Stratified Homotopy Types" by Mäder and Waas.