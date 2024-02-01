# Colloidal-Patchy-Paticle-Analysis
C++ code written for cluster analysis for the self assembly of colloidal particles - written for computational chemistry PhD

Overview
This program is designed to analyze clusters in a simulation environment, specifically focusing on calculating gyration tensors and related quantities like asphericity, cylindricity, and anisotropy. It takes input data regarding the size and positions of the largest clusters per frame from simulation files and computes various physical properties of these clusters.

Configuration and Usage
Input Files
The program requires two main input files:

lrgst_x_clstr.dat: Contains the sizes of the largest clusters per frame.
pos_x_clstr.dat: Contains the positions of clusters per frame. (x, y, z)
These files should be formatted correctly and placed in the program's working directory or at the specified paths.

Adjustable Variables
The program allows you to adjust several key parameters:

npart: Number of particles.
n_frames: Number of frames in the simulation.
density: Density of the particles.
These can be modified in the program's source code to match the specific details of your simulation.

Output Files
The program generates several output files, each containing different types of data:

TEST.dat: General output for testing and verification.
TRAJ_of_GY: Trajectory of the gyration radius over time.
REFINED_ASPHERICITY: Calculated asphericity values.
REFINED_SPHERICITY2: Alternative asphericity calculations.
REFINED_ACYLINDRICITY: Cylindricity values.
REFINED_AISOTOPY: Anisotropy calculations.
PrincipalMomentX/Y/Z: Principal moments along X, Y, and Z axes.
Main Functionalities
File Reading: The program reads cluster sizes and positions from the input files.

Gyration Tensor Calculation: For each cluster, a gyration tensor is computed using its size and position data.

Physical Properties Calculation: Utilizing the gyration tensor, the program calculates physical properties like the radius of gyration, asphericity, cylindricity, and anisotropy.

Output Generation: Results are written to various output files for further analysis.

Running the Program
To run the program, compile the source code with a C++ compiler that supports C++11 or later and Eigen library. Ensure that the input files are in the correct format and location before executing the program.

Dependencies
The program relies on the following libraries:

Eigen: For matrix and vector operations.
Boost: Specifically, the Graph library for handling graph-related data structures.
Make sure these libraries are properly installed and configured in your development environment.

Notes
The program assumes a specific format for input files. Ensure your data matches the expected format.
The program's output is extensive and meant for users familiar with cluster analysis in particle simulations.
Modify the source code to tailor the analysis to specific requirements or to handle different input formats.
This README provides a basic guide to using and understanding the program. For detailed information about specific functions or calculations, refer to the comments within the source code.





