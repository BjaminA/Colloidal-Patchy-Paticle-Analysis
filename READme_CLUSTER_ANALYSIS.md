Overview
This Cluster Analysis Program is a comprehensive tool designed for analysing particle clusters in a simulation environment concerning octahedrally functionalised patches with 'Red' bonds corresponding to equatorial bonds and 'Blue' corresponding to axial. It reads simulation data, builds a graph representation of particle interactions, and calculates various metrics such as cluster size, bonding types, and chain properties. The program is especially focused on analyzing clusters based on their geometric and bonding characteristics.

Features
Graph-Based Representation: Constructs a graph where each particle is a vertex, and edges represent particle interactions.
Cluster Detection: Identifies clusters using a breadth-first search algorithm.
Chain Analysis: Evaluate each cluster to identify and report the longest sequence of specific bonding types (e.g., blue-blue connections).
Metrics Calculation: Computes various metrics like total edges, proportion of different bond types, average connections, and asphericity.
Configuration and Usage
Adjustable Parameters
The program allows the customisation of several key parameters:

n_part: Number of particles in the simulation.
n_frames: Number of frames to be analyzed.
density: Particle density in the simulation.
Lambda, Delta, patch_no, boxl, max_neighbours: Parameters related to the physical properties and interaction criteria of particles.
These parameters can be adjusted in the source code to match the simulation's specifics.

Input Files
The program requires two main input files:

pos.dat: Contains the positions of particles (x, y, z).
ortn.dat: Contains the orientations of particles (quaternions).
Output Files
The program outputs several files, including:

CLUSTERS.dat: Information about the identified clusters.
CLUSTER_COLOURS.dat, RED_PROPORTION.dat, BLUE_PROPORTION.dat, REDnBLUE.dat: Metrics related to the types of bonds in the clusters corresponding to equatorial or axial patches (red and blue respectively).
CHAIN_DATA.dat, CHAIN_COUNT.dat: Information about chain-like structures in the clusters (returning quantities relating to the longest axial-axial (blue-blue) path in a cluster).
Running the Program
Compile the source code with a C++ compiler supporting C++11 or later. Ensure input files are correctly formatted and placed in the appropriate directory. Execute the compiled program to perform the analysis.

Dependencies
Eigen: For matrix and vector operations.
Boost: Specifically, the Graph library for managing graph-related data structures.
Notes
It's crucial to ensure that the input files are correctly formatted and in the right location.
The program is designed for users familiar with particle simulations and cluster analysis.
The source code can be modified for specific requirements or different input formats.
