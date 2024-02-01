#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <optional>
#include <list>

// Including necessary headers for matrix operations, graph data structures, file IO, and standard utilities.

using namespace Eigen; // Use the Eigen namespace for matrix and vector operations.
using namespace std; // Standard namespace.
using namespace boost; // Boost library namespace for graph operations.
// Avoid potential namespace conflicts by specifying which map and vector to use.
using std::map;
using std::vector;

// Define constants for the number of particles, frames, and density.
constexpr int npart = 2000, n_frames = 500;
constexpr double density = 0.3;
// Calculate the box length based on the number of particles and density.
double boxl = std::cbrt(npart / density);
// Total rows calculated for a matrix holding all positions over time, assuming each particle has 6 position values per frame.
int tot_rows = npart * n_frames * 6;
int total_lines = 25000; // Total number of lines to process, purpose needs clarification.

// File streams for input and output.
std::ifstream clstr_size("./lrgst_x_clstr.dat"); // File containing sizes of the largest clusters.
std::ifstream clstr_pos("./pos_x_clstr.dat"); // File containing positions of clusters.
// Output files for various measurements and properties.
std::ofstream Test_file("TEST.dat");
std::ofstream GYRATION_OVERTIME("TRAJ_of_GY");
std::ofstream SPHER("REFINED_ASPHERICITY");
std::ofstream SPHER2("AREFINED_SPHERICITY2");
std::ofstream CYLIND("REFINED_ACYLINDRICITY");
std::ofstream ANISOTROPY("REFINED_AISOTOPY");
std::ofstream PRINCMOMX("PrincipalMomentX");
std::ofstream PRINCMOMY("PrincipalMomentY");
std::ofstream PRINCMOMZ("PrincipalMomentZ");

// Type definitions for matrices using Eigen.
typedef Matrix<double, Dynamic,-1,RowMajor>MatrixDR; // General dynamic matrix.
typedef Matrix<double, Dynamic,3,RowMajor>Matrix3DR; // Matrix specifically for 3D positions.
std::vector<Matrix3DR> cluster_positions; // Vector to store matrices of cluster positions.
std::vector<std::pair<int, Matrix3DR>> GYRATION_TENSORS; // Vector to store pairs of cluster size and its gyration tensor.


void checkFileExistence() {
    // Correctly open the file by specifying its name in the ifstream constructor.
    std::ifstream clstr_size("./lrgst_x_clstr.dat"); // Example file path/name.

    // Check if the file was successfully opened.
    if (clstr_size.good()) {
        // File exists, so write to an output file as a test.
        std::ofstream output_file("test.dat");
        if (output_file.is_open()) {
            output_file << "found";
            output_file.close();
        } else {
            std::cout << "Failed to open output file." << std::endl;
        }
    } else {
        std::cout << "Input file does not exist or cannot be opened." << std::endl;
    }
    // Close the input file stream.
    clstr_size.close();
}




Eigen::Vector3d ORDER_PARAM(Eigen::Vector3d &principalMoments) {
    // Sort the principal moments in ascending order for consistency.
    if (principalMoments(0) > principalMoments(1)) std::swap(principalMoments(0), principalMoments(1));
    if (principalMoments(1) > principalMoments(2)) std::swap(principalMoments(1), principalMoments(2));
    if (principalMoments(0) > principalMoments(1)) std::swap(principalMoments(0), principalMoments(1));

    // Return the sorted vector of principal moments.
    return principalMoments;
}






void Radius_of_GYRATIONS(std::vector<std::pair<int, Matrix3DR>>& GYRATION_TENSORS, int& no_of_cluster)
{
    std::vector<double> radii; // Vector to store the radius of gyrations for each cluster

    //itterate all clusters
    for (int clusterIndex = 0; clusterIndex < no_of_cluster; clusterIndex++)
    {
        const std::pair<int, Matrix3DR>& GyrationPair = GYRATION_TENSORS[clusterIndex];
        const Matrix3DR& Gyration = GyrationPair.second;
        int size = GyrationPair.first;

        // if (size > 3 )
        // {

        // Get the real parts of eigenvalues
        Eigen::EigenSolver<Eigen::Matrix3d> eigensolver(Gyration);
        Eigen::Vector3cd eigenvalues = eigensolver.eigenvalues();
        Eigen::Vector3d principalMoments = eigenvalues.real(); 
        ORDER_PARAM(principalMoments);

        //calculate radius of gyration
        double Radius_GY = 0;
        for (int i = 0; i < 3; i++)
        {
            double comp = (principalMoments(i));
            Radius_GY += comp;
        }
        radii.push_back(sqrt(Radius_GY)); // Store the radius of gyrations for the current cluster
        Radius_GY = sqrt(Radius_GY);

        //Print Radius of GY to file         
        GYRATION_OVERTIME << Radius_GY << endl; 
        //print pincripal moments
        cout<< principalMoments(2) << " " << principalMoments(1) << " " << principalMoments(0) << endl;

        //store principal moments
        PRINCMOMX << principalMoments(0) << endl;
        PRINCMOMY << principalMoments(1) << endl;
        PRINCMOMZ << principalMoments(2) << endl;
        
        //calculate parameters 
        double asphericity = (principalMoments(2) - (0.5*((principalMoments(0))+ (principalMoments(1))))); 
        double asphericity2 = (((3.0/2.0)*(principalMoments(2))) - (pow(Radius_GY,2)/2)); 
        double acylindicity = ((principalMoments(1)) - principalMoments(0));
        double anisotropy = (((pow(asphericity,2)) + ((3.0/4.0)*(pow(acylindicity,2))))/pow(Radius_GY,4));

        //Test_file << "Cluster " << (clusterIndex) << " " << "SIZE: " << size << " " <<  principalMoments(2) << " " << principalMoments(1) << " " << principalMoments(0) << " - Radius of Gyration: " << radii[clusterIndex] << " APHERICITY2: " << asphericity2 << std::endl;

        //alternative way to calculate anisotropy
            // double numerator = 0;
            // double deNOM = 0;
            // for (int i = 0; i < 3; i++) 
            // {
            //     numerator += pow(principalMoments(i),2);
            //     deNOM += principalMoments(i);

            // }


        //double anisotropy = (((3.0/2.0)*(numerator/(pow(deNOM,2))) - 0.5));
        //ANISOTROPY_DEBUG << anisotropy << " " << " numerator: " << numerator << " denom: " << deNOM << " MOMENTS: " << principalMoments(2) << " " << principalMoments(1) << " " << principalMoments(0) <<  endl; 

        //store parameters
        SPHER << asphericity << endl;
        SPHER2 << asphericity2 << endl;
        CYLIND << acylindicity << endl;
        ANISOTROPY << anisotropy << endl;
    }
}





// Calculates the gyration tensor for each cluster from their particle positions.
void gyration_tens(std::vector<Matrix3DR>& cluster_positions, int& no_of_cluster, std::vector<std::pair<int, Matrix3DR>>& GYRATION_TENSORS) {
    // Iterate over each cluster to compute its gyration tensor.
    for (int clusterIndex = 0; clusterIndex < no_of_cluster; clusterIndex++) {
        // Initialize a 3x3 matrix to store the gyration tensor, starting with zeros.
        Eigen::Matrix3d gyrationTensor = Eigen::Matrix3d::Zero();
        Eigen::Matrix3d Sxy_TOT = Eigen::Matrix3d::Zero();
        // Retrieve the position matrix for the current cluster.
        const Matrix3DR& cluster_pos = cluster_positions[clusterIndex];
        // The size of the cluster is determined by the number of rows in its position matrix.
        int size = cluster_pos.rows();
        cout << "size: " << size << endl;
        
        // Double loop over all particle pairs to accumulate the gyration tensor components.
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                // Nested loop over the tensor dimensions (3D space).
                for (int a = 0; a < 3; a++) {
                    for (int b = 0; b < 3; b++) {
                        // Calculate the separation components, applying periodic boundary conditions.
                        double Sxx = cluster_pos(i, a) - cluster_pos(j, a);
                        Sxx = Sxx - (boxl * (lround(Sxx / boxl))); // Correct for periodic boundary conditions.
                        double Syy = cluster_pos(i, b) - cluster_pos(j, b);
                        Syy = Syy - (boxl * (lround(Syy / boxl))); // Correct for periodic boundary conditions.
                        
                        // Compute the product of separation components to contribute to the tensor.
                        double Sxy = Sxx * Syy;
                        gyrationTensor(a, b) += Sxy;
                    }
                }
            }
        }

        // Normalize the gyration tensor by the total number of pairwise distances squared.
        gyrationTensor = gyrationTensor / (2.0 * pow(size, 2));
        
        // Store the computed gyration tensor along with the cluster size.
        GYRATION_TENSORS.push_back(std::make_pair(size, gyrationTensor));
    }

    // After calculating gyration tensors for all clusters, analyze them further.
    Radius_of_GYRATIONS(GYRATION_TENSORS, no_of_cluster);
}


int main() {
    int line_number = 0; // Track the number of lines read from clstr_size file
    int no_of_cluster = 0; // Count of clusters with size > 0

    // Attempt to move to the beginning of the clstr_size file
    clstr_size.seekg(0);

    // Ensure the file is open before proceeding
    if (clstr_size.is_open()) {
        string line; // To hold each line read from the file
        // Read sizes of clusters line by line
        while (getline(clstr_size, line)) {
            line_number++; // Increment line counter
            int value = stoi(line); // Convert line to integer representing cluster size

            // Process only if cluster size is positive
            if (value > 0) {
                no_of_cluster++; // Increment cluster count
                // Create a matrix to hold the positions of the current cluster
                Eigen::MatrixXd cluster_pos(value, 3);

                // Populate the position matrix from clstr_pos file
                for (int i = 0; i < value; i++) {
                    for (int j = 0; j < 3; j++) { // Simplified with a loop for readability
                        clstr_pos >> cluster_pos(i, j);
                    }
                }
                // Add the populated matrix to the vector of cluster positions
                cluster_positions.push_back(cluster_pos);
                // Optional: Write the cluster position matrix to Test_file for verification
                Test_file << "Cluster size: " << value << std::endl; // Added for clarity in output
                Test_file << cluster_pos << std::endl;
                Test_file << "END MATRIX" << std::endl;
            }
        }
        // Close the clstr_size file after reading
        clstr_size.close();
        // Output the total number of clusters processed
        cout << "Cluster count: " << no_of_cluster << endl;
        // Calculate and analyze gyration tensors for all clusters
        gyration_tens(cluster_positions, no_of_cluster, GYRATION_TENSORS);
    }
return 0;

}










