//load libaries
#include <Eigen/Dense>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/histogram.hpp>
#include <boost/histogram/ostream.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <optional>
#include <list>
#include <queue>
#include <unordered_set>
//include namespaces
using namespace Eigen;
using namespace std;
using namespace boost;
using std::map;
using std::vector;



// Define the vertex and edge properties
struct VertexProperties {
    int id;
    int connections;
    int RED_BLUE; 
    //V_CC_Map v_cc_map;
   // vector<vector<int>> cc_vec;
};

struct EdgeProperties {
    int type;
    int RED;
    int BLUE;
    int REDnBLUE;
};

//Desribe graph traits 
typedef adjacency_list<vecS, vecS, undirectedS, VertexProperties, EdgeProperties> Graph;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::edge_descriptor Edge;
typedef std::map<Graph::vertex_descriptor, int> V_CC_Map;

//Create cluster struct 
struct Cluster {
    int no;
    int size;
    std::vector<int> vertices;
    int chain;
    int sheet;

    Cluster() : size(0) {}
};

//instantiate global variables
constexpr int n_part = 2000, n_frames = 5000;
constexpr double density = 0.5;
double Lambda = 1.2;
double Delta = (20 * (M_PI / 180));
int patch_no = 6;
double boxl = std::cbrt(n_part/density);
int max_neighbours = 12;

//read files
std::ifstream pos_file("../pos.dat");
std::ifstream ortn_file("../ortn.dat");

//output files
ofstream CLUSTER_file("CLUSTERS.dat");
ofstream CLUSTER_COLOURS("CLUSTER_COLOURS.dat");
ofstream RED("RED_PROPORTION.dat");
ofstream BLUE("BLUE_PROPORTION.dat");
ofstream SIZE("CLUSTER_SIZE.dat");
ofstream REDnBLUE("REDnBLUE.dat");
ofstream CHAIN_DATA("CHAIN_DATA.dat");
ofstream CHAIN_COUNT("CHAIN_COUNT.dat");

//Instantiate matricies
typedef Matrix<double, Dynamic,3,RowMajor>Matrix3DR;
typedef Matrix<double, Dynamic,4,RowMajor>Matrix4DR;
typedef Matrix<double, Dynamic,Dynamic,RowMajor>MatrixDR;
Matrix3DR positions(n_part,3);
Matrix4DR orientations(n_part,4);
MatrixDR Neighbours(n_part,20);
MatrixDR R_hat(n_part,20);
MatrixDR Relative_ortn(n_part, 20);
MatrixDR RIJ_matrix(n_part, n_part);
Matrix3DR PatchSites(6,3);
MatrixDR HETRO_STORE(n_part, 5);
vector<vector<int>> vertex_matrix;

//Create patch position matrix
void DefineSites(Matrix3DR & PatchSites) {
  PatchSites.row(0) << 1.0, 0.0, 0.0; //Blue (axial)
  PatchSites.row(1) << 0.0, 1.0, 0.0; //Red (equatorial)
  PatchSites.row(2) << 0.0, 0.0, 1.0; //Red (equatorial)
  PatchSites.row(3) << -1.0, 0.0, 0.0; //Blue (axial)
  PatchSites.row(4) << 0.0, -1.0, 0.0; //Red (equatorial)
  PatchSites.row(5) << 0.0, 0.0, -1.0; //Red (equatorial)
}

//Nearest image convention
RowVector3d get_sep_vect(const Vector3d &i, const Vector3d &j)
// Take the coordinates of particles i and j and return their nearest image separation vector
{
    RowVector3d r_ij = (j) - (i);
    // Implement nearest-image convention
    r_ij(0) = r_ij(0) - (boxl * (lround(r_ij(0) / boxl)));
    r_ij(1) = r_ij(1) - (boxl * (lround(r_ij(1) / boxl))); 
    r_ij(2) = r_ij(2) - (boxl * (lround(r_ij(2) / boxl)));
    return (r_ij);
}

//Read pos file
void read_in_positions(const int &n_part, Matrix3DR &positions)
// append the positions in the pos file to the positions matrix
{
    for (int i = 0; i < n_part; i++)
    {
        pos_file >> positions(i, 0);
        pos_file >> positions(i, 1);
        pos_file >> positions(i, 2);
    }
}

//Read ortn file
void read_in_orientations(const int &n_part, Matrix4DR &orientations)
// append the positions in the pos file to the positions matrix
{
    for (int i = 0; i < n_part; i++)
    {
        ortn_file >> orientations(i, 0);
        ortn_file >> orientations(i, 1);
        ortn_file >> orientations(i, 2);
        ortn_file >> orientations(i, 3);
    }
}

//initalise raltive ortn
void init_REl_Orientation(const int &n_part, MatrixDR & Relative_ortn, const int&  max_neighbours)
{
 // initialise Relative_ortn array
        for (int i = 0; i < n_part; i++)
        {
            for (int j = 0; j < max_neighbours; j++)
            {
                {
                   Relative_ortn(i, j) = -1;
                }
            }
        }
}


// Calculates the coordination number (degree) of the first vertex in the graph.
// The coordination number is the number of edges connected to a vertex.
// Note: This function currently only calculates and returns the coordination number
// for the first vertex it encounters and then exits, which may not be the intended behavior.
int calculate_cn(const Graph& graph) {
    Graph::vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
        int cn = 0;
        Vertex v = *vi;
        graph_traits<Graph>::adjacency_iterator vj, vj_end;
        for (tie(vj, vj_end) = adjacent_vertices(v, graph); vj != vj_end; ++vj) {
            cn++;
        }
        return cn; // Returns after calculating for the first vertex only
    }  
    return 0; // Only reached if the graph is empty
}


// Calculates the total number of edges in the graph by summing the 'connections' property
// of each vertex and dividing by 2 (since each edge is counted twice, once for each vertex it connects).
double Edge_number(const Graph& graph) {
    double total_connections = 0.0;
    graph_traits<Graph>::vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
        total_connections += (graph[*vi].connections);
    }
    return total_connections / 2;
}



// Calculates the average degree of vertices in the graph.
// The average degree is the total number of edges connected to vertices, divided by the number of vertices.
// Each edge contributes to the degree of two vertices, hence the division by 2 for each connection count.
double average_connections(const Graph& graph) {
    double total_connections = 0.0;
    int num_vertices = 0;
    graph_traits<Graph>::vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
        total_connections += (graph[*vi].connections / 2);
        num_vertices++;
    }
    return total_connections / num_vertices;
}


// Calculates the average degree of vertices in the graph by summing the degrees (out_degree)
// of all vertices and dividing by the total number of vertices
double connection_byedges(const Graph& graph){
    double edge_number = 0.0;
    double ave_by_edges = 0.0;
    graph_traits<Graph>::vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
        Vertex v = *vi;
        int degree = out_degree(v, graph);
        edge_number+= degree;  
    }
    return ave_by_edges = edge_number/4000;

}


// This function counts the number of homogeneous (same-type) bonding connections in a graph.
// It specifically looks for edges where the 'type' property is set to 1, indicating a specific
// type of bond between vertices, and returns the count of such edges.
// Parameters:
// - Graph& g: A reference to the graph being analysed.
// Returns:
// - int: The number of edges in the graph that have a 'type' property value of 1.
int cumulative_type(Graph& g) {
    int HOMO = 0; // Counter for edges with 'type' == 1, indicating homogeneous bonding.
    // Loop over all edges in the graph to check their 'type' property.
    auto edge_pair = edges(g);
    for (auto it = edge_pair.first; it != edge_pair.second; ++it) {
        // Access the properties of the current edge.
        EdgeProperties& edge_prop = g[*it];
        // Increment HOMO if the edge's 'type' is 1 (homogeneous bonding).
        if (edge_prop.type == 1) {
            HOMO += 1;
        }
        // Note: The cumulative_total variable and its related operations have been removed
        // because they were not used beyond calculation. This simplifies the function
        // to focus solely on counting homogeneous connections.
    }
    // Return the count of homogeneous connections found in the graph.
    return HOMO;
}



// Obtain all clusters as per the graph structure.

std::vector<Cluster> findClusters(const Graph& graph) {
    std::vector<Cluster> clusters; // Stores all clusters found in the graph.
    std::unordered_set<int> visited; // Keeps track of visited vertices to avoid revisiting.

    // Iterate over all vertices in the graph to initiate BFS from each unvisited vertex.
    auto vertices = boost::vertices(graph);
    for (auto it = vertices.first; it != vertices.second; ++it) {
        int vertex = *it;
        // Skip this vertex if it has already been visited.
        if (visited.find(vertex) != visited.end()) {
            continue;
        }

        std::queue<int> queue; // Queue used for BFS.
        queue.push(vertex);
        visited.insert(vertex); // Mark the starting vertex as visited.

        Cluster cluster; // Initialise a new cluster.
        while (!queue.empty()) {
            int currentVertex = queue.front();
            queue.pop();
            cluster.vertices.push_back(currentVertex); // Add current vertex to the cluster.

            // Iterate over the neighbors of the current vertex.
            auto neighbors = boost::adjacent_vertices(currentVertex, graph);
            for (auto nIt = neighbors.first; nIt != neighbors.second; ++nIt) {
                int neighbor = *nIt;
                // If this neighbor hasn't been visited, mark it as visited and add to the queue.
                if (visited.find(neighbor) == visited.end()) {
                    visited.insert(neighbor);
                    queue.push(neighbor);
                }
            }
        }

        // After finishing BFS for a cluster, set its size and add it to the list of clusters.
        cluster.size = cluster.vertices.size();
        clusters.push_back(cluster);
    }

    // The variable CLUSTER_COUNT and the second Cluster instance seem unnecessary
    // because each cluster's size and contents are already determined within the loop.
    // If each cluster needs a unique identifier, consider adding IDs during cluster creation.

    return clusters; // Return the list of all found clusters.
}


// Evaluates each cluster in the given vector, calculating metrics like total edges,
// RED, BLUE, REDnBLUE edge counts, and determining chain clusters based on BLUE proportion.
void evaluateClusters(std::vector<Cluster>& clusters, const Graph& graph) {
    int CLUSTER_COUNT = 0; // Counter for assigning a unique ID to each cluster
    int totCHAINS = 0; // Counter for the total number of clusters qualifying as chains



    for (auto& cluster : clusters) {
        if (cluster.size > 4) { // Only evaluate clusters larger than a threshold size
            int totalEdges = 0;
            int totalRED = 0;
            int totalBLUE = 0;
            int totalREDnBLUE = 0; // Initialize to 0 for correct accumulation

            CLUSTER_COUNT++;
            cluster.no = CLUSTER_COUNT; // Assign a unique ID to the cluster

            // Iterate over the vertices in the cluster to count edge types
            for (int vertex : cluster.vertices) {
                auto outEdges = boost::out_edges(vertex, graph);
                for (auto eIt = outEdges.first; eIt != outEdges.second; ++eIt) {
                    Edge edge = *eIt;
                    const EdgeProperties& edgeProps = graph[edge];

                    // Accumulate counts for each edge type
                    totalEdges++;
                    totalRED += edgeProps.RED;
                    totalBLUE += edgeProps.BLUE;
                    totalREDnBLUE += edgeProps.REDnBLUE;
                }
            }

            // Calculate proportions of each edge type
            double redProportion = (totalEdges > 0) ? static_cast<double>(totalRED) / totalEdges : 0;
            double blueProportion = (totalEdges > 0) ? static_cast<double>(totalBLUE) / totalEdges : 0;
            double REDnBLUEProportion = (totalEdges > 0) ? static_cast<double>(totalREDnBLUE) / totalEdges : 0;

 

            // Output the evaluated metrics for the cluster
            CLUSTER_COLOURS << "no: " << cluster.no << std::endl;
            CLUSTER_COLOURS << "Cluster size: " << cluster.size << std::endl;
            CLUSTER_COLOURS << "Total edges: " << totalEdges << std::endl;
            CLUSTER_COLOURS << "RED proportion: " << redProportion << std::endl;
            CLUSTER_COLOURS << "RED tot: " << totalRED << std::endl;
            CLUSTER_COLOURS << "BLUE tot: " << totalBLUE << std::endl;
            CLUSTER_COLOURS << "BLUE proportion: " << blueProportion << std::endl << std::endl;
            SIZE << cluster.size << std::endl;
            BLUE << blueProportion << std::endl;
            RED << redProportion << std::endl;
            REDnBLUE << REDnBLUEProportion << std::endl;
        }
    }
}



// Helper function to recursively count consecutive blue-blue connections from a given vertex.
// This function implements a depth-first search (DFS) to explore blue-blue connections.
//
// Parameters:
// - graph, currentVertex, prevVertex, visited: A set of visited edges to avoid revisiting the same edge.
//
// Returns:
// - The maximum number of consecutive blue-blue connections found in any path from the current vertex.
int countBlueBlueConnections(const Graph& graph, Vertex currentVertex, Vertex prevVertex, std::set<Edge>& visited) {
    // Initialize the count for blue-blue connections starting from this vertex.
    int blueBlueCount = 0;

    // Iterate through all outgoing edges of the current vertex.
    for (auto eIt = boost::out_edges(currentVertex, graph); eIt.first != eIt.second; ++eIt.first) {
        Edge edge = *eIt.first;
        const EdgeProperties& edgeProps = graph[edge];

        // Check for a blue-blue connection.
        if (edgeProps.BLUE == 1) {
            Vertex targetVertex = target(edge, graph);

            // Avoid cycles by skipping the previous vertex.
            if (targetVertex != prevVertex && visited.find(edge) == visited.end()) {
                visited.insert(edge); // Mark the edge as visited.

                // Recursively explore from the target vertex, updating the max count.
                int pathCount = 1 + countBlueBlueConnections(graph, targetVertex, currentVertex, visited);
                blueBlueCount = std::max(blueBlueCount, pathCount);
            }
        }
    }

    return blueBlueCount; // Return the max count of consecutive blue-blue connections.
}




// Evaluates clusters to identify and report the longest sequence of blue-blue connections for each.
//
// Parameters:
// - clusters , graph
void evaluateCHAINS(std::vector<Cluster>& clusters, const Graph& graph) {
    // Iterate through each cluster to evaluate blue-blue connections.
    for (auto& cluster : clusters) {
        int largestBlueLine = 0; // Tracks the longest sequence of blue-blue connections in this cluster.

        // Only consider clusters above a certain size for detailed analysis.
        if (cluster.size > 4) {
            for (int vertex : cluster.vertices) {
                std::set<Edge> visited; // Track visited edges to avoid double-counting.

                // Explore outgoing edges from this vertex for blue-blue connections.
                for (auto eIt = boost::out_edges(vertex, graph); eIt.first != eIt.second; ++eIt.first) {
                    Edge edge = *eIt.first;
                    const EdgeProperties& edgeProps = graph[edge];

                    // Check for a blue-blue connection and initiate a DFS search if found.
                    if (edgeProps.BLUE == 1) {
                        Vertex targetVertex = target(edge, graph);
                        int blueBlueCount = countBlueBlueConnections(graph, targetVertex, vertex, visited);
                        largestBlueLine = std::max(largestBlueLine, blueBlueCount);
                    }
                }
            }

            // Output the findings for this cluster.
            // Assuming CHAIN_COUNT and CHAIN_DATA are std::ofstream objects opened elsewhere.
            CHAIN_COUNT << "Cluster: " << cluster.no << " of size: " << cluster.size 
                        << " has the longest blue-blue connection sequence of: " << largestBlueLine << std::endl;
            CHAIN_DATA << largestBlueLine << std::endl;
        }
    }
}


int main ()
{
    //itterate over all frames of simulation
    for (int frame=0;frame<n_frames;frame++)
    {
     
     read_in_positions(n_part, positions);
     read_in_orientations(n_part, orientations); 

        //fill Matrix with patches 
        DefineSites(PatchSites);

        // Create the graph
        Graph graph(n_part);

        // Add vertex properties
        for (int i = 0; i < n_part; i++) {
            Vertex v = vertex(i, graph);
            graph[v].id = i;
            graph[v].connections = 0;
            graph[v].RED_BLUE = 0;
        }

    //loop over all particles
    for (int i=0; i < n_part; i++)
    {
        Vertex v_i = vertex(i, graph);
        for (int j = 0; j < n_part; j++)
        {
            if (i !=j)
            {       //obtain notmalised distance 
                    double dist = get_sep_vect(positions.row(i), positions.row(j)).norm();
                    Vector3d norm_dist = get_sep_vect(positions.row(i), positions.row(j)).normalized();
                    if (dist <= Lambda)
                    {
                        Vertex v_j = vertex(j, graph);
                        int aligned = 0;
                        //find realtive orientation via quaternoid
                        Quaterniond qi(orientations.row(i).x(), orientations.row(i).y(),
                                    orientations.row(i).z(), orientations.row(i).w());
                        Quaterniond qj(orientations.row(j).x(), orientations.row(j).y(),
                                    orientations.row(j).z(), orientations.row(j).w());
 
                        for (int ni = 0; ni < patch_no; ni++) 
                                {
                                // check particle i patch in to see if within half angle 
                                Vector3d patch = PatchSites.row(ni);
                                Vector3d in_site = qi * patch;
                                //cout << in_site - patch << '\n';
                                bool in_fail = (-(in_site.dot(norm_dist)) <= cos(Delta));
                                if (in_fail) {
                                continue;
                                }

                            for (int nj = 0; nj < patch_no; nj++) 
                                                    {
                                    // check particle j patch jn to see if within half angle 
                                    Vector3d jn_site = qj * PatchSites.row(nj);
                                    bool jn_fail = (norm_dist.dot(jn_site) <= cos(Delta));
                                                if (jn_fail) {
                                                    continue;
                                                    } 
                                                    else {
                                                        //If bool false for both chekcs - particles are bonded
                                                        aligned ++; 
                                                        if (aligned != 0){
                                                        //Draw edge
                                                        Edge e; bool success;
                                                        tie(e, success) = add_edge(v_i, v_j, graph);
                                                            if (success) {

                                                                graph[v_i].connections++;
                                                                graph[v_j].connections++;  

                                                                if ((ni == 0 || ni == 3) && (nj == 0 || nj == 3)){
                                                          
                                                                    graph[e].BLUE = 1;
                                                                    graph[e].type = 1;
                                                                    char I = 'B';
                                                                    char J = 'B';
                                                                } 
                                                                else if ((ni == 1 || ni == 2 || ni == 4 || ni == 5) && (nj == 1 || nj == 2 || nj == 3 || nj == 4 || nj == 5)){
                                                                    graph[e].RED = 1;
                                                                    graph[e].type = 1;
                                                                    char I = 'R';
                                                                    char J = 'R';
                                                            }
                                                                else {
                                                                    graph[e].type = 0;
                                                                    graph[e].REDnBLUE = 1;
                                                                }                                                   
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


            int HOMO = cumulative_type(graph);
            double edge_number = Edge_number(graph);
            float HOMO_bonds = (HOMO/edge_number);


            // Printing graph function
            graph_traits<Graph>::edge_iterator ei, ei_end;
            for (tie(ei, ei_end) = edges(graph); 
                ei != ei_end; ++ei) {
                Vertex src = source(*ei, graph);
                Vertex tgt = target(*ei, graph);
                //Graph_file << graph[src].id << " -- " << graph[tgt].id << " connections: " << graph[src].connections/2 << " Type: " << graph[*ei].type << endl;
            }

            //Connection printing function
            graph_traits<Graph>::vertex_iterator vi, vi_end;
            for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi) {
                Vertex v = *vi;
                //Verticies_file << graph[v].id << " connections: " << graph[v].connections/2 << endl;
            }

            //call function to calculate connectivity of verticies/particles
            calculate_cn(graph);
            double avg_conn = (average_connections(graph));
            connection_byedges(graph);
            

            //store clusters as vector of Cluster struct
            std::vector<Cluster> clusters = findClusters(graph);
            // Displaying clusters function
            int tot = 0;
            for (const auto& cluster : clusters) {
                CLUSTER_file << "Cluster size: " << cluster.size << std::endl;
                CLUSTER_file << "Vertices: ";
                for (int vertex : cluster.vertices) {
                    CLUSTER_file << vertex << " ";
                }
                CLUSTER_file << std::endl;
                tot += cluster.size; 
                CLUSTER_file << tot << endl;
            }
            evaluateClusters(clusters, graph);
            evaluateCHAINS(clusters, graph);
    }
}




