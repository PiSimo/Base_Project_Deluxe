#include <iostream>
#include <mfem.hpp>
#include <queue>
#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/surface_point.h>

#include "include/Agent.h"
#include "include/Space.h"
#include "include/utils.h"

using namespace  std;
using namespace geometrycentral::surface;


// Hash function for Vertex
struct VertexHash {
    std::size_t operator()(const Vertex& v) const {
        return std::hash<size_t>()(v.getIndex());
    }
};

std::vector<Vertex> findVerticesWithinRadius(SurfacePoint sp, double radius, Space *space) {
    std::vector<Vertex> withinRadius;

    // Priority queue for Dijkstra's algorithm
    std::priority_queue<std::pair<double, Vertex>, std::vector<std::pair<double, Vertex>>, std::greater<>> pq;

    // Distance map with custom hash function
    std::unordered_map<Vertex, double, VertexHash> distance;

    // Initialize with the closest vertex to SurfacePoint
    Vertex startVertex = sp.nearestVertex();
    pq.push({0.0, startVertex});
    distance[startVertex] = 0.0;

    while (!pq.empty()) {
        auto [dist, v] = pq.top();
        pq.pop();

        // If the distance exceeds the radius, stop processing
        if (dist > radius) continue;

        withinRadius.push_back(v);

        // Traverse neighbors
        for (Halfedge he : v.outgoingHalfedges()) {
            Vertex neighbor = he.tipVertex();
            double edgeLength = space->gc_geometry->edgeLengths[he.edge()];
            double newDist = dist + edgeLength;

            if (distance.find(neighbor) == distance.end() || newDist < distance[neighbor]) {
                distance[neighbor] = newDist;
                pq.push({newDist, neighbor});
            }
        }
    }

    return withinRadius;
}

void writeVerticesToVTK(const std::string& filename, const std::vector<Vertex>& vertices, Space* space) {
    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }

    // VTK Header
    outFile << "# vtk DataFile Version 3.0\n";
    outFile << "Vertex visualization\n";
    outFile << "ASCII\n";
    outFile << "DATASET POLYDATA\n";

    // Write vertex positions
    outFile << "POINTS " << vertices.size() << " float\n";
    for (const Vertex& v : vertices) {
        geometrycentral::Vector3 pos = space->gc_geometry->inputVertexPositions[v];
        outFile << pos.x << " " << pos.y << " " << pos.z << "\n";
    }

    // Write point connectivity (as vertices, without edges or faces)
    outFile << "VERTICES " << vertices.size() << " " << 2 * vertices.size() << "\n";
    for (size_t i = 0; i < vertices.size(); i++) {
        outFile << "1 " << i << "\n";
    }

    outFile.close();
    std::cout << "VTK file saved to " << filename << "\n";
}

int main(int argc, char *argv[]) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 2.0 * M_PI);

    string mesh_file = "../input/mesh/alveolar_sac_mesh";
    mfem::Mesh mfem_mesh(mesh_file+".msh");
    Space space(mfem_mesh);
    space.loadGCMeshFromFile(mesh_file+".stl");

    vector<double> local_coords{1./3., 1./3., 1./3.}; //barycentric

    double radius=0.5;
    Agent agent_1(&space, 100, local_coords, radius);
    Agent agent_2(&space, 200, local_coords, radius);

    double speed=2;

    vector<double> velocity_1{0, speed};
    vector<double> velocity_2{0, -speed};

    agent_1.setLocalVelocity(velocity_1);
    agent_2.setLocalVelocity(velocity_2);

    vector<geometrycentral::Vector3> trajectory_1;
    vector<geometrycentral::Vector3> trajectory_2;

    double T=10;
    double dt=0.1;
    for (int i = 0;i <int(T/dt); i++) {
        agent_1.move(dt);
        agent_2.move(dt);
        if (i%10==0) {
            agent_1.rotateVelocityDirection(dis(gen));
            agent_2.rotateVelocityDirection(dis(gen));
        }
        utils::saveAgentTrajectoryWithRadius({agent_1.getGlobalPosition()},"../output/alveolar_trajectory_agent_1_step"+to_string(i+1)+".vtk",radius);
        vector<Vertex> vertexes = findVerticesWithinRadius(agent_1.gc_position, radius+0.1, &space);
        writeVerticesToVTK("../output/alveolar_neighbourhood_agent_1_step_"+to_string(i+1)+".vtk", vertexes, &space);
    }
    trajectory_1.push_back(agent_1.getGlobalPosition());
    //trajectory_2.push_back(agent_2.getGlobalPosition());

    //utils::saveAgentTrajectoryWithRadius(trajectory_1, "../output/plane_trajectory_agent_1.vtk", radius);
    //utils::saveAgentTrajectoryWithRadius(trajectory_2, "../output/plane_trajectory_agent_2.vtk", radius);

    vector<Vertex> vertexes = findVerticesWithinRadius(agent_1.gc_position, radius, &space);
    //writeVerticesToVTK("../output/neighbourhood_agent_1.vtk", vertexes, &space);
    ofstream out("../output/alveolar_manifold.vtk");
    mfem_mesh.PrintVTK(out);
    out.close();

    return 0;
}