#include <iostream>
#include <mfem.hpp>
#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <geometrycentral/surface/trace_geodesic.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/surface_point.h>

#include "include/Agent.h"
#include "include/Space.h"
#include "include/utils.h"

using namespace  std;

int main(int argc, char *argv[]) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 2.0 * M_PI);


    string mesh_file = "../input/mesh/alveolar_sac_mesh";
    mfem::Mesh mfem_mesh(mesh_file+".msh");
    Space space(mfem_mesh);
    space.loadGCMeshFromFile(mesh_file+".stl");

    vector<double> local_coords{1./3., 1./3., 1./3.}; //barycentric

    Agent agent(&space, 3000, local_coords, 1);

    double speed=10;
    vector<double> velocity{speed, 0};
    agent.setLocalVelocity(velocity);

    double T = 1000;
    double dt = 0.1;
    double persistance_time = 4;

    int n_steps = int(T/dt);
    int n_persistance = int(persistance_time/dt);

    std::vector<geometrycentral::Vector3> trajectory;
    trajectory.reserve(n_steps);

    geometrycentral::Vector3 pos;
    trajectory.push_back(agent.getGlobalPosition());


    for (int step=1;step<n_steps;step++) {
        agent.move(dt);
        trajectory.push_back(agent.getGlobalPosition());
        cout<<"Step:"<<step<<endl;

        if (step%n_persistance==0) {
            double angle = dis(gen);
            agent.rotateVelocityDirection(angle);
        }
    }

    utils::saveAgentTrajectory(trajectory, "../output/alveolar_trajectory.vtk");
    ofstream out("../output/alveolar_manifold.vtk");
    mfem_mesh.PrintVTK(out);
    out.close();

    return 0;
}