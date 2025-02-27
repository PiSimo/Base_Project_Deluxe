#include <iostream>
#include <mfem.hpp>
#include <queue>
#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/surface_point.h>
#include <geometrycentral/surface/flip_geodesics.h>

#include "include/Agent.h"
#include "include/Space.h"
#include "include/utils.h"

using namespace  std;
using namespace geometrycentral::surface;


int main(int argc, char *argv[]) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 2.0 * M_PI);

    string mesh_file = "plane_mesh";
    mfem::Mesh mfem_mesh("../input/mesh/"+mesh_file+".msh");
    Space space(mfem_mesh);
    space.loadGCMeshFromFile("../input/mesh/"+mesh_file+".stl");

    vector<double> local_coords{1./3., 1./3., 1./3.}; //barycentric

    double radius=0.4;
    Agent agent_01(&space, 100, local_coords, radius, 1);
    Agent agent_02(&space, 200, local_coords, radius, 2);

    double speed=2;

    vector<double> velocity_1{0, speed};
    vector<double> velocity_2{0, -speed};

    agent_01.setLocalVelocity(velocity_1);
    agent_02.setLocalVelocity(velocity_2);

    vector<geometrycentral::Vector3> trajectory_1;
    vector<geometrycentral::Vector3> trajectory_2;

    vector<Agent> agents{agent_01, agent_02};

    vector<int> occupation_matrix(space.gc_mesh->nFaces(), 0);

    double T=40;
    double dt=0.1;
    for (int i = 0;i <int(T/dt); i++) {
        std::fill(occupation_matrix.begin(), occupation_matrix.end(), 0);
        cout<<" step:"<<i<<endl;
        for (int n_agent =0; n_agent <agents.size(); n_agent++) {
            agents[n_agent].move(dt);
            if (i%10==0) {
                agents[n_agent].rotateVelocityDirection(dis(gen));
            }


            unordered_set<Face> faces = agents[n_agent].findFacesWithinRadius(agents[n_agent].gc_position, agents[n_agent].radius, &space);

            for (Face f : faces) {
                int face_index = f.getIndex();
                if (occupation_matrix[face_index] == 0) {
                    occupation_matrix[face_index] = agents[n_agent].agent_index;
                }
                else {
                    Vertex start = agents[n_agent].gc_position.nearestVertex();
                    Vertex end = agents[(n_agent+1)%2].gc_position.nearestVertex();

                    std::unique_ptr<FlipEdgeNetwork> edgeNetwork = FlipEdgeNetwork::constructFromDijkstraPath(*space.gc_mesh, *space.gc_geometry, start, end);
                    edgeNetwork->iterativeShorten();

                    edgeNetwork->posGeom = space.gc_geometry.get();
                    std::vector<std::vector<geometrycentral::Vector3>> polyline = edgeNetwork->getPathPolyline3D();

                    // local-velocities
                    double totalLength = 0.0;

                    vector<geometrycentral::Vector3> percorso;
                    // There may be multiple polylines, so iterate over each
                    for (const auto& path : polyline) {
                        // For each consecutive pair of points, add the distance
                        for (size_t i = 0; i + 1 < path.size(); i++) {
                            geometrycentral::Vector3 p0 = path[i];
                            geometrycentral::Vector3 p1 = path[i + 1];
                            double segmentLength = (p1 - p0).norm();
                            totalLength += segmentLength;
                            percorso.push_back(p0);

                            percorso.push_back(p1);
                        }
                    }

                    utils::saveAgentTrajectory(percorso, "../output/shortest_"+to_string(i)+"_path.vtk");


                    // 3) Compare path length to sum of radii
                    double lengthDifference = (agents[n_agent].radius + agents[(n_agent+1)%2].radius)-totalLength;
                    cout<<"   toatal_length:"<<totalLength<<endl;
                    cout<<"   sum of radii:"<<agents[n_agent].radius + agents[(n_agent+1)%2].radius<<endl;
                    if (lengthDifference <0 )break;
                    const auto& path = polyline[0];

                    // Velocity for agents[n_agent] (start of the path)
                    geometrycentral::Vector3 velocityAgent1 =
                        (percorso[1] - percorso[0]).normalize();

                    // Velocity for agents[(n_agent+1)%2] (end of the path)
                    geometrycentral::Vector3 velocityAgent2 =
                        (percorso[percorso.size()-1] - percorso[percorso.size()-2]).normalize();

                    geometrycentral::Vector2 vel_agents_first,vel_agents_second;
                    space.convertGlobalToLocalVector(agents[n_agent].gc_position, velocityAgent1, vel_agents_first);
                    vel_agents_first = -vel_agents_first.normalize();
                    geometrycentral::Vector2 old_vel_first = agents[n_agent].gc_local_velocity_direction;
                    agents[n_agent].gc_local_velocity_direction = vel_agents_first.normalize();

                    space.convertGlobalToLocalVector(agents[(n_agent+1)%2].gc_position, velocityAgent2, vel_agents_second);
                    vel_agents_second    = -vel_agents_second.normalize();
                    geometrycentral::Vector2 old_vel_second =  agents[(n_agent+1)%2].gc_local_velocity_direction;
                    agents[(n_agent+1)%2].gc_local_velocity_direction = vel_agents_second.normalize();

                    agents[n_agent].move(lengthDifference/(2*speed));
                    agents[(n_agent+1)%2].move(lengthDifference/(2*speed));

                    agents[n_agent].gc_local_velocity_direction = old_vel_first;
                    agents[(n_agent+1)%2].gc_local_velocity_direction = old_vel_second;

                    cout<<"   collision!"<<" overlap-size"<<lengthDifference<<endl;

                    break;
                }
            }

            faces = agents[n_agent].findFacesWithinRadius(agents[n_agent].gc_position, agents[n_agent].radius, &space);

            // visualisation:
            utils::saveAgentTrajectoryWithRadius({agents[n_agent].getGlobalPosition()},
                                "../output/"+mesh_file+"_trajectory_agent_"+to_string(n_agent)+"_step"+to_string(i+1)+".vtk",
                                        agents[n_agent].radius);
            utils::writeFacesToVTK("../output/"+mesh_file+"_neighbourhood_agent_"+to_string(n_agent)+"_step_"+to_string(i+1)+".vtk", faces, &space);
        }


    }

    ofstream out("../output/"+mesh_file+"_manifold.vtk");
    mfem_mesh.PrintVTK(out);
    out.close();

    return 0;
}