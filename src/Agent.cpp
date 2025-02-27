//
// Created by spiccioni on 2/20/25.
//
#include "../include/Agent.h"

using namespace std;
using namespace geometrycentral;
using namespace geometrycentral::surface;

Agent::Agent(Space *s, int element_id, std::vector<double> local_barycentrinc_coords, double radius, int agent_index):
    space(s) {
    element_id_CM=element_id;
    this->radius=radius;
    this->agent_index=agent_index;

    //setting position on the GC mesh:
    // TODO: check that they use the same element indexes
    geometrycentral::surface::Face element= s->gc_mesh->face(element_id);
    gc_position = geometrycentral::surface::SurfacePoint(element, geometrycentral::Vector3{local_barycentrinc_coords[0], local_barycentrinc_coords[1], local_barycentrinc_coords[2]});
}

int Agent::move(double dt) {

    geometrycentral::surface::TraceGeodesicResult result =traceGeodesic(*space->gc_geometry, gc_position, gc_local_velocity_direction*speed*dt, options);
    gc_position = result.endPoint;
    gc_local_velocity_direction = result.endingDir;

    //TODO: if you have boundary conditions you can put them here:
    if (result.hitBoundary) {
        double remaining_length = speed*dt - result.length;
        Vector3 global_position = getGlobalPosition();

        // save the velocity direction before the corrective-shift:
        Vector3 global_vel_direction;
        space->convertLocaToGlobalVector(gc_position, gc_local_velocity_direction, global_vel_direction);

        {
            // BC: x-direction periodic boundary conditions:
            // NOTE: this is a quick and dirty solution that doesn't require much coding
            //       if your mesh has some more complicated structures, you can think of creating a table containing element-to-element
            //       pairings to encode the periodicity, otherwise you could also perform the shift in global coordinates and then loop
            //       over the elements to check on which one we are...
            /*
            double periodic_shift=0.0;
            if (global_position[0] >= 10.0) {
                periodic_shift = -10.0;
            }else if (global_position[0] <= 0.0) {
                periodic_shift = 10.0;
            }
            if (periodic_shift != 0.0) {
                // you have to find the new element and the new local coordinates:
                Vector3 shift{periodic_shift,0,0};
                Vector2 local_shift;

                space->convertGlobalToLocalVector(gc_position, shift, local_shift);
                result = traceGeodesic(*space->gc_geometry, gc_position, local_shift, options);
                gc_position = result.endPoint;

                space->convertGlobalToLocalVector(gc_position, global_vel_direction, gc_local_velocity_direction);
            }
            */
        }

        // REFLECTIVE BOUNDARY CONDITIONS:
        // BC:x-direction reflective boundary conditions:
        if (global_position[0] >= 10.0 || global_position[0] <= 0.0) {
            global_vel_direction[0] *= -1.0;
            space->convertGlobalToLocalVector(gc_position, global_vel_direction, gc_local_velocity_direction);
        }

        // BC: y-direction reflective boundary conditions:
        if (global_position[1] >= 10.0 || global_position[1] <= 0.0) {
            global_vel_direction[1] *= -1.0;
            space->convertGlobalToLocalVector(gc_position, global_vel_direction, gc_local_velocity_direction);
        }
        
        // move for how much is left:
        result = traceGeodesic(*space->gc_geometry, gc_position, gc_local_velocity_direction*remaining_length, options);
        gc_position = result.endPoint;
        gc_local_velocity_direction = result.endingDir;
    }

    return 0;
}

std::unordered_set<Face> Agent::findFacesWithinRadius(SurfacePoint sp, double radius, Space *space) {
    std::unordered_set<Face> withinRadiusFaces;

    // Priority queue for Dijkstra's algorithm
    std::priority_queue<std::pair<double, Vertex>, std::vector<std::pair<double, Vertex>>, std::greater<>> pq;

    // Distance map with custom hash function
    std::unordered_map<Vertex, double, utils::VertexHash> distance;

    // Initialize with the closest vertex to SurfacePoint
    Vertex startVertex = sp.nearestVertex();
    pq.push({0.0, startVertex});
    distance[startVertex] = 0.0;

    int iter_counter = 0;
    while (!pq.empty() && iter_counter < utils::MAX_NUMBER_OF_MESH_SEARCH_ITERS) {
        auto [dist, v] = pq.top();
        pq.pop();

        // If the distance exceeds the radius, stop processing
        if (dist > radius) continue;

        // Add adjacent faces of the vertex
        for (Face f : v.adjacentFaces()) {
            withinRadiusFaces.insert(f);
        }

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
        iter_counter++;
    }
    return withinRadiusFaces;
}

void Agent::rotateVelocityDirection(double angle) {
    gc_local_velocity_direction = gc_local_velocity_direction.rotate(angle);
}

void Agent::setLocalVelocity(std::vector<double> velocity) {
    double norm = 0.0;
    for (int d =0;d<velocity.size();d++) {
        norm += velocity[d]*velocity[d];
        gc_local_velocity_direction[d] = velocity[d];
    }
    norm = sqrt(norm);

    // update the private variables:
    speed = norm;
    gc_local_velocity_direction = gc_local_velocity_direction/norm;
}

geometrycentral::Vector2 Agent::getLocalVelocity() {
    return gc_local_velocity_direction*speed;
}

geometrycentral::Vector3 Agent::getGlobalPosition() {
    return gc_position.interpolate(space->gc_geometry->inputVertexPositions);
}

