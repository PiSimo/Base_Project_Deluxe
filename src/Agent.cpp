//
// Created by spiccioni on 2/20/25.
//
#include "../include/Agent.h"

using namespace std;

Agent::Agent(Space *s, int element_id, std::vector<double> local_barycentrinc_coords, double radius):
    space(s) {
    element_id_CM=element_id;
    this->radius=radius;

    //setting position on the GC mesh:
    // TODO: check that they use the same element indexes
    geometrycentral::surface::Face element= s->gc_mesh->face(element_id);
    gc_position = geometrycentral::surface::SurfacePoint(element, geometrycentral::Vector3{local_barycentrinc_coords[0], local_barycentrinc_coords[1], local_barycentrinc_coords[2]});
}

int Agent::move(double dt) {

    geometrycentral::surface::TraceGeodesicResult result =traceGeodesic(*space->gc_geometry, gc_position, gc_local_velocity_direction*speed*dt, options);
    gc_position = result.endPoint;
    gc_local_velocity_direction = result.endingDir;

    return 0;
}

void Agent::rotateVelocityDirection(double angle) {
    double c = cos(angle);
    double s = sin(angle);

    gc_local_velocity_direction.x = c*gc_local_velocity_direction.x - s*gc_local_velocity_direction.y;
    gc_local_velocity_direction.y = s*gc_local_velocity_direction.x + c*gc_local_velocity_direction.y;
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

