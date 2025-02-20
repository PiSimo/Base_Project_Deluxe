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

    if (result.hitBoundary) {
        //TODO: if you have boundary conditions you can put them here:
        geometrycentral::Vector3 global_position = getGlobalPosition();

        bool periodic_transformation = false;
        // BC: x-direction periodic boundary conditions:
        double shift=0.0;
        if (global_position[0] >= 10.0) {
            //global_position[0] = global_position[0] - 10.0;
            shift = -10.0;
            periodic_transformation = true;
        }else if (global_position[0] <= 0.0) {
            //global_position[0] = global_position[0] + 10.0;
            shift = 10.0;
            periodic_transformation = true;
        }
        if (periodic_transformation) {
            // you have to find the new element and the new local coordinates:
            geometrycentral::Vector2 dir{0,1};
            result =traceGeodesic(*space->gc_geometry, gc_position, dir*shift, options);
            gc_position = result.endPoint;
            gc_local_velocity_direction = result.endingDir;
        }

        // BC: y-direction reflective boundary conditions:
        if (global_position[1] >= 10.0 || global_position[1] <= 0.0) {
            gc_local_velocity_direction[1] *= -1;
        }


    }

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

