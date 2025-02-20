//
// Created by spiccioni on 2/20/25.
//

#include "Space.h"

#ifndef AGENT_H
#define AGENT_H

class Agent {
public:
    Agent(Space *s, int element_id, std::vector<double> local_barycentrinc_coords, double radius);

    // methods:
    int move(double dt);
    void rotateVelocityDirection(double angle);
    void setLocalVelocity(std::vector<double> velocity);
    geometrycentral::Vector2 getLocalVelocity();
    geometrycentral::Vector3 getGlobalPosition();


    geometrycentral::surface::SurfacePoint gc_position;

private:
    // NOTE: CM stands for center of mass (CENTER)
    int element_id_CM=-1;
    double speed = 0.0;
    double radius = 0.0;
    Space *space;

    // GC variables:
    geometrycentral::Vector2 gc_local_velocity_direction;
    geometrycentral::surface::TraceOptions options = geometrycentral::surface::defaultTraceOptions;
};
#endif //AGENT_H
