//
// Created by spiccioni on 2/20/25.
//
#include "Space.h"

#ifndef UTILS_H
#define UTILS_H

namespace utils {
    void saveAgentTrajectory(std::vector<geometrycentral::Vector3> trajectory, std::string filename);
    void saveAgentTrajectoryWithRadius(const std::vector<geometrycentral::Vector3>& trajectory,
                                   const std::string& filename,
                                   float agentRadius);
}
#endif //UTILS_H
