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
    void writeFacesToVTK(const std::string& filename, const std::unordered_set<geometrycentral::surface::Face>& faces, Space* space);

    struct VertexHash {
        std::size_t operator()(const geometrycentral::surface::Vertex& v) const {
            return std::hash<size_t>()(v.getIndex());
        }
    };
    const int MAX_NUMBER_OF_MESH_SEARCH_ITERS = 1e5;
}
#endif //UTILS_H
