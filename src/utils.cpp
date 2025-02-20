//
// Created by spiccioni on 2/20/25.
//
#include "../include/utils.h"

void utils::saveAgentTrajectory(std::vector<geometrycentral::Vector3> trajectory, std::string filename) {
    std::ofstream vtk_file(filename);
    if (!vtk_file)
    {
        std::cerr << "Cannot open " << filename << " for writing." << std::endl;
        return;
    }

    int num_agents = trajectory.size();

    // Write the header for a legacy VTK file.
    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << "Agent positions\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET POLYDATA\n";

    // Write the POINTS section.
    vtk_file << "POINTS " << num_agents << " float\n";
    for (const geometrycentral::Vector3 pos: trajectory)
    {
        vtk_file << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
    }

    // Optionally, define each point as a vertex cell.
    vtk_file << "VERTICES " << num_agents << " " << (2*num_agents) << "\n";
    for (int i = 0; i < num_agents; i++)
    {
        vtk_file << "1 " << i << "\n";
    }
    vtk_file.close();

}
