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

void utils::saveAgentTrajectoryWithRadius(const std::vector<geometrycentral::Vector3>& trajectory,
                                   const std::string& filename,
                                   float agentRadius)
{
    int nSegments = 16;
    std::ofstream vtk_file(filename);
    if (!vtk_file) {
        std::cerr << "Cannot open " << filename << " for writing.\n";
        return;
    }

    // Number of agents
    int numAgents = static_cast<int>(trajectory.size());

    // For each agent, we store (nSegments+1) points:
    //   1 center + nSegments around the perimeter
    // Therefore total points:
    int totalPoints = numAgents * (nSegments + 1);

    // Each disk is made of nSegments triangles
    // So total polygons = numAgents * nSegments
    // Each triangle is stored in the VTK POLYGONS section as:
    //   3 vertexIndex0 vertexIndex1 vertexIndex2  (4 integers per triangle)
    int totalPolygons = numAgents * nSegments;
    // The connectivity array has 4* totalPolygons entries:
    //   (3 vertices + 1 count) per polygon
    int totalConnectivityEntries = 4 * totalPolygons;

    // Write header for legacy VTK
    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << "Agent disks\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET POLYDATA\n";

    // Write POINTS
    vtk_file << "POINTS " << totalPoints << " float\n";

    // We will store the (x,y,z) coordinates of all disk vertices
    // Keep track of them so we can refer to indices later
    //   Alternatively, we can write them directly, but we need the same ordering
    //   for building the polygons section.

    // We'll build them in memory, but you can also stream them out directly.
    std::vector<geometrycentral::Vector3> allPoints;
    allPoints.reserve(totalPoints);

    // Precompute the angle step for the circle
    float dTheta = 2.0f * static_cast<float>(M_PI) / static_cast<float>(nSegments);

    // Generate all points
    for (int i = 0; i < numAgents; i++) {
        // The center
        const auto& center = trajectory[i];
        // Add center point
        allPoints.push_back(center);

        // Add perimeter points
        for (int seg = 0; seg < nSegments; seg++) {
            float angle = seg * dTheta;
            float xOff = agentRadius * std::cos(angle);
            float yOff = agentRadius * std::sin(angle);
            // Place the circle in XY-plane around "center"
            // If you want to orient it differently, you'd apply a rotation or normal-based transform here.
            geometrycentral::Vector3 perimeterPos(center[0] + xOff,
                                                 center[1] + yOff,
                                                 center[2]);
            allPoints.push_back(perimeterPos);
        }
    }

    // Now write all the points to the file
    for (const auto& pt : allPoints) {
        vtk_file << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
    }

    // Write the POLYGONS section
    vtk_file << "POLYGONS " << totalPolygons << " " << totalConnectivityEntries << "\n";

    // For each agent's disk, we have nSegments triangles
    // Index offset for each agent's set of points
    //   #points per agent = (nSegments + 1)
    for (int i = 0; i < numAgents; i++) {
        int startIndex = i * (nSegments + 1);

        int centerIndex = startIndex;
        int firstPerimeterIndex = startIndex + 1;

        // For each segment, define a triangle: center, perimeter[j], perimeter[(j+1)%nSegments]
        for (int seg = 0; seg < nSegments; seg++) {
            int p0 = firstPerimeterIndex + seg;
            int p1 = firstPerimeterIndex + ((seg + 1) % nSegments);

            vtk_file << "3 " << centerIndex << " " << p0 << " " << p1 << "\n";
        }
    }

    vtk_file.close();
}

