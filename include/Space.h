//
// Created by spiccioni on 2/20/25.
//
#include <iostream>
#include <mfem.hpp>
#include <geometrycentral/numerical/linear_algebra_types.h>
#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <geometrycentral/surface/trace_geodesic.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/surface_point.h>

#ifndef SPACE_H
#define SPACE_H

class Space {
public:
    Space(mfem::Mesh &mesh);

    void saveToVTK(const std::string &filename, mfem::GridFunction &field,
                                          const std::string &field_name);

    void setMFEMMesh(mfem::Mesh &mesh);
    void loadGCMeshFromFile(std::string filename);

    mfem::FiniteElementCollection* finite_element_collection;    // Type of finite elements (H1, DG, etc.)
    mfem::FiniteElementSpace* finite_element_space;

    std::unique_ptr<geometrycentral::surface::ManifoldSurfaceMesh> gc_mesh;
    std::unique_ptr<geometrycentral::surface::VertexPositionGeometry> gc_geometry;
private:
    //MFEM Mesh:
    int manifold_dimension;
    int space_dimension;
    int order = 1;

    mfem::Mesh* mfem_mesh;
};

#endif //SPACE_H
