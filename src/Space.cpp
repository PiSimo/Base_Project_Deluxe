//
// Created by spiccioni on 2/20/25.
//
#include "../include/Space.h"

using namespace mfem;
using namespace std;


Space::Space(mfem::Mesh &mesh) {
    setMFEMMesh(mesh);
}

void Space::setMFEMMesh(mfem::Mesh &mesh) {
    // update the mesh (or define)
    mfem_mesh = &mesh;

    // line=1,disk=2,...
    manifold_dimension = mfem_mesh->Dimension();
    space_dimension = mfem_mesh->SpaceDimension();

    // define the functional spaces
    finite_element_collection = new mfem::H1_FECollection(order, manifold_dimension); // type of functional space to use (default H1 of order 1)
    finite_element_space = new mfem::FiniteElementSpace(mfem_mesh, finite_element_collection); // generate the basis/test functions for our specific mesh
}

void Space::loadGCMeshFromFile(std::string filename) {
    std::tie(gc_mesh, gc_geometry) = geometrycentral::surface::readManifoldSurfaceMesh(filename);
}

void Space::saveToVTK(const std::string &filename, mfem::GridFunction &field,
                      const std::string &field_name) {
    std::ofstream ofs(filename);
    // First write the geometry (points and connectivity)
    mfem_mesh->PrintVTK(ofs, 0, false);
    // Then append the field data (the actual solution values)
    field.SaveVTK(ofs, field_name, 0);
    ofs.close();
}
