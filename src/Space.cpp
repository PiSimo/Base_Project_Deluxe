//
// Created by spiccioni on 2/20/25.
//
#include "../include/Space.h"

using namespace mfem;
using namespace std;

using namespace geometrycentral;
using namespace geometrycentral::surface;

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
    std::tie(gc_mesh, gc_geometry) = readManifoldSurfaceMesh(filename);
}

std::vector<Vector3> Space::gc_getSurfaceTangentBasis(SurfacePoint surface_point, bool include_normal=false) {
    Vector3 e1, e2, normal; // Local tangent space basis

    if (surface_point.type == SurfacePointType::Vertex) {
        // Vertex case
        geometrycentral::surface::Vertex v = surface_point.vertex;
        normal = gc_geometry->vertexNormals[v];

        // Compute a local tangent basis using adjacent edges
        for (Halfedge he : v.outgoingHalfedges()) {
            Vector3 candidate = (gc_geometry->inputVertexPositions[he.twin().vertex()] -
                                 gc_geometry->inputVertexPositions[v]).normalize();
            if (fabs(dot(candidate, normal)) < 1e-6) {  // Ensure perpendicularity
                e1 = candidate;
                break;
            }
        }
        e2 = cross(normal, e1).normalize();

    } else if (surface_point.type == SurfacePointType::Edge) {
        // Edge case
        Edge e = surface_point.edge;
        Halfedge he = e.halfedge();
        Vector3 v0 = gc_geometry->inputVertexPositions[he.vertex()];
        Vector3 v1 = gc_geometry->inputVertexPositions[he.twin().vertex()];

        e1 = (v1 - v0).normalize();  // Tangent along the edge

        // Compute normal from one of the adjacent faces
        Face f = he.face();
        normal = f.isBoundaryLoop() ? Vector3{0, 0, 0} : gc_geometry->faceNormal(f);

        e2 = cross(normal, e1).normalize();  // Perpendicular to the edge

    } else if (surface_point.type == SurfacePointType::Face) {
        // Face case
        Face f = surface_point.face;
        Halfedge he = f.halfedge();
        Vector3 v0 = gc_geometry->inputVertexPositions[he.vertex()];
        Vector3 v1 = gc_geometry->inputVertexPositions[he.next().vertex()];
        Vector3 v2 = gc_geometry->inputVertexPositions[he.next().next().vertex()];

        e1 = (v1 - v0).normalize();
        normal = gc_geometry->faceNormal(f);
        e2 = cross(normal, e1).normalize();
    } else {
        throw std::runtime_error("Invalid SurfacePoint type.");
    }


    if (include_normal) {
        return {e1,e2,normal};
    }

    return {e1, e2};
}


void Space::convertGlobalToLocalVector(SurfacePoint s_point, Vector3 global_vec, Vector2 &local_v) {
    std::vector<Vector3> tangent_basis = gc_getSurfaceTangentBasis(s_point);

    // Project the global vector onto the local basis
    double localX = dot(global_vec, tangent_basis[0]);
    double localY = dot(global_vec, tangent_basis[1]);

    local_v.x = localX;
    local_v.y = localY;
}

void Space::convertLocaToGlobalVector(SurfacePoint s_point, Vector2 local_vec, Vector3 &global_vec) {
    std::vector<Vector3> tangent_basis = gc_getSurfaceTangentBasis(s_point);

    global_vec = local_vec.x*tangent_basis[0] + local_vec.y*tangent_basis[1];

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
