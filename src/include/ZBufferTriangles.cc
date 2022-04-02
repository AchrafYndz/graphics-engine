#include "ZBufferTriangles.h"
#include "Logic3D.h"

std::vector<Face> triangulate(const Face &face) {
    std::vector<Face> triangles;
    for (int i=1; i<=face.point_indexes.size(); i++) {
        triangles.emplace_back(std::vector<int> {0, i, i+1});
    }
    return triangles;
}

void
draw_zbuf_trag(ZBuffer &zbuffer, img::EasyImage &img, Vector3D const &A, Vector3D const &B, Vector3D const &C, double d,
               double dx, double dy, Color color) {
    // Project the triangle
    Point2D AProjected = doProjection(A, d);
    Point2D BProjected = doProjection(B, d);
    Point2D CProjected = doProjection(C, d);
    // Determine which pixels belong to the triangle
    int y_min = ;
    int y_max = ;
    // Calculate the 1/z value
}