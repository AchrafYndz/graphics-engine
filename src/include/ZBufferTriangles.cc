#include <cmath>
#include <algorithm>
#include <limits>

#include "ZBufferTriangles.h"
#include "Logic3D.h"

std::vector<Face> triangulate(const Face &face) {
    std::vector<Face> triangles;
    for (int i = 1; i <= face.point_indexes.size(); i++) {
        triangles.emplace_back(std::vector<int>{0, i, i + 1});
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
    int y_min = round(std::min({AProjected.y, BProjected.y, CProjected.y}) + 0.5);
    int y_max = round(std::max({AProjected.y, BProjected.y, CProjected.y}) - 0.5);

    int x_L_AB = std::numeric_limits<int>::max();
    int x_L_AC = std::numeric_limits<int>::max();
    int x_L_BC = std::numeric_limits<int>::max();

    int x_R_AB = -std::numeric_limits<int>::max();
    int x_R_AC = -std::numeric_limits<int>::max();
    int x_R_BC = -std::numeric_limits<int>::max();

    for (int y = y_min; y <= y_max; y++) {
        int x_L = round(std::min({x_L_AB, x_L_AC, x_L_BC}) + 0.5);
        int x_R = round(std::max({x_R_AB, x_R_AC, x_R_BC}) - 0.5);
    }
    // Calculate the 1/z value
    double x_G = (AProjected.x + BProjected.x + CProjected.x)/3;
    double y_G = (AProjected.y + BProjected.y + CProjected.y)/3;
    double z_GInverse = 1/3*A.z + 1/3*B.z + 1/3*C.z;


}