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
draw_zbuf_trag(ZBuffer &zbuffer, img::EasyImage &image, Vector3D const &A, Vector3D const &B, Vector3D const &C,
               double d,
               double dx, double dy, img::Color color) {
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

    Vector3D P = A;
    Vector3D Q = B;
    if ((-P.y) * (-Q.y) <= 0 && P.y != Q.y) {
        double x_I = Q.x + (P.x - Q.x) * (-Q.y) / (P.y - Q.y);
        x_L_AB = x_I;
        x_R_AB = x_I;
    }
    P = A;
    Q = C;
    if ((-P.y) * (-Q.y) <= 0 && P.y != Q.y) {
        double x_I = Q.x + (P.x - Q.x) * (-Q.y) / (P.y - Q.y);
        x_L_AC = x_I;
        x_R_AC = x_I;
    }
    P = B;
    Q = C;
    if ((-P.y) * (-Q.y) <= 0 && P.y != Q.y) {
        double x_I = Q.x + (P.x - Q.x) * (-Q.y) / (P.y - Q.y);
        x_L_BC = x_I;
        x_R_BC = x_I;
    }

    int x_L = round(std::min({x_L_AB, x_L_AC, x_L_BC}) + 0.5);
    int x_R = round(std::max({x_R_AB, x_R_AC, x_R_BC}) - 0.5);

    // Calculate the 1/z values
    double x_G = (AProjected.x + BProjected.x + CProjected.x) / 3;
    double y_G = (AProjected.y + BProjected.y + CProjected.y) / 3;
    double z_GReciprocal = (1 / A.z + 1 / B.z + 1 / C.z) / 3;

    // Determine dzdx and dzdy values
    Vector3D u = B - A;
    Vector3D v = C - A;
    // scalar product u*v
    Vector3D w = Vector3D::vector(u.y * u.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x);

    double k = w.x * A.x + w.y * A.y + w.z * A.z;

    double dzdx = -w.x / k * d;
    double dzdy = -w.y / k * d;

    // zbuffer
    for (int y = y_min; y <= y_max; y++) {
        for (int x = x_L; x <= x_R; x++) {
            double zReciprocal = z_GReciprocal + (x - x_G) * dzdx + (y - y_G) * dzdy;
            if (zReciprocal < zbuffer[x][y]) {
                (image)(x, y) = color;
                zbuffer[x][y] = zReciprocal;
            }
        }
    }
}