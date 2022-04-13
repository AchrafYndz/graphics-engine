#include <cmath>
#include <algorithm>
#include <limits>
#include <cassert>

#include "ZBufferTriangles.h"

std::vector<Face> triangulate(const Face &face) {
    std::vector<Face> triangles;
    for (int i = 1; i <= face.point_indexes.size() - 2; i++) {
        triangles.emplace_back(
                std::vector<int>{face.point_indexes[0], face.point_indexes[i], face.point_indexes[i + 1]});
    }
    return triangles;
}

void
draw_zbuf_trag(ZBuffer &zbuffer, img::EasyImage &image, Vector3D const &A, Vector3D const &B, Vector3D const &C,
               double d,
               double dx, double dy, img::Color color) {
    // Project the triangle
    Point2D AProjected = doProjection(A, d, dx, dy);
    Point2D BProjected = doProjection(B, d, dx, dy);
    Point2D CProjected = doProjection(C, d, dx, dy);

    // Calculations needed for the 1/z values
    double x_G = (AProjected.x + BProjected.x + CProjected.x) / 3;
    double y_G = (AProjected.y + BProjected.y + CProjected.y) / 3;
    double z_GReciprocal = (1 / (3 * A.z)) + (1 / (3 * B.z)) + (1 / (3 * C.z));

//     Determine dzdx and dzdy values
    Vector3D u = B - A;
    Vector3D v = C - A;
//     scalar product u*v
    Vector3D w = Vector3D::cross(u, v);

    double k = w.x * A.x + w.y * A.y + w.z * A.z;

    double dzdx = -w.x / k * d;
    double dzdy = -w.y / k * d;

    // Determine which pixels belong to the triangle
    int y_min = lround(std::min({AProjected.y, BProjected.y, CProjected.y}) + 0.5);
    int y_max = lround(std::max({AProjected.y, BProjected.y, CProjected.y}) - 0.5);

    for (int y_I = y_min; y_I <= y_max; y_I++) {
        int x_L_AB = std::numeric_limits<int>::max();
        int x_L_AC = std::numeric_limits<int>::max();
        int x_L_BC = std::numeric_limits<int>::max();

        int x_R_AB = -std::numeric_limits<int>::max();
        int x_R_AC = -std::numeric_limits<int>::max();
        int x_R_BC = -std::numeric_limits<int>::max();
        // AB
        Point2D P = AProjected;
        Point2D Q = BProjected;
        if ((y_I - P.y) * (y_I - Q.y) <= 0 && P.y != Q.y) {
            double x_I = Q.x + ((P.x - Q.x) * (y_I - Q.y) / (P.y - Q.y));
            x_L_AB = lround(x_I);
            x_R_AB = lround(x_I);
        }
        // AC
        Q = CProjected;
        if ((y_I - P.y) * (y_I - Q.y) <= 0 && P.y != Q.y) {
            double x_I = Q.x + ((P.x - Q.x) * (y_I - Q.y) / (P.y - Q.y));
            x_L_AC = lround(x_I);
            x_R_AC = lround(x_I);
        }
        // BC
        P = BProjected;
        if ((y_I - P.y) * (y_I - Q.y) <= 0 && P.y != Q.y) {
            double x_I = Q.x + ((P.x - Q.x) * (y_I - Q.y) / (P.y - Q.y));
            x_L_BC = lround(x_I);
            x_R_BC = lround(x_I);
        }
        int x_L = lround(std::min({x_L_AB, x_L_AC, x_L_BC}) + 0.5);
        int x_R = lround(std::max({x_R_AB, x_R_AC, x_R_BC}) - 0.5);
        for (int x = x_L; x <= x_R; x++) {
            double zReciprocal = (1.0001 * z_GReciprocal) + ((x - x_G) * dzdx) + ((y_I - y_G) * dzdy);
            if (zReciprocal < zbuffer[x][y_I]) {
                (image)(x, y_I) = color;
                zbuffer[x][y_I] = zReciprocal;
            }
        }
    }
}

void
getImageSpecs(const Lines2D &lines, const int size, double &d, double &dx, double &dy, double &width, double &height) {
    double xmin = lines.front().p1.x;
    double xmax = lines.front().p1.x;
    double ymin = lines.front().p2.y;
    double ymax = lines.front().p2.y;
    // determine max and min
    for (Line2D line: lines) {
        if (line.p1.x < xmin) xmin = line.p1.x;
        if (line.p2.x < xmin) xmin = line.p2.x;
        if (line.p1.y < ymin) ymin = line.p1.y;
        if (line.p2.y < ymin) ymin = line.p2.y;

        if (line.p1.x > xmax) xmax = line.p1.x;
        if (line.p2.x > xmax) xmax = line.p2.x;
        if (line.p1.y > ymax) ymax = line.p1.y;
        if (line.p2.y > ymax) ymax = line.p2.y;
    }

    double xrange = std::abs(xmax - xmin);
    double yrange = std::abs(ymax - ymin);

    double imagex = size * xrange / std::max(xrange, yrange);
    double imagey = size * yrange / std::max(xrange, yrange);

    d = 0.95 * imagex / xrange;

    double DCx = d * (xmin + xmax) / 2;
    double DCy = d * (ymin + ymax) / 2;

    dx = imagex / 2 - DCx;
    dy = imagey / 2 - DCy;

    width = lround(imagex);
    height = lround(imagey);
}

Point2D doProjection(const Vector3D point, const double d, const double dx, const double dy) {
    return {-(point.x * d / point.z) + dx, -(point.y * d / point.z) + dy};
}
