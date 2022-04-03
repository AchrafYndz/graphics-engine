#define _USE_MATH_DEFINES // Windows only

#include "Figures3D.h"
#include "ZBufferTriangles.h"

#include <cmath>

Figure
createCube(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, bool toTriangulate) {
    // Create points
    std::vector<Vector3D> points;
    points.push_back(Vector3D::point(1, -1, -1));
    points.push_back(Vector3D::point(-1, 1, -1));
    points.push_back(Vector3D::point(1, 1, 1));
    points.push_back(Vector3D::point(-1, -1, 1));
    points.push_back(Vector3D::point(1, 1, -1));
    points.push_back(Vector3D::point(-1, -1, -1));
    points.push_back(Vector3D::point(1, -1, 1));
    points.push_back(Vector3D::point(-1, 1, 1));

    // Create faces
    std::vector<Face> faces;
    faces.emplace_back(std::vector<int>{0, 4, 2, 6});
    faces.emplace_back(std::vector<int>{4, 1, 7, 2});
    faces.emplace_back(std::vector<int>{1, 5, 3, 7});
    faces.emplace_back(std::vector<int>{5, 0, 6, 3});
    faces.emplace_back(std::vector<int>{6, 2, 7, 3});
    faces.emplace_back(std::vector<int>{0, 5, 1, 4});

    // Triangulation
    if (toTriangulate) {
        std::vector<Face> triangles;
        for (const Face &face: faces) {
            std::vector<Face> newTriangles = triangulate(face);
            triangles.insert(triangles.end(), newTriangles.begin(), newTriangles.end());
        }
        faces = triangles;
    }

    // Create figure
    return {points, faces, color, center, scale, angleX, angleY, angleZ};
}

Figure createTetrahedron(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, bool toTriangulate) {
    // Create points
    std::vector<Vector3D> points;
    points.push_back(Vector3D::point(1, -1, -1));
    points.push_back(Vector3D::point(-1, 1, -1));
    points.push_back(Vector3D::point(1, 1, 1));
    points.push_back(Vector3D::point(-1, -1, 1));
    // Create faces
    std::vector<Face> faces;
    faces.emplace_back(std::vector<int>{0, 1, 2});
    faces.emplace_back(std::vector<int>{1, 3, 2});
    faces.emplace_back(std::vector<int>{0, 3, 1});
    faces.emplace_back(std::vector<int>{0, 2, 3});

    // Triangulation
    if (toTriangulate) {
        std::vector<Face> triangles;
        for (const Face &face: faces) {
            std::vector<Face> newTriangles = triangulate(face);
            triangles.insert(triangles.end(), newTriangles.begin(), newTriangles.end());
        }
        faces = triangles;
    }

    // Create figure
    return {points, faces, color, center, scale, angleX, angleY, angleZ};
}

Figure createOctahedron(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, bool toTriangulate) {
    // Create points
    std::vector<Vector3D> points;
    points.push_back(Vector3D::point(1, 0, 0));
    points.push_back(Vector3D::point(0, 1, 0));
    points.push_back(Vector3D::point(-1, 0, 0));
    points.push_back(Vector3D::point(0, -1, 0));
    points.push_back(Vector3D::point(0, 0, -1));
    points.push_back(Vector3D::point(0, 0, 1));

    // Create faces
    std::vector<Face> faces;
    faces.emplace_back(std::vector<int>{0, 1, 5});
    faces.emplace_back(std::vector<int>{1, 2, 5});
    faces.emplace_back(std::vector<int>{2, 3, 5});
    faces.emplace_back(std::vector<int>{3, 0, 5});
    faces.emplace_back(std::vector<int>{1, 0, 4});
    faces.emplace_back(std::vector<int>{2, 1, 4});
    faces.emplace_back(std::vector<int>{3, 2, 4});
    faces.emplace_back(std::vector<int>{0, 3, 4});

    // Triangulation
    if (toTriangulate) {
        std::vector<Face> triangles;
        for (const Face &face: faces) {
            std::vector<Face> newTriangles = triangulate(face);
            triangles.insert(triangles.end(), newTriangles.begin(), newTriangles.end());
        }
        faces = triangles;
    }

    // Create figure
    return {points, faces, color, center, scale, angleX, angleY, angleZ};
}

Figure createIcosahedron(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, bool toTriangulate) {
    // Create points
    std::vector<Vector3D> points;

    for (int i = 1; i <= 12; i++) {
        if (i == 1) points.push_back(Vector3D::point(0, 0, sqrt(5) / 2));
        else if (2 <= i && i <= 6)
            points.push_back(Vector3D::point(cos((i - 2) * 2 * M_PI / 5), sin((i - 2) * 2 * M_PI / 5), 0.5));
        else if (7 <= i && i <= 11)
            points.push_back(
                    Vector3D::point(cos(M_PI / 5 + (i - 7) * 2 * M_PI / 5), sin(M_PI / 5 + (i - 7) * 2 * M_PI / 5),
                                    -0.5));
        else points.push_back(Vector3D::point(0, 0, -sqrt(5) / 2));
    }

    // Create faces
    std::vector<Face> faces;
    std::vector<int> point_indexes0 = {0, 1, 2};
    faces.emplace_back(std::vector<int>{0, 1, 2});
    faces.emplace_back(std::vector<int>{0, 2, 3});
    faces.emplace_back(std::vector<int>{0, 3, 4});
    faces.emplace_back(std::vector<int>{0, 4, 5});
    faces.emplace_back(std::vector<int>{0, 5, 1});
    faces.emplace_back(std::vector<int>{1, 6, 2});
    faces.emplace_back(std::vector<int>{2, 6, 7});
    faces.emplace_back(std::vector<int>{2, 7, 3});
    faces.emplace_back(std::vector<int>{3, 7, 8});
    faces.emplace_back(std::vector<int>{3, 8, 4});
    faces.emplace_back(std::vector<int>{4, 8, 9});
    faces.emplace_back(std::vector<int>{4, 9, 5});
    faces.emplace_back(std::vector<int>{5, 9, 10});
    faces.emplace_back(std::vector<int>{5, 10, 1});
    faces.emplace_back(std::vector<int>{1, 10, 6});
    faces.emplace_back(std::vector<int>{11, 7, 6});
    faces.emplace_back(std::vector<int>{11, 8, 7});
    faces.emplace_back(std::vector<int>{11, 9, 8});
    faces.emplace_back(std::vector<int>{11, 10, 9});
    faces.emplace_back(std::vector<int>{11, 6, 10});

    // Triangulation
    if (toTriangulate) {
        std::vector<Face> triangles;
        for (const Face &face: faces) {
            std::vector<Face> newTriangles = triangulate(face);
            triangles.insert(triangles.end(), newTriangles.begin(), newTriangles.end());
        }
        faces = triangles;
    }

    // Create figure
    return {points, faces, color, center, scale, angleX, angleY, angleZ};
}

Figure createDodecahedron(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, bool toTriangulate) {
    // Create points
    std::vector<Vector3D> points;

    Figure icosahedron = createIcosahedron(color, center, scale, angleX, angleY, angleZ, toTriangulate);

    for (const Face &triangle: icosahedron.faces) {
        Vector3D point1 = icosahedron.points[triangle.point_indexes[0]];
        Vector3D point2 = icosahedron.points[triangle.point_indexes[1]];
        Vector3D point3 = icosahedron.points[triangle.point_indexes[2]];
        points.push_back((point1 + point2 + point3) / 3);
    }

    // Create faces
    std::vector<Face> faces;

    faces.emplace_back(std::vector<int>{0, 1, 2, 3, 4});
    faces.emplace_back(std::vector<int>{0, 5, 6, 7, 1});
    faces.emplace_back(std::vector<int>{1, 7, 8, 9, 2});
    faces.emplace_back(std::vector<int>{2, 9, 10, 11, 3});
    faces.emplace_back(std::vector<int>{3, 11, 12, 13, 4});
    faces.emplace_back(std::vector<int>{4, 13, 14, 5, 0});
    faces.emplace_back(std::vector<int>{19, 18, 17, 16, 15});
    faces.emplace_back(std::vector<int>{19, 14, 13, 12, 18});
    faces.emplace_back(std::vector<int>{18, 12, 11, 10, 17});
    faces.emplace_back(std::vector<int>{17, 10, 9, 8, 16});
    faces.emplace_back(std::vector<int>{16, 8, 7, 6, 15});
    faces.emplace_back(std::vector<int>{15, 6, 5, 14, 19});

    // Triangulation
    if (toTriangulate) {
        std::vector<Face> triangles;
        for (const Face &face: faces) {
            std::vector<Face> newTriangles = triangulate(face);
            triangles.insert(triangles.end(), newTriangles.begin(), newTriangles.end());
        }
        faces = triangles;
    }

    // Create figure
    return {points, faces, color, center, scale, angleX, angleY, angleZ};
}

Figure
createSphere(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, const int n, bool toTriangulate) {
    Figure icosahedron = createIcosahedron(color, center, scale, angleX, angleY, angleZ, toTriangulate);

    for (int _ = 0; _ < n; _++) {
        std::vector<Face> faces;
        std::vector<Vector3D> points;
        int triangleCounter = 0;
        for (Face triangle: icosahedron.faces) {
            // Fetch A, B and C
            Vector3D A = icosahedron.points[triangle.point_indexes[0]];
            Vector3D B = icosahedron.points[triangle.point_indexes[1]];
            Vector3D C = icosahedron.points[triangle.point_indexes[2]];

            // Calculate D, E and F
            Vector3D D = (A + B) / 2;
            Vector3D E = (A + C) / 2;
            Vector3D F = (B + C) / 2;

            // Add new points to temp points vector
            points.push_back(A); // 0 + 6*triangleCounter
            points.push_back(B); // 1 + 6*triangleCounter
            points.push_back(C); // 2 + 6*triangleCounter
            points.push_back(D); // 3 + 6*triangleCounter
            points.push_back(E); // 4 + 6*triangleCounter
            points.push_back(F); // 5 + 6*triangleCounter

            // Add new faces to temp faces vector
            faces.emplace_back(
                    std::vector<int>{6 * triangleCounter, 3 + 6 * triangleCounter, 4 + 6 * triangleCounter}); // ADE
            faces.emplace_back(
                    std::vector<int>{1 + 6 * triangleCounter, 5 + 6 * triangleCounter, 3 + 6 * triangleCounter}); // BFD
            faces.emplace_back(
                    std::vector<int>{2 + 6 * triangleCounter, 4 + 6 * triangleCounter, 5 + 6 * triangleCounter}); // CEF
            faces.emplace_back(
                    std::vector<int>{3 + 6 * triangleCounter, 5 + 6 * triangleCounter, 4 + 6 * triangleCounter}); // DFE
            triangleCounter++;
        }
        // Replace points and faces
        icosahedron.points = points;
        // Triangulation
        if (toTriangulate) {
            std::vector<Face> triangles;
            for (const Face &face: faces) {
                std::vector<Face> newTriangles = triangulate(face);
                triangles.insert(triangles.end(), newTriangles.begin(), newTriangles.end());
            }
            icosahedron.faces = triangles;
        }
        else icosahedron.faces = faces;

    }

    for (Vector3D &point: icosahedron.points) {
        // Calculate radius
        double r = sqrt(pow(point.x, 2) + pow(point.y, 2) + pow(point.z, 2));

        // Replace coordinates
        point /= r;
    }

    return icosahedron;
}

Figure createCone(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, const int n,
                  const double h, bool toTriangulate) {
    // Create points
    std::vector<Vector3D> points;
    for (int i = 0; i <= n; i++) {
        if (i == n) points.push_back(Vector3D::point(0, 0, h));
        else points.push_back(Vector3D::point(cos(2 * i * M_PI / n), sin(2 * i * M_PI / n), 0));
    }

    // Create faces
    std::vector<Face> faces;
    for (int i = 0; i <= n; i++) {
        if (i != n) faces.emplace_back(std::vector<int>{i, (i + 1) % n, n});
        else {
            std::vector<int> point_indexes;
            for (int j = n - 1; j >= 0; j--) point_indexes.push_back(j);
            faces.emplace_back(point_indexes);
        }
    }

    // Triangulation
    if (toTriangulate) {
        std::vector<Face> triangles;
        for (const Face &face: faces) {
            std::vector<Face> newTriangles = triangulate(face);
            triangles.insert(triangles.end(), newTriangles.begin(), newTriangles.end());
        }
        faces = triangles;
    }

    return {points, faces, color, center, scale, angleX, angleY, angleZ};
}

Figure
createCylinder(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, const int n,
               const double h, bool toTriangulate) {
    // Create points
    std::vector<Vector3D> points;

    points.reserve(2 * n);
    // Bottom surface
    for (int i = 0; i < n; i++) points.push_back(Vector3D::point(cos(2 * i * M_PI / n), sin(2 * i * M_PI / n), 0));
    // Top surface
    for (int i = n; i < 2 * n; i++) points.push_back(Vector3D::point(cos(2 * i * M_PI / n), sin(2 * i * M_PI / n), h));

    // Create faces
    std::vector<Face> faces;

    faces.reserve(n + 2);
    // draw bottom and top surface outside for loop
    for (int i = 0; i <= n; i++) {
        if (i == 0) {
            // bottom surface
            std::vector<int> point_indexes_bottom;
            for (int j = n - 1; j >= 0; j--) point_indexes_bottom.push_back(j);
            faces.emplace_back(point_indexes_bottom);
        }
        if (i == n) {
            // top surface
            std::vector<int> point_indexes_top;
            for (int j = 2 * n - 1; j >= n; j--) point_indexes_top.push_back(j);
            faces.emplace_back(point_indexes_top);
        } else if (i == n - 1) faces.emplace_back(std::vector<int>{i + 1, n + i, i, (i + 1) % n});
        else faces.emplace_back(std::vector<int>{n + i + 1, n + i, i, i + 1});
    }

    // Triangulation
    if (toTriangulate) {
        std::vector<Face> triangles;
        for (const Face &face: faces) {
            std::vector<Face> newTriangles = triangulate(face);
            triangles.insert(triangles.end(), newTriangles.begin(), newTriangles.end());
        }
        faces = triangles;
    }

    return {points, faces, color, center, scale, angleX, angleY, angleZ};
}

Figure
createTorus(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, const double r,
            const double R, const int n, const int m, bool toTriangulate) {
    // Create points
    std::vector<Vector3D> points;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            // Calculate u and v
            double u = 2 * i * M_PI / n;
            double v = 2 * j * M_PI / m;

            // Determine x, y and z using the parametric equations
            double x_uv = (R + r * cos(v)) * cos(u);
            double y_uv = (R + r * cos(v)) * sin(u);
            double z_uv = r * sin(v);

            points.push_back(Vector3D::point(x_uv, y_uv, z_uv));
        }
    }

    // Create faces
    std::vector<Face> faces;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            faces.emplace_back(
                    std::vector<int>{m * i + j, ((i + 1) % n) * m + j, ((i + 1) % n) * m + ((j + 1) % m),
                                     i * m + ((j + 1) % m)});
        }
    }

    // Triangulation
    if (toTriangulate) {
        std::vector<Face> triangles;
        for (const Face &face: faces) {
            std::vector<Face> newTriangles = triangulate(face);
            triangles.insert(triangles.end(), newTriangles.begin(), newTriangles.end());
        }
        faces = triangles;
    }

    return {points, faces, color, center, scale, angleX, angleY, angleZ};
}


