#include "Fractals3D.h"
#include "Logic3D.h"
#include "Figures3D.h"

void generateFractal(Figure &fig, Figures3D &fractal, const int nr_iterations, const double scale) {
    fractal.push_back(fig);
    for (int i = 0; i < nr_iterations; i++) {
        Figures3D newFractals;
        for (Figure fracFig: fractal) {
            for (int j = 0; j < fig.points.size(); j++) {
                Figure copyFig = fracFig;
                const Matrix scaleMatrix = scaleFigure(1 / scale);
                applyTransformation(copyFig, scaleMatrix);
                Matrix translationMatrix = translate(fracFig.points[j] - copyFig.points[j]);
                applyTransformation(copyFig, translationMatrix);
                newFractals.push_back(copyFig);
            }
        }
        fractal = newFractals;
    }
}

Figure createBuckyBall(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ,
                       bool toTriangulate) {
    Figure buckyBall = createIcosahedron(color, center, scale, angleX, angleY, angleZ, toTriangulate);

    std::vector<Vector3D> points;
    std::vector<Face> faces;

    for (Face face: buckyBall.faces) {
        // Generate the points
        for (int i = 0; i < 3; i++) {
            points.push_back(buckyBall.points[face.point_indexes[i % 3]] +
                             (((buckyBall.points[face.point_indexes[(i + 1) % 3]]) -
                               (buckyBall.points[face.point_indexes[i % 3]])) * (1.0 / 3)));

            points.push_back((buckyBall.points[face.point_indexes[i % 3]]) +
                             (((buckyBall.points[face.point_indexes[(i + 1) % 3]]) -
                               (buckyBall.points[face.point_indexes[i % 3]])) * (2.0 / 3)));
        }
        // Generate the hexagon faces
        std::vector<int> hexagon;
        for (int i = 6; i > 0; i--) {
            hexagon.push_back(points.size() - i);
        }
        faces.emplace_back(hexagon);
    }
    // Add pentagon faces
    faces.emplace_back(Face({0, 5, 11, 17, 23}));
    faces.emplace_back(Face({27, 81, 30, 2, 1}));
    faces.emplace_back(Face({8, 4, 3, 33, 41}));
    faces.emplace_back(Face({9, 45, 53, 14, 10}));
    faces.emplace_back(Face({16, 15, 57, 65, 20}));
    faces.emplace_back(Face({77, 26, 22, 21, 69}));
    faces.emplace_back(Face({32, 31, 116, 94, 38}));
    faces.emplace_back(Face({40, 39, 91, 50, 44}));
    faces.emplace_back(Face({52, 51, 106, 62, 56}));
    faces.emplace_back(Face({68, 64, 63, 103, 74}));
    faces.emplace_back(Face({109, 86, 80, 76, 75}));
    faces.emplace_back(Face({95, 108, 102, 107, 90}));

    buckyBall.points = points;
    buckyBall.faces = faces;

    return buckyBall;
}
