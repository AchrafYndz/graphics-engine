#include "Fractals3D.h"
#include "Logic3D.h"
#include "Figures3D.h"
#include "ZBufferTriangles.h"

void generateFractal(Figure &fig, Figures3D &fractal, const int nrIterations, const double scale) {
    fractal.push_back(fig);
    for (int i = 0; i < nrIterations; i++) {
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

void createMengerSponge(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ,
                        bool toTriangulate, int nrIterations, Figures3D &cubes) {
    // Start with a cube
    Figure mengerSponge = createCube(color, center, scale, angleX, angleY, angleZ, false);

    // Push the initial cube onto the list of figures
    cubes.push_back(mengerSponge);

    for (int i = 0; i < nrIterations; i++) {
        Figures3D newSponges;

        Matrix scaleMatrix = scaleFigure(1.0 / 3);
        for (Figure cube: cubes) {
            // put corner cubes
            for (int j = 0; j < 8; j++) {
                Figure copySponge = cube;
                applyTransformation(copySponge, scaleMatrix);
                Matrix translationMatrix = translate(cube.points[j] - copySponge.points[j]);
                applyTransformation(copySponge, translationMatrix);
                newSponges.push_back(copySponge);
            }
            // put side cubes
            for (Face square: cube.faces) {
                for (int p = 0; p < 3; p++) {
                    Figure copySponge = cube;
                    applyTransformation(copySponge, scaleMatrix);

                    Vector3D midPointOuter =
                            (cube.points[square.point_indexes[p]] + cube.points[square.point_indexes[p + 1]]) / 2;

                    Vector3D midPointInner =
                            (copySponge.points[square.point_indexes[p]] +
                             copySponge.points[square.point_indexes[p + 1]]) / 2;

                    Matrix translationMatrix = translate(midPointOuter - midPointInner);
                    applyTransformation(copySponge, translationMatrix);
                    newSponges.push_back(copySponge);
                }
            }
        }
        cubes = newSponges;
    }

    // Triangulation
    if (toTriangulate) {
        for (Figure& cube: cubes) {
            std::vector<Face> triangles;
            for (const Face &face: cube.faces) {
                std::vector<Face> newTriangles = triangulate(face);
                triangles.insert(triangles.end(), newTriangles.begin(), newTriangles.end());
            }
            cube.faces = triangles;
        }
    }
}