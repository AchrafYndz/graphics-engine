#ifndef ENGINE_FIGURES3D_H
#define ENGINE_FIGURES3D_H


#include "../lib/Figure.h"

Figure createCube(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, bool toTriangulate);

Figure createTetrahedron(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, bool toTriangulate);

Figure createOctahedron(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, bool toTriangulate);

Figure createIcosahedron(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, bool toTriangulate);

Figure createDodecahedron(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, bool toTriangulate);

Figure
createSphere(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, int n, bool toTriangulate);

Figure createCone(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, int n,
                  double h, bool toTriangulate);

Figure
createCylinder(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, int n,
               double h, bool toTriangulate);

Figure
createTorus(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, double r,
            double R, int n, int m, bool toTriangulate);

#endif //ENGINE_FIGURES3D_H
