#ifndef ENGINE_FIGURES3D_H
#define ENGINE_FIGURES3D_H


#include "../lib/Figure.h"

Figure createCube(Color ambientReflection, Vector3D &center, double scale, double angleX, double angleY, double angleZ, bool toTriangulate);

Figure createTetrahedron(Color ambientReflection, Vector3D &center, double scale, double angleX, double angleY, double angleZ, bool toTriangulate);

Figure createOctahedron(Color ambientReflection, Vector3D &center, double scale, double angleX, double angleY, double angleZ, bool toTriangulate);

Figure createIcosahedron(Color ambientReflection, Vector3D &center, double scale, double angleX, double angleY, double angleZ, bool toTriangulate);

Figure createDodecahedron(Color ambientReflection, Vector3D &center, double scale, double angleX, double angleY, double angleZ, bool toTriangulate);

Figure
createSphere(Color ambientReflection, Vector3D &center, double scale, double angleX, double angleY, double angleZ, int n, bool toTriangulate);

Figure createCone(Color ambientReflection, Vector3D &center, double scale, double angleX, double angleY, double angleZ, int n,
                  double h, bool toTriangulate);

Figure
createCylinder(Color ambientReflection, Vector3D &center, double scale, double angleX, double angleY, double angleZ, int n,
               double h, bool toTriangulate);

Figure
createTorus(Color ambientReflection, Vector3D &center, double scale, double angleX, double angleY, double angleZ, double r,
            double R, int n, int m, bool toTriangulate);

#endif //ENGINE_FIGURES3D_H
