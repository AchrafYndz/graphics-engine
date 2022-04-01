#ifndef ENGINE_FIGURES3D_H
#define ENGINE_FIGURES3D_H


#include "../lib/Figure.h"

Figure createCube(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ);

Figure createTetrahedron(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ);

Figure createOctahedron(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ);

Figure createIcosahedron(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ);

Figure createDodecahedron(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ);

Figure
createSphere(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, int n);

Figure createCone(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, int n,
                  double h);

Figure
createCylinder(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, int n,
               double h);

Figure
createTorus(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, double r,
            double R, int n, int m);

#endif //ENGINE_FIGURES3D_H
