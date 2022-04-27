#ifndef ENGINE_FRACTALS3D_H
#define ENGINE_FRACTALS3D_H

#include "../lib/Figure.h"

void generateFractal(Figure& fig, Figures3D& fractal, int nr_iterations, double scale);

Figure createBuckyBall(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ,
                        bool toTriangulate);

void createMengerSponge(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ,
                          bool toTriangulate, int nrIterations, Figures3D& sponges);

#endif //ENGINE_FRACTALS3D_H