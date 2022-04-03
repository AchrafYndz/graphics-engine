#ifndef ENGINE_ZBUFFERTRIANGLES_H
#define ENGINE_ZBUFFERTRIANGLES_H

#include <vector>
#include "../lib/Face.h"
#include "ZBufferWireframes.h"
#include "../../util/Vector3D/vector3d.h"
#include "../lib/Color.h"
#include "LSystems2D.h"


std::vector<Face> triangulate(const Face &face);

void
draw_zbuf_trag(ZBuffer &zbuffer, img::EasyImage &image, Vector3D const &A, Vector3D const &B, Vector3D const &C, double d,
               double dx, double dy, img::Color color);

void getImageSpecs(const Lines2D &lines, int size, double &d, double &dx, double &dy, double &width, double &height);

Point2D doProjection(const Vector3D point, const double d, const double dx, const double dy);

#endif //ENGINE_ZBUFFERTRIANGLES_H
