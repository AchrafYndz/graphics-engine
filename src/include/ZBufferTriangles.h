#ifndef ENGINE_ZBUFFERTRIANGLES_H
#define ENGINE_ZBUFFERTRIANGLES_H

#include <vector>
#include "../lib/Face.h"
#include "ZBufferWireframes.h"
#include "../../util/Vector3D/vector3d.h"
#include "../lib/Color.h"

std::vector<Face> triangulate(const Face &face);

void
draw_zbuf_trag(ZBuffer &zbuffer, img::EasyImage &image, Vector3D const &A, Vector3D const &B, Vector3D const &C, double d,
               double dx, double dy, img::Color color);



#endif //ENGINE_ZBUFFERTRIANGLES_H
