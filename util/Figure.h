#ifndef ENGINE_FIGURE_H
#define ENGINE_FIGURE_H

#include <vector>
#include <list>

#include "Color.h"
#include "Vector3D/vector3d.h"

class Face;

class Figure {
public:
    std::vector<Vector3D> points;
    std::vector<Face> faces;
    Color color;
    // Transformation matrices
    Matrix scaleFigure(const double scale);
    Matrix rotateX(const double angle);
    Matrix rotateY(const double angle);
    Matrix rotateZ(const double angle);
    Matrix translate(const Vector3D& vector);
};

typedef std::list<Figure> Figures3D;

#endif //ENGINE_FIGURE_H
