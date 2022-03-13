#ifndef ENGINE_FIGURE_H
#define ENGINE_FIGURE_H

#include <vector>
#include <list>

#include "Color.h"
#include "Vector3D/vector3d.h"
#include "Face.h"

class Figure {
public:
    std::vector<Vector3D> points;
    std::vector<Face> faces;
    Vector3D center;
    Color color;
    double rotateAngleX;
    double rotateAngleY;
    double rotateAngleZ;
    Figure(std::vector<Vector3D> points_, std::vector<Face> faces_, Color color_, Vector3D center_, double angleX_, double angleY_, double angleZ_): points(points_), faces(faces_), color(color_), center(center_), rotateAngleX(angleX_), rotateAngleY(angleY_), rotateAngleZ(angleY_){};
};

typedef std::list<Figure> Figures3D;

#endif //ENGINE_FIGURE_H
