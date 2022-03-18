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
    double scale;
    Figure(std::vector<Vector3D>& points_, std::vector<Face>& faces_, Color color_, Vector3D& center_, double scale_, double angleX_, double angleY_, double angleZ_): points(points_), faces(faces_), color(color_), center(center_), scale(scale_), rotateAngleX(angleX_), rotateAngleY(angleY_), rotateAngleZ(angleZ_){};

};

typedef std::list<Figure> Figures3D;

#endif //ENGINE_FIGURE_H
