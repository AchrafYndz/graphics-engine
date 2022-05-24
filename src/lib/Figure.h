#ifndef ENGINE_FIGURE_H
#define ENGINE_FIGURE_H

#include <vector>
#include <list>

#include "Color.h"
#include "../../util/Vector3D/vector3d.h"
#include "Face.h"
#include "../../util/Easy_image/easy_image.h"

class Figure {
public:
    std::vector<Vector3D> points;
    std::vector<Face> faces;
    Vector3D center;

    Color ambientReflection;
    Color diffuseReflection;
    Color specularReflection;

    double reflectionCoefficient;

    double rotateAngleX;
    double rotateAngleY;
    double rotateAngleZ;
    double scale;

    Figure(std::vector<Vector3D> &points_, std::vector<Face> &faces_, Color ambientReflection_, Vector3D &center_,
           double scale_,
           double angleX_, double angleY_, double angleZ_) : points(points_), faces(faces_),
                                                             ambientReflection(ambientReflection_),
                                                             center(center_), scale(scale_), rotateAngleX(angleX_),
                                                             rotateAngleY(angleY_), rotateAngleZ(angleZ_) {};

    img::Color getEzColor();
};

typedef std::list<Figure> Figures3D;

#endif //ENGINE_FIGURE_H
