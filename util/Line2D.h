#ifndef ENGINE_LINE2D_H
#define ENGINE_LINE2D_H

#include "Point2D.h"
#include "Color.h"
#include "../easy_image.h"

class Line2D {
public:
    Point2D p1;
    Point2D p2;
    Color color;
    Line2D(Point2D point1, Point2D point2, Color col): p1(point1), p2(point2), color(col) {};
    img::Color getEzColor();
};


#endif //ENGINE_LINE2D_H
