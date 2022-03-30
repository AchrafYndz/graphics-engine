#ifndef ENGINE_LINE2D_H
#define ENGINE_LINE2D_H

#include "Point2D.h"
#include "Color.h"
#include "Easy_image/easy_image.h"

class Line2D {
public:
    Point2D p1;
    Point2D p2;
    Color color;
    double z1;
    double z2;
    Line2D(Point2D p1_, Point2D p2_, Color color_): p1(p1_), p2(p2_), color(color_) {};
    Line2D(Point2D p1_, Point2D p2_, Color color_, double z1_, double z2_): p1(p1_), p2(p2_), color(color_), z1(z1_), z2(z2_) {};
    img::Color getEzColor();
};


#endif //ENGINE_LINE2D_H
