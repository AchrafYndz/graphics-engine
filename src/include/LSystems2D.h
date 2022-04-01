#ifndef ENGINE_LSYSTEMS2D_H
#define ENGINE_LSYSTEMS2D_H

#include <list>
#include <stack>

#include "../../util/Parser/l_parser.h"
#include "../lib/Line2D.h"

using Lines2D = std::list<Line2D>;

struct Brackets2D {
    double x;
    double y;
    double angle;

    Brackets2D(double x_, double y_, double angle_) : x(x_), y(y_), angle(angle_) {};
};

void draw2DLSystemHelper(const LParser::LSystem2D &l_system, Lines2D &lines, const Color col, int &recursionDepth,
                         const unsigned int maxRecursion, std::string currentString, double &currentAngle,
                         const double angleIncrement, double &x0, double &y0, std::stack<Brackets2D> &bracketStack);

Lines2D draw2DLSystem(const LParser::LSystem2D &l_system, Color col);


#endif //ENGINE_LSYSTEMS2D_H
