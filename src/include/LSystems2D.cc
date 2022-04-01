#define _USE_MATH_DEFINES // Windows only

#include "LSystems2D.h"

#include <algorithm>
#include <cmath>

void draw2DLSystemHelper(const LParser::LSystem2D &l_system, Lines2D &lines, const Color col, int &recursionDepth,
                         const unsigned int maxRecursion, std::string currentString, double &currentAngle,
                         const double angleIncrement, double &x0, double &y0, std::stack<Brackets2D> &bracketStack) {
    if (recursionDepth == maxRecursion) {
        // Make the lines
        double x1;
        double y1;
        for (char c: currentString) {
            if (c == '+') currentAngle += angleIncrement;
            else if (c == '-') currentAngle -= angleIncrement;
            else if (c == '(') bracketStack.push(Brackets2D(x0, y0, currentAngle));
            else if (c == ')') {
                Brackets2D brackets = bracketStack.top();
                x0 = brackets.x;
                y0 = brackets.y;
                currentAngle = brackets.angle;
                bracketStack.pop();
            } else if (l_system.draw(c)) {
                x1 = x0 + cos(currentAngle);
                y1 = y0 + sin(currentAngle);
                lines.push_back(Line2D(Point2D(x0, y0), Point2D(x1, y1), col));
                x0 = x1;
                y0 = y1;
            } else if (find(l_system.get_alphabet().begin(), l_system.get_alphabet().end(), c) !=
                       l_system.get_alphabet().end()) {
                x1 = x0 + cos(currentAngle);
                y1 = y0 + sin(currentAngle);
                x0 = x1;
                y0 = y1;
            }
        }
        recursionDepth--;
    } else {
        for (char c: currentString) {
            if (c == '+') currentAngle += angleIncrement;
            else if (c == '-') currentAngle -= angleIncrement;
            else if (c == '(') bracketStack.push(Brackets2D(x0, y0, currentAngle));
            else if (c == ')') {
                Brackets2D brackets = bracketStack.top();
                x0 = brackets.x;
                y0 = brackets.y;
                currentAngle = brackets.angle;
                bracketStack.pop();
            } else if (find(l_system.get_alphabet().begin(), l_system.get_alphabet().end(), c) !=
                       l_system.get_alphabet().end()) {
                recursionDepth++;
                draw2DLSystemHelper(l_system, lines, col, recursionDepth, maxRecursion, l_system.get_replacement(c),
                                    currentAngle, angleIncrement, x0, y0, bracketStack);
            }
        }
        recursionDepth--;
    }
}

Lines2D draw2DLSystem(const LParser::LSystem2D &l_system, Color col) {
    Lines2D lines;
    // Call recursive function
    std::stack<Brackets2D> bracketStack;
    unsigned int Iterations = l_system.get_nr_iterations();
    std::string const &Initiator = l_system.get_initiator();
    double startingAngle = l_system.get_starting_angle() / 180 * M_PI;
    double currentAngle = startingAngle;
    double angleIncrement = l_system.get_angle() / 180 * M_PI;
    double x0 = 0;
    double y0 = 0;
    int recursionDepth = 0;
    draw2DLSystemHelper(l_system, lines, col, recursionDepth, Iterations, Initiator, currentAngle, angleIncrement,
                        x0, y0, bracketStack);
    return lines;
}