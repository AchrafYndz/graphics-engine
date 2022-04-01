#define _USE_MATH_DEFINES // Windows only


#include "LSystems3D.h"

#include <cmath>
#include <algorithm>

void draw3DLSystemHelper(const LParser::LSystem3D &l_system, std::vector<Vector3D> &points,
                                     std::vector<Face> &faces, const Color col, int &recursionDepth,
                                     const unsigned int maxRecursion, const std::string &currentString, double &angle,
                                     Vector3D &position, Vector3D &H, Vector3D &L, Vector3D &U,
                                     std::stack<Brackets3D> &bracketStack) {
    if (recursionDepth == maxRecursion) {
        // Make the lines
        for (char c: currentString) {
            if (c == '+') {
                Vector3D prevH = H;
                H = H * cos(angle) + L * sin(angle);
                L = -prevH * sin(angle) + L * cos(angle);
            } else if (c == '-') {
                Vector3D prevH = H;
                H = H * cos(-angle) + L * sin(-angle);
                L = -prevH * sin(-angle) + L * cos(-angle);
            } else if (c == '^') {
                Vector3D prevH = H;
                H = H * cos(angle) + U * sin(angle);
                U = -prevH * sin(angle) + U * cos(angle);
            } else if (c == '&') {
                Vector3D prevH = H;
                H = H * cos(-angle) + U * sin(-angle);
                U = -prevH * sin(-angle) + U * cos(-angle);
            } else if (c == '\\') {
                Vector3D prevL = L;
                L = L * cos(angle) - U * sin(angle);
                U = prevL * sin(angle) + U * cos(angle);
            } else if (c == '/') {
                Vector3D prevL = L;
                L = L * cos(-angle) - U * sin(-angle);
                U = prevL * sin(-angle) + U * cos(-angle);
            } else if (c == '|') {
                H = -H;
                L = -L;
            } else if (c == '(') bracketStack.push(Brackets3D(position, H, L, U));
            else if (c == ')') {
                Brackets3D brackets = bracketStack.top();
                position = brackets.position;
                H = brackets.H;
                L = brackets.L;
                U = brackets.U;
                bracketStack.pop();
            } else if (l_system.draw(c)) {
                points.push_back(position);
                position += H;
                points.push_back(position);
                faces.emplace_back(
                        std::vector<int>{static_cast<int>(points.size()) - 2, static_cast<int>(points.size()) - 1});
            } else if (find(l_system.get_alphabet().begin(), l_system.get_alphabet().end(), c) !=
                       l_system.get_alphabet().end()) {
                position += H;
            }
        }
        recursionDepth--;
    } else {
        for (char c: currentString) {
            if (c == '+') {
                Vector3D prevH = H;
                H = H * cos(angle) + L * sin(angle);
                L = -prevH * sin(angle) + L * cos(angle);
            } else if (c == '-') {
                Vector3D prevH = H;
                H = H * cos(-angle) + L * sin(-angle);
                L = -prevH * sin(-angle) + L * cos(-angle);
            } else if (c == '^') {
                Vector3D prevH = H;
                H = H * cos(angle) + U * sin(angle);
                U = -prevH * sin(angle) + U * cos(angle);
            } else if (c == '&') {
                Vector3D prevH = H;
                H = H * cos(-angle) + U * sin(-angle);
                U = -prevH * sin(-angle) + U * cos(-angle);
            } else if (c == '\\') {
                Vector3D prevL = L;
                L = L * cos(angle) - U * sin(angle);
                U = prevL * sin(angle) + U * cos(angle);
            } else if (c == '/') {
                Vector3D prevL = L;
                L = L * cos(-angle) - U * sin(-angle);
                U = prevL * sin(-angle) + U * cos(-angle);
            } else if (c == '|') {
                H = -H;
                L = -L;
            } else if (c == '(') bracketStack.push(Brackets3D(position, H, L, U));
            else if (c == ')') {
                Brackets3D brackets = bracketStack.top();
                position = brackets.position;
                H = brackets.H;
                L = brackets.L;
                U = brackets.U;
                bracketStack.pop();
            } else if (find(l_system.get_alphabet().begin(), l_system.get_alphabet().end(), c) !=
                       l_system.get_alphabet().end()) {
                recursionDepth++;
                draw3DLSystemHelper(l_system, points, faces, col, recursionDepth, maxRecursion,
                                    l_system.get_replacement(c),
                                    angle, position, H, L, U, bracketStack);
            }
        }
        recursionDepth--;
    }
}

Figure draw3DLSystem(const LParser::LSystem3D &l_system, Vector3D &center, Color color, const double &scale,
                     const double &angleX, const double &angleY, const double &angleZ) {
    // Call recursive function
    std::stack<Brackets3D> bracketStack;
    Vector3D position = Vector3D::point(0, 0, 0);
    std::vector<Vector3D> points;
    std::vector<Face> faces;
    Vector3D H = Vector3D::vector(1, 0, 0);
    Vector3D L = Vector3D::vector(0, 1, 0);
    Vector3D U = Vector3D::vector(0, 0, 1);
    unsigned int Iterations = l_system.get_nr_iterations();
    std::string const &Initiator = l_system.get_initiator();
    double angle = l_system.get_angle();
    angle = angle / 180 * M_PI;
    int recursionDepth = 0;
    draw3DLSystemHelper(l_system, points, faces, color, recursionDepth, Iterations, Initiator, angle,
                        position, H, L, U, bracketStack);
    return {points, faces, color, center, scale, angleX, angleY, angleZ};
}

