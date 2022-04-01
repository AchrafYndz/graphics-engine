#ifndef ENGINE_LSYSTEMS3D_H
#define ENGINE_LSYSTEMS3D_H

#include "../../util/Vector3D/vector3d.h"
#include "../../util/Parser/l_parser.h"
#include "../lib/Color.h"
#include "../lib/Figure.h"

#include <vector>
#include <stack>

struct Brackets3D {
    Vector3D position;

    Vector3D H;
    Vector3D L;
    Vector3D U;

    Brackets3D(const Vector3D &position_, const Vector3D &H_, const Vector3D &L_, const Vector3D &U_) : position(
            position_), H(H_), L(L_), U(U_) {};
};

void
draw3DLSystemHelper(const LParser::LSystem3D &l_system, std::vector<Vector3D> &points, std::vector<Face> &faces,
                    const Color col, int &recursionDepth,
                    const unsigned int maxRecursion, const std::string &currentString, double &angle,
                    Vector3D &position,
                    Vector3D &H, Vector3D &L, Vector3D &U, std::stack<Brackets3D> &bracketStack);

Figure draw3DLSystem(const LParser::LSystem3D &l_system, Vector3D &center, Color color, const double &scale,
                     const double &angleX, const double &angleY, const double &angleZ);


#endif //ENGINE_LSYSTEMS3D_H
