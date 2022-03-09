#include "Figure.h"
#include <cmath>

Matrix Figure::scaleFigure(const double scale) {
    Matrix scaleMatrix;
    for (int i=0; i<3; i++) {
        scaleMatrix(i, i) = scale;
    }
    return scaleMatrix;
}

Matrix Figure::rotateX(const double angle) {
    Matrix rotationMatrix;
    rotationMatrix(1, 1) = cos(angle);
    rotationMatrix(1, 2) = sin(angle);
    rotationMatrix(2, 1) = -sin(angle);
    rotationMatrix(1, 2) = cos(angle);
    return rotationMatrix;
}

Matrix Figure::rotateY(const double angle) {
    Matrix rotationMatrix;
    rotationMatrix(0, 0) = cos(angle);
    rotationMatrix(0, 2) = -sin(angle);
    rotationMatrix(2, 0) = sin(angle);
    rotationMatrix(2, 2) = cos(angle);
    return rotationMatrix;
}

Matrix Figure::rotateZ(const double angle) {
    Matrix rotationMatrix;
    rotationMatrix(0, 0) = cos(angle);
    rotationMatrix(0, 1) = sin(angle);
    rotationMatrix(1, 0) = -sin(angle);
    rotationMatrix(1, 1) = cos(angle);
    return rotationMatrix;
}

Matrix Figure::translate(const Vector3D &vector) {
    Matrix translationMatrix;
    translationMatrix(3, 0) = vector.x;
    translationMatrix(3, 1) = vector.y;
    translationMatrix(3, 2) = vector.z;
    return translationMatrix;
}
