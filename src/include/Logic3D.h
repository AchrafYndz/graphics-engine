#ifndef ENGINE_LOGIC3D_H
#define ENGINE_LOGIC3D_H


#include "../../util/Easy_image/easy_image.h"
#include "LSystems2D.h"
#include "../../util/Vector3D/vector3d.h"
#include "../lib/Figure.h"

img::EasyImage draw2DLines(const Lines2D &lines, int size, img::Color &bg_col, bool zBuffer);

Matrix scaleFigure(double scale);

Matrix rotateX(double angle);

Matrix rotateY(double angle);

Matrix rotateZ(double angle);

Matrix translate(const Vector3D &vector);

void toPolar(const Vector3D &point, double &theta, double &phi, double &r);

Matrix eyePointTrans(const Vector3D &eyepoint);

Point2D doProjection(const Vector3D &point, const int d);

void applyTransformation(Figure &fig, const Matrix &m);

Lines2D doProjection(Figures3D &figs, const Vector3D &eyepoint);



#endif //ENGINE_LOGIC3D_H
