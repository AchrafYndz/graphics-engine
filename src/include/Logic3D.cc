#define _USE_MATH_DEFINES // Windows only

#include "Logic3D.h"
#include "ZBufferWireframes.h"

#include <cmath>
#include <fstream>

img::EasyImage draw2DLines(const Lines2D &lines, const int size, img::Color &bg_col, const bool zBuffer) {
    double xmin = lines.front().p1.x;
    double xmax = lines.front().p1.y;
    double ymin = lines.front().p2.x;
    double ymax = lines.front().p2.y;
    // determine max and min
    for (Line2D line: lines) {
        if (line.p1.x < xmin) xmin = line.p1.x;
        if (line.p2.x < xmin) xmin = line.p2.x;
        if (line.p1.y < ymin) ymin = line.p1.y;
        if (line.p2.y < ymin) ymin = line.p2.y;

        if (line.p1.x > xmax) xmax = line.p1.x;
        if (line.p2.x > xmax) xmax = line.p2.x;
        if (line.p1.y > ymax) ymax = line.p1.y;
        if (line.p2.y > ymax) ymax = line.p2.y;
    }

    double xrange = std::abs(xmax - xmin);
    double yrange = std::abs(ymax - ymin);

    double imagex = size * xrange / std::max(xrange, yrange);
    double imagey = size * yrange / std::max(xrange, yrange);

    img::EasyImage image(lround(imagex), lround(imagey), bg_col);

    double d = 0.95 * imagex / xrange;

    double DCx = d * (xmin + xmax) / 2;
    double DCy = d * (ymin + ymax) / 2;

    double dx = imagex / 2 - DCx;
    double dy = imagey / 2 - DCy;

    ZBuffer zbuf(image.get_width(), image.get_height());
    for (Line2D line: lines) {
        if (zBuffer)
            draw_zbuf_line(zbuf, image, lround(line.p1.x * d + dx), lround(line.p1.y * d + dy), line.z1,
                           lround(line.p2.x * d + dx),
                           lround(line.p2.y * d + dy), line.z2, line.getEzColor());
        else
            image.draw_line(lround(line.p1.x * d + dx), lround(line.p1.y * d + dy),
                            lround(line.p2.x * d + dx),
                            lround(line.p2.y * d + dy), line.getEzColor());
    }
    std::ofstream fout("out.bmp", std::ios::binary);
    fout << image;
    fout.close();
    return image;
}

Matrix scaleFigure(const double scale) {
    Matrix scaleMatrix;
    for (int i = 1; i < 4; i++) {
        scaleMatrix(i, i) = scale;
    }
    return scaleMatrix;
}

Matrix rotateX(const double angle) {
    Matrix rotationMatrix;
    rotationMatrix(2, 2) = cos(angle);
    rotationMatrix(2, 3) = sin(angle);
    rotationMatrix(3, 2) = -sin(angle);
    rotationMatrix(3, 3) = cos(angle);
    return rotationMatrix;
}

Matrix rotateY(const double angle) {
    Matrix rotationMatrix;
    rotationMatrix(1, 1) = cos(angle);
    rotationMatrix(1, 3) = -sin(angle);
    rotationMatrix(3, 1) = sin(angle);
    rotationMatrix(3, 3) = cos(angle);
    return rotationMatrix;
}

Matrix rotateZ(const double angle) {
    Matrix rotationMatrix;
    rotationMatrix(1, 1) = cos(angle);
    rotationMatrix(1, 2) = sin(angle);
    rotationMatrix(2, 1) = -sin(angle);
    rotationMatrix(2, 2) = cos(angle);
    return rotationMatrix;
}

Matrix translate(const Vector3D &vector) {
    Matrix translationMatrix;
    translationMatrix(4, 1) = vector.x;
    translationMatrix(4, 2) = vector.y;
    translationMatrix(4, 3) = vector.z;
    return translationMatrix;
}

void toPolar(const Vector3D &point, double &theta, double &phi, double &r) {
    r = sqrt(pow(point.x, 2) + pow(point.y, 2) + pow(point.z, 2));
    theta = atan2(point.y, point.x);
    phi = acos(point.z / r);
}

Matrix eyePointTrans(const Vector3D &eyepoint) {
    double theta;
    double phi;
    double r;
    toPolar(eyepoint, theta, phi, r);
    Matrix eyePointMatrix;
    eyePointMatrix(1, 1) = -sin(theta);
    eyePointMatrix(1, 2) = -cos(theta) * cos(phi);
    eyePointMatrix(1, 3) = cos(theta) * sin(phi);
    eyePointMatrix(2, 1) = cos(theta);
    eyePointMatrix(2, 2) = -sin(theta) * cos(phi);
    eyePointMatrix(2, 3) = sin(theta) * sin(phi);
    eyePointMatrix(3, 2) = sin(phi);
    eyePointMatrix(3, 3) = cos(phi);
    eyePointMatrix(4, 3) = -r;
    return eyePointMatrix;
}

Point2D doProjection(const Vector3D &point, const double d) {
    return {-point.x * d / point.z, -point.y * d / point.z};
}

void applyTransformation(Figure &fig, const Matrix &m) {
    for (Vector3D &point: fig.points) {
        point *= m;
    }
}

Lines2D doProjection(Figures3D &figs, const Vector3D &eyepoint) {
    Lines2D projection;
    for (Figure &fig: figs) {
        Matrix m = scaleFigure(fig.scale) * rotateX(fig.rotateAngleX) * rotateY(fig.rotateAngleY) *
                   rotateZ(fig.rotateAngleZ) *
                   translate(fig.center) * eyePointTrans(eyepoint);
        applyTransformation(fig, m);
        for (Face &face: fig.faces) {
            for (int i = 0; i < face.point_indexes.size(); i++) {
                if (i != face.point_indexes.size() - 1) {
                    Vector3D p0 = fig.points[face.point_indexes[i]];
                    Vector3D p1 = fig.points[face.point_indexes[i + 1]];
                    Point2D x = doProjection(p0, 1);
                    Point2D y = doProjection(p1, 1);
                    projection.push_back(Line2D(x, y, fig.color, p0.z, p1.z));
                } else {
                    Vector3D p0 = fig.points[face.point_indexes[i]];
                    Vector3D p1 = fig.points[face.point_indexes[0]];
                    Point2D x = doProjection(p0, 1);
                    Point2D y = doProjection(p1, 1);
                    projection.push_back(Line2D(x, y, fig.color, p0.z, p1.z));
                }
            }
        }
    }
    return projection;
}