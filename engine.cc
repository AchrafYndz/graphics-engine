#define _USE_MATH_DEFINES

#include "util/Easy_image/easy_image.h"
#include "util/Ini_config/ini_configuration.h"
#include "util/Line2D.h"
#include "util/Parser/l_parser.h"
#include "util/Vector3D/vector3d.h"
#include "util/Figure.h"
#include "util/Face.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cmath>
#include <list>
#include <stack>


using Lines2D = std::list<Line2D>;

struct Brackets {
    double x;
    double y;
    double angle;

    Brackets(double x_, double y_, double angle_) : x(x_), y(y_), angle(angle_) {};
};

void ColorRectangle(img::EasyImage &img) {
    for (unsigned int i = 0; i < 256; i++) {
        for (unsigned int j = 0; j < 256; j++) {
            img(i, j).red = i;
            img(i, j).green = j;
            img(i, j).blue = (i + j) % 256;
        }
    }
}

void
Blocks(img::EasyImage &img, int Wi, int Hi, int Nx, int Ny, std::vector<double> color1, std::vector<double> color2) {
    int Wb = Wi / Nx;
    int Hb = Hi / Ny;
    for (int px = 0; px < Wi; px++) {
        for (int py = 0; py < Hi; py++) {
            int Bx = px / Wb;
            int By = py / Hb;
            if ((Bx + By) % 2 == 0) {
                img(px, py).red = color1[0] * 255;
                img(px, py).green = color1[1] * 255;
                img(px, py).blue = color1[2] * 255;
            } else {
                img(px, py).red = color2[0] * 255;
                img(px, py).green = color2[1] * 255;
                img(px, py).blue = color2[2] * 255;
            }
        }
    }
}

void QuarterCircle(img::EasyImage &img, int Hi, int Wi, int N, std::vector<double> lineColor,
                   std::vector<double> backgroundColor, std::string figure = "") {
    int Hs = Hi / (N - 1);
    int Ws = Wi / (N - 1);
    img::Color li_col;
    li_col.red = lineColor[0] * 255;
    li_col.green = lineColor[1] * 255;
    li_col.blue = lineColor[2] * 255;

    if (figure == "eye") {
        for (int n = 0; n < N - 1; n++) {
            img.draw_line(Ws * n, Hi - 1, 0, Hs * n, li_col);
        }
        for (int n = 0; n < N - 1; n++) {
            img.draw_line(Hi - 1, Ws * n, Hs * n, 0, li_col);
        }
    } else if (figure == "diamond") {
        Wi = Wi / 2;
        Hi = Hi / 2;
        Ws = Wi / (N - 1);
        Hs = Hi / (N - 1);
        for (int n = 0; n < N - 1; n++) {
            img.draw_line(Ws * n + Wi, Hi - 1, Wi, Hs * n, li_col);
        }
        for (int n = 0; n < N - 1; n++) {
            img.draw_line(Hi - 1, Ws * n + Hi, Hs * n, Hi, li_col);
        }
        for (int n = 0; n < N - 1; n++) {
            img.draw_line(Hi - 1, Ws * n, Hs * (N - n - 2), Hi, li_col);
        }
        for (int n = 0; n < N - 1; n++) {
            img.draw_line(Hi - 1, Ws * n + Wi, Hs * (N - n - 2) + Hi, Hi - 1, li_col);
        }
    } else {
        for (int n = 0; n < N - 1; n++) {
            img.draw_line(Ws * n, Hi - 1, 0, Hs * n, li_col);
        }
    }
}


void
Eye(img::EasyImage &img, int Hi, int Wi, int N, std::vector<double> lineColor, std::vector<double> backgroundColor) {
    QuarterCircle(img, Hi, Wi, N, lineColor, backgroundColor, "eye");
}

void Diamond(img::EasyImage &img, int Hi, int Wi, int N, std::vector<double> lineColor,
             std::vector<double> backgroundColor) {
    QuarterCircle(img, Hi, Wi, N, lineColor, backgroundColor, "diamond");
}

img::EasyImage draw2DLines(const Lines2D &lines, const int size, img::Color bg_col) {
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

    for (Line2D line: lines) {
        image.draw_line(lround(line.p1.x * d + dx), lround(line.p1.y * d + dy), lround(line.p2.x * d + dx),
                        lround(line.p2.y * d + dy), line.getEzColor());
    }
    std::ofstream fout("out.bmp", std::ios::binary);
    fout << image;
    fout.close();
    return image;
}

void drawLSystemHelper(const LParser::LSystem2D &l_system, Lines2D &lines, const Color col, int &recursionDepth,
                       const unsigned int maxRecursion, std::string currentString, double &currentAngle,
                       const double angleIncrement, double &x0, double &y0, std::stack<Brackets> &bracketStack) {
    if (recursionDepth == maxRecursion) {
        // Make the lines
        double x1;
        double y1;
        for (char c: currentString) {
            if (c == '+') currentAngle += angleIncrement;
            else if (c == '-') currentAngle -= angleIncrement;
            else if (c == '(') bracketStack.push(Brackets(x0, y0, currentAngle));
            else if (c == ')') {
                Brackets brackets = bracketStack.top();
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
            }
        }
        recursionDepth--;
    } else {
        for (char c: currentString) {
            if (c == '+') currentAngle += angleIncrement;
            else if (c == '-') currentAngle -= angleIncrement;
            else if (c == '(') bracketStack.push(Brackets(x0, y0, currentAngle));
            else if (c == ')') {
                Brackets brackets = bracketStack.top();
                x0 = brackets.x;
                y0 = brackets.y;
                currentAngle = brackets.angle;
                bracketStack.pop();
            } else if (l_system.draw(c)) {
                recursionDepth++;
                drawLSystemHelper(l_system, lines, col, recursionDepth, maxRecursion, l_system.get_replacement(c),
                                  currentAngle, angleIncrement, x0, y0, bracketStack);
            }
        }
        recursionDepth--;
    }
}

Lines2D drawLSystem(const LParser::LSystem2D &l_system, Color col) {
    Lines2D lines;
    // Call recursive function
    std::stack<Brackets> bracketStack;
    unsigned int Iterations = l_system.get_nr_iterations();
    std::string const &Initiator = l_system.get_initiator();
    double startingAngle = l_system.get_starting_angle() / 180 * M_PI;
    double currentAngle = startingAngle;
    double angleIncrement = l_system.get_angle() / 180 * M_PI;
    double x0 = 0;
    double y0 = 0;
    int recursionDepth = 0;
    drawLSystemHelper(l_system, lines, col, recursionDepth, Iterations, Initiator, currentAngle, angleIncrement,
                      x0, y0, bracketStack);
    return lines;
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

Point2D doProjection(const Vector3D &point, const int d) {
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
                    projection.push_back(Line2D(x, y, fig.color));
                } else {
                    Vector3D p0 = fig.points[face.point_indexes[i]];
                    Vector3D p1 = fig.points[face.point_indexes[0]];
                    Point2D x = doProjection(p0, 1);
                    Point2D y = doProjection(p1, 1);
                    projection.push_back(Line2D(x, y, fig.color));
                }
            }
        }
    }
    return projection;
}

// 3D Figures
Figure createCube(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ) {
    // Create points
    std::vector<Vector3D> points;

    Vector3D p0 = Vector3D::point(1, -1, -1);
    points.push_back(p0);

    Vector3D p1 = Vector3D::point(-1, 1, -1);
    points.push_back(p1);

    Vector3D p2 = Vector3D::point(1, 1, 1);
    points.push_back(p2);

    Vector3D p3 = Vector3D::point(-1, -1, 1);
    points.push_back(p3);

    Vector3D p4 = Vector3D::point(1, 1, -1);
    points.push_back(p4);

    Vector3D p5 = Vector3D::point(-1, -1, -1);
    points.push_back(p5);

    Vector3D p6 = Vector3D::point(1, -1, 1);
    points.push_back(p6);

    Vector3D p7 = Vector3D::point(-1, 1, 1);
    points.push_back(p7);

    // Create faces
    std::vector<Face> faces;

    std::vector<int> point_indexes0 = {0, 4, 2, 6};
    Face face0(point_indexes0);
    faces.push_back(face0);

    std::vector<int> point_indexes1 = {4, 1, 7, 2};
    Face face1(point_indexes1);
    faces.push_back(face1);

    std::vector<int> point_indexes2 = {1, 5, 3, 7};
    Face face2(point_indexes2);
    faces.push_back(face2);

    std::vector<int> point_indexes3 = {5, 0, 6, 3};
    Face face3(point_indexes3);
    faces.push_back(face3);

    std::vector<int> point_indexes4 = {6, 2, 7, 3};
    Face face4(point_indexes4);
    faces.push_back(face4);

    std::vector<int> point_indexes5 = {0, 5, 1, 4};
    Face face5(point_indexes5);
    faces.push_back(face5);

    // Create figure
    Figure cube(points, faces, color, center, scale, angleX, angleY, angleZ);
    return cube;
}

Figure createTetrahedron(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ) {
    // Create points
    std::vector<Vector3D> points;
    Vector3D p0 = Vector3D::point(1, -1, -1);
    points.push_back(p0);
    Vector3D p1 = Vector3D::point(-1, 1, -1);
    points.push_back(p1);
    Vector3D p2 = Vector3D::point(1, 1, 1);
    points.push_back(p2);
    Vector3D p3 = Vector3D::point(-1, -1, 1);
    points.push_back(p3);

    // Create faces
    std::vector<Face> faces;
    std::vector<int> point_indexes0 = {0, 1, 2};
    Face face0(point_indexes0);
    faces.push_back(face0);
    std::vector<int> point_indexes1 = {1, 3, 2};
    Face face1(point_indexes1);
    faces.push_back(face1);
    std::vector<int> point_indexes2 = {0, 3, 1};
    Face face2(point_indexes2);
    faces.push_back(face2);
    std::vector<int> point_indexes3 = {0, 2, 3};
    Face face3(point_indexes3);
    faces.push_back(face3);

    // Create figure
    Figure tetrahedron(points, faces, color, center, scale, angleX, angleY, angleZ);
    return tetrahedron;
}

Figure createOctahedron(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ) {
    // Create points
    std::vector<Vector3D> points;

    Vector3D p0 = Vector3D::point(1, 0, 0);
    points.push_back(p0);

    Vector3D p1 = Vector3D::point(0, 1, 0);
    points.push_back(p1);

    Vector3D p2 = Vector3D::point(-1, 0, 0);
    points.push_back(p2);

    Vector3D p3 = Vector3D::point(0, -1, 0);
    points.push_back(p3);

    Vector3D p4 = Vector3D::point(0, 0, -1);
    points.push_back(p4);

    Vector3D p5 = Vector3D::point(0, 0, 1);
    points.push_back(p5);

    // Create faces
    std::vector<Face> faces;

    std::vector<int> point_indexes0 = {0, 1, 5};
    Face face0(point_indexes0);
    faces.push_back(face0);

    std::vector<int> point_indexes1 = {1, 2, 5};
    Face face1(point_indexes1);
    faces.push_back(face1);

    std::vector<int> point_indexes2 = {2, 3, 5};
    Face face2(point_indexes2);
    faces.push_back(face2);

    std::vector<int> point_indexes3 = {3, 0, 5};
    Face face3(point_indexes3);
    faces.push_back(face3);

    std::vector<int> point_indexes4 = {1, 0, 4};
    Face face4(point_indexes4);
    faces.push_back(face4);

    std::vector<int> point_indexes5 = {2, 1, 4};
    Face face5(point_indexes5);
    faces.push_back(face5);

    std::vector<int> point_indexes6 = {3, 2, 4};
    Face face6(point_indexes6);
    faces.push_back(face6);

    std::vector<int> point_indexes7 = {0, 3, 4};
    Face face7(point_indexes7);
    faces.push_back(face7);

    // Create figure
    Figure octahedron(points, faces, color, center, scale, angleX, angleY, angleZ);
    return octahedron;
}

Figure createIcosahedron(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ) {
    // Create points
    std::vector<Vector3D> points;

    for (int i = 1; i <= 12; i++) {
        if (i == 1) points.push_back(Vector3D::point(0, 0, sqrt(5) / 2));
        else if (2 <= i && i <= 6)
            points.push_back(Vector3D::point(cos((i - 2) * 2 * M_PI / 5), sin((i - 2) * 2 * M_PI / 5), 0.5));
        else if (7 <= i && i <= 11)
            points.push_back(
                    Vector3D::point(cos(M_PI / 5 + (i - 7) * 2 * M_PI / 5), sin(M_PI / 5 + (i - 7) * 2 * M_PI / 5),
                                    -0.5));
        else points.push_back(Vector3D::point(0, 0, -sqrt(5) / 2));
    }

    // Create faces
    std::vector<Face> faces;

    std::vector<int> point_indexes0 = {0, 1, 2};
    Face face0(point_indexes0);
    faces.push_back(face0);

    std::vector<int> point_indexes1 = {0, 2, 3};
    Face face1(point_indexes1);
    faces.push_back(face1);

    std::vector<int> point_indexes2 = {0, 3, 4};
    Face face2(point_indexes2);
    faces.push_back(face2);

    std::vector<int> point_indexes3 = {0, 4, 5};
    Face face3(point_indexes3);
    faces.push_back(face3);

    std::vector<int> point_indexes4 = {0, 5, 1};
    Face face4(point_indexes4);
    faces.push_back(face4);

    std::vector<int> point_indexes5 = {1, 6, 2};
    Face face5(point_indexes5);
    faces.push_back(face5);

    std::vector<int> point_indexes6 = {2, 6, 7};
    Face face6(point_indexes6);
    faces.push_back(face6);

    std::vector<int> point_indexes7 = {2, 7, 3};
    Face face7(point_indexes7);
    faces.push_back(face7);

    std::vector<int> point_indexes8 = {3, 7, 8};
    Face face8(point_indexes8);
    faces.push_back(face8);

    std::vector<int> point_indexes9 = {3, 8, 4};
    Face face9(point_indexes9);
    faces.push_back(face9);

    std::vector<int> point_indexes10 = {4, 8, 9};
    Face face10(point_indexes10);
    faces.push_back(face10);

    std::vector<int> point_indexes11 = {4, 9, 5};
    Face face11(point_indexes11);
    faces.push_back(face11);

    std::vector<int> point_indexes12 = {5, 9, 10};
    Face face12(point_indexes12);
    faces.push_back(face12);

    std::vector<int> point_indexes13 = {5, 10, 1};
    Face face13(point_indexes13);
    faces.push_back(face13);

    std::vector<int> point_indexes14 = {1, 10, 6};
    Face face14(point_indexes14);
    faces.push_back(face14);

    std::vector<int> point_indexes15 = {11, 7, 6};
    Face face15(point_indexes15);
    faces.push_back(face15);

    std::vector<int> point_indexes16 = {11, 8, 7};
    Face face16(point_indexes16);
    faces.push_back(face16);

    std::vector<int> point_indexes17 = {11, 9, 8};
    Face face17(point_indexes17);
    faces.push_back(face17);

    std::vector<int> point_indexes18 = {11, 10, 9};
    Face face18(point_indexes18);
    faces.push_back(face18);

    std::vector<int> point_indexes19 = {11, 6, 10};
    Face face19(point_indexes19);
    faces.push_back(face19);

    // Create figure
    Figure octahedron(points, faces, color, center, scale, angleX, angleY, angleZ);
    return octahedron;
}

Figure createDodecahedron(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ) {
    // Create points
    std::vector<Vector3D> points;

    Figure icosahedron = createIcosahedron(color, center, scale, angleX, angleY, angleZ);

    for (const Face &triangle: icosahedron.faces) {
        Vector3D point1 = icosahedron.points[triangle.point_indexes[0]];
        Vector3D point2 = icosahedron.points[triangle.point_indexes[1]];
        Vector3D point3 = icosahedron.points[triangle.point_indexes[2]];
        points.push_back((point1 + point2 + point3) / 3);
    }

    // Create faces
    std::vector<Face> faces;

    std::vector<int> point_indexes0 = {0, 1, 2, 3, 4};
    Face face0(point_indexes0);
    faces.push_back(face0);

    std::vector<int> point_indexes1 = {0, 5, 6, 7, 1};
    Face face1(point_indexes1);
    faces.push_back(face1);

    std::vector<int> point_indexes2 = {1, 7, 8, 9, 2};
    Face face2(point_indexes2);
    faces.push_back(face2);

    std::vector<int> point_indexes3 = {2, 9, 10, 11, 3};
    Face face3(point_indexes3);
    faces.push_back(face3);

    std::vector<int> point_indexes4 = {3, 11, 12, 13, 4};
    Face face4(point_indexes4);
    faces.push_back(face4);

    std::vector<int> point_indexes5 = {4, 13, 14, 5, 0};
    Face face5(point_indexes5);
    faces.push_back(face5);

    std::vector<int> point_indexes6 = {19, 18, 17, 16, 15};
    Face face6(point_indexes6);
    faces.push_back(face6);

    std::vector<int> point_indexes7 = {19, 14, 13, 12, 18};
    Face face7(point_indexes7);
    faces.push_back(face7);

    std::vector<int> point_indexes8 = {18, 12, 11, 10, 17};
    Face face8(point_indexes8);
    faces.push_back(face8);

    std::vector<int> point_indexes9 = {17, 10, 9, 8, 16};
    Face face9(point_indexes9);
    faces.push_back(face9);

    std::vector<int> point_indexes10 = {16, 8, 7, 6, 15};
    Face face10(point_indexes10);
    faces.push_back(face10);

    std::vector<int> point_indexes11 = {15, 6, 5, 14, 19};
    Face face11(point_indexes11);
    faces.push_back(face11);

    // Create figure
    Figure dodecahedron(points, faces, color, center, scale, angleX, angleY, angleZ);
    return dodecahedron;
}

Figure
createSphere(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, const int n) {
    Figure icosahedron = createIcosahedron(color, center, scale, angleX, angleY, angleZ);

    for (int _ = 0; _ < n; _++) {
        std::vector<Face> faces;
        std::vector<Vector3D> points;
        int triangleCounter = 0;
        for (Face triangle: icosahedron.faces) {
            // Fetch A, B and C
            Vector3D A = icosahedron.points[triangle.point_indexes[0]];
            Vector3D B = icosahedron.points[triangle.point_indexes[1]];
            Vector3D C = icosahedron.points[triangle.point_indexes[2]];

            // Calculate D, E and F
            Vector3D D = (A + B) / 2;
            Vector3D E = (A + C) / 2;
            Vector3D F = (B + C) / 2;

            // Add new points to temp points vector
            points.push_back(A); // 0 + 6*triangleCounter
            points.push_back(B); // 1 + 6*triangleCounter
            points.push_back(C); // 2 + 6*triangleCounter
            points.push_back(D); // 3 + 6*triangleCounter
            points.push_back(E); // 4 + 6*triangleCounter
            points.push_back(F); // 5 + 6*triangleCounter

            // Add new faces to temp faces vector
            faces.emplace_back(
                    std::vector<int>{6 * triangleCounter, 3 + 6 * triangleCounter, 4 + 6 * triangleCounter}); // ADE
            faces.emplace_back(
                    std::vector<int>{1 + 6 * triangleCounter, 5 + 6 * triangleCounter, 3 + 6 * triangleCounter}); // BFD
            faces.emplace_back(
                    std::vector<int>{2 + 6 * triangleCounter, 4 + 6 * triangleCounter, 5 + 6 * triangleCounter}); // CEF
            faces.emplace_back(
                    std::vector<int>{3 + 6 * triangleCounter, 5 + 6 * triangleCounter, 4 + 6 * triangleCounter}); // DFE
            triangleCounter++;
        }
        // Replace points and faces
        icosahedron.points = points;
        icosahedron.faces = faces;
    }

    for (Vector3D &point: icosahedron.points) {
        // Calculate radius
        double r = sqrt(pow(point.x, 2) + pow(point.y, 2) + pow(point.z, 2));

        // Replace coordinates
        point /= r;
    }
    return icosahedron;
}

Figure createCone(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, const int n,
                  const double h) {
    // Create points
    std::vector<Vector3D> points;
    for (int i = 0; i <= n; i++) {
        if (i == n) points.push_back(Vector3D::point(0, 0, h));
        else points.push_back(Vector3D::point(cos(2 * i * M_PI / n), sin(2 * i * M_PI / n), 0));
    }

    // Create faces
    std::vector<Face> faces;
    for (int i = 0; i <= n; i++) {
        if (i != n) faces.emplace_back(std::vector<int>{i, (i + 1) % n, n});
        else {
            std::vector<int> point_indexes;
            for (int j = n - 1; j >= 0; j--) point_indexes.push_back(j);
            faces.emplace_back(point_indexes);
        }
    }
    Figure cone = Figure(points, faces, color, center, scale, angleX, angleY, angleZ);
    return cone;
}

Figure
createCylinder(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, const int n,
               const double h) {
    // Create points
    std::vector<Vector3D> points;

    // Bottom surface
    points.reserve(2 * n);
    for (int i = 0; i < n; i++) points.push_back(Vector3D::point(cos(2 * i * M_PI / n), sin(2 * i * M_PI / n), 0));

    // Top surface
    for (int i = n; i < 2 * n; i++) points.push_back(Vector3D::point(cos(2 * i * M_PI / n), sin(2 * i * M_PI / n), h));

    // Create faces
    std::vector<Face> faces;

//    faces.reserve(n + 2);
    for (int i = 0; i <= n; i++) {
        if (i == n) {
            std::vector<int> point_indexes_bottom;
            for (int j = n - 1; j >= 0; j--) point_indexes_bottom.push_back(j);
            faces.emplace_back(point_indexes_bottom);
            std::vector<int> point_indexes_top;
            for (int j = 2 * n - 1; j >= n; j--) point_indexes_top.push_back(j);
            faces.emplace_back(point_indexes_top);
        } else if (i == n - 1) faces.emplace_back(std::vector<int>{i + 1, n + i, i, (i + 1) % n});
        else faces.emplace_back(std::vector<int>{n + i + 1, n + i, i, i + 1});
    }


    Figure cylinder = Figure(points, faces, color, center, scale, angleX, angleY, angleZ);
    return cylinder;
}

Figure
createTorus(Color color, Vector3D &center, double scale, double angleX, double angleY, double angleZ, const double r,
            const double R, const int n, const int m) {
    // Create points
    std::vector<Vector3D> points;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            // Calculate u and v
            double u = 2 * i * M_PI / n;
            double v = 2 * j * M_PI / m;

            // Determine x, y and z using the parametric equations
            double x_uv = (R + r * cos(v)) * cos(u);
            double y_uv = (R + r * cos(v)) * sin(u);
            double z_uv = r * sin(v);

            points.push_back(Vector3D::point(x_uv, y_uv, z_uv));
        }
    }

    // Create faces
    std::vector<Face> faces;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            faces.emplace_back(
                    std::vector<int>{m * i + j, ((i + 1) % n) * m + j, ((i + 1) % n) * m + ((j + 1) % m),
                                     i * m + ((j + 1) % m)});
        }
    }

    Figure torus = Figure(points, faces, color, center, scale, angleX, angleY, angleZ);
    return torus;
}

img::EasyImage generate_image(const ini::Configuration &configuration) {
//     ############################# INTRO #############################
//    img::EasyImage image((int) configuration["ImageProperties"]["width"],
//                         (int) configuration["ImageProperties"]["height"]);
//    ColorRectangle(image);
//    Blocks(image, configuration["ImageProperties"]["width"], configuration["ImageProperties"]["height"], configuration["BlockProperties"]["nrXBlocks"], configuration["BlockProperties"]["nrYBlocks"], configuration["BlockProperties"]["colorWhite"], configuration["BlockProperties"]["colorBlack"]);
//    QuarterCircle(image, configuration["ImageProperties"]["height"], configuration["ImageProperties"]["width"],
//                  configuration["LineProperties"]["nrLines"], configuration["LineProperties"]["lineColor"],
//                  configuration["LineProperties"]["backgroundcolor"]);
//    Eye(image, configuration["ImageProperties"]["height"], configuration["ImageProperties"]["width"],
//        configuration["LineProperties"]["nrLines"], configuration["LineProperties"]["lineColor"],
//        configuration["LineProperties"]["backgroundcolor"]);
//    Diamond(image, configuration["ImageProperties"]["height"], configuration["ImageProperties"]["width"],
//            configuration["LineProperties"]["nrLines"], configuration["LineProperties"]["lineColor"],
//            configuration["LineProperties"]["backgroundcolor"]);

//    ############################# 2D L-systems #############################
//    LParser::LSystem2D l_system;
//    std::ifstream input_stream(configuration["2DLSystem"]["inputfile"]);
//    input_stream >> l_system;
//    input_stream.close();
//    std::vector<double> color = configuration["2DLSystem"]["color"];
//    std::vector<double> bg_col = configuration["General"]["backgroundcolor"];
//    Color c(color[0], color[1], color[2]);
//    img::Color bg(bg_col[0] * 255, bg_col[1] * 255, bg_col[2] * 255);
//    img::EasyImage image = draw2DLines(drawLSystem(l_system, c), configuration["General"]["size"], bg);

//    ############################# 3D Line drawings #############################
//    std::vector<double> bg_col = configuration["General"]["backgroundcolor"];
//    img::Color bg(bg_col[0] * 255, bg_col[1] * 255, bg_col[2] * 255);
//    int size = configuration["General"]["size"];
//    Figures3D figures;
//    int nrFigures = configuration["General"]["nrFigures"];
//    for (int i = 0; i < nrFigures; i++) {
//        std::vector<Vector3D> points;
//        int nrPoints = configuration["Figure" + std::to_string(i)]["nrPoints"];
//        for (int indexp = 0; indexp < nrPoints; indexp++) {
//            std::vector<double> point = configuration["Figure" + std::to_string(i)]["point" + std::to_string(indexp)];
//            Vector3D p = Vector3D::point(point[0], point[1], point[2]);
//            points.push_back(p);
//        }
//
//        std::vector<Face> faces;
//        int nrLines = configuration["Figure" + std::to_string(i)]["nrLines"];
//        for (int indexl = 0; indexl < nrLines; indexl++) {
//            std::vector<int> indexes = configuration["Figure" + std::to_string(i)]["line" + std::to_string(indexl)];
//            Face f = Face(indexes);
//            faces.push_back(f);
//        }
//        std::vector<double> col = configuration["Figure" + std::to_string(i)]["color"];
//        Color color(col[0], col[1], col[2]);
//        std::vector<double> centerFetch = configuration["Figure" + std::to_string(i)]["center"];
//        Vector3D center = Vector3D::point(centerFetch[0], centerFetch[1], centerFetch[2]);
//        double scale = configuration["Figure" + std::to_string(i)]["scale"];
//
//        // Get rotation angles
//        double degreeX = configuration["Figure" + std::to_string(i)]["rotateX"];
//        double angleX = degreeX / 180 * M_PI;
//        double degreeY = configuration["Figure" + std::to_string(i)]["rotateY"];
//        double angleY = degreeY / 180 * M_PI;
//        double degreeZ = configuration["Figure" + std::to_string(i)]["rotateZ"];
//        double angleZ = degreeZ / 180 * M_PI;
//        Figure f(points, faces, color, center, scale, angleX, angleY, angleZ);
//        figures.push_back(f);
//    }
//    std::vector<double> eyepoint_ = configuration["General"]["eye"];
//    Vector3D eyepoint = Vector3D::point(eyepoint_[0], eyepoint_[1], eyepoint_[2]);
//    img::EasyImage image = draw2DLines(doProjection(figures, eyepoint), size, bg);

//    ############################# 3D Figures #############################
    std::vector<double> bg_col = configuration["General"]["backgroundcolor"];
    img::Color bg(bg_col[0] * 255, bg_col[1] * 255, bg_col[2] * 255);
    int size = configuration["General"]["size"];
    Figures3D figures;
    int nrFigures = configuration["General"]["nrFigures"];
    bool l_sys = false;
    img::EasyImage image;
    for (int i = 0; i < nrFigures; i++) {
        std::vector<double> col = configuration["Figure" + std::to_string(i)]["color"];
        Color color(col[0], col[1], col[2]);
        std::vector<double> centerFetch = configuration["Figure" + std::to_string(i)]["center"];
        Vector3D center = Vector3D::point(centerFetch[0], centerFetch[1], centerFetch[2]);
        double scale = configuration["Figure" + std::to_string(i)]["scale"];

        // Get rotation angles
        double degreeX = configuration["Figure" + std::to_string(i)]["rotateX"];
        double angleX = degreeX / 180 * M_PI;
        double degreeY = configuration["Figure" + std::to_string(i)]["rotateY"];
        double angleY = degreeY / 180 * M_PI;
        double degreeZ = configuration["Figure" + std::to_string(i)]["rotateZ"];
        double angleZ = degreeZ / 180 * M_PI;
        std::string type = configuration["Figure" + std::to_string(i)]["type"];
        if (type == "Cube") figures.push_back(createCube(color, center, scale, angleX, angleY, angleZ));
        else if (type == "Tetrahedron")
            figures.push_back(createTetrahedron(color, center, scale, angleX, angleY, angleZ));
        else if (type == "Icosahedron")
            figures.push_back(createIcosahedron(color, center, scale, angleX, angleY, angleZ));
        else if (type == "Octahedron")
            figures.push_back(createOctahedron(color, center, scale, angleX, angleY, angleZ));
        else if (type == "Dodecahedron")
            figures.push_back(createDodecahedron(color, center, scale, angleX, angleY, angleZ));
        else if (type == "Cone") {
            double height = configuration["Figure" + std::to_string(i)]["height"];
            int n = configuration["Figure" + std::to_string(i)]["n"];
            figures.push_back(createCone(color, center, scale, angleX, angleY, angleZ, n, height));
        }
        else if (type == "Cylinder") {
            double height = configuration["Figure" + std::to_string(i)]["height"];
            int n = configuration["Figure" + std::to_string(i)]["n"];
            figures.push_back(createCylinder(color, center, scale, angleX, angleY, angleZ, n, height));
        }
        else if (type == "Sphere") {
            int n = configuration["Figure" + std::to_string(i)]["n"];
            figures.push_back(createSphere(color, center, scale, angleX, angleY, angleZ, n));
        }
        else if (type == "Torus") {
            double r = configuration["Figure" + std::to_string(i)]["r"];
            double R = configuration["Figure" + std::to_string(i)]["R"];
            int n = configuration["Figure" + std::to_string(i)]["m"];
            int m = configuration["Figure" + std::to_string(i)]["n"];
            figures.push_back(createTorus(color, center, scale, angleX, angleY, angleZ, r, R, n, m));
        }
        else if (type == "3DLSystem") {
            l_sys = true;
            LParser::LSystem2D l_system;
            std::ifstream input_stream(configuration["Figure" + std::to_string(i)]["inputfile"]);
            input_stream >> l_system;
            input_stream.close();
            image = draw2DLines(drawLSystem(l_system, color), configuration["General"]["size"], bg);
        }
    }
    std::vector<double> eyepoint_ = configuration["General"]["eye"];
    Vector3D eyepoint = Vector3D::point(eyepoint_[0], eyepoint_[1], eyepoint_[2]);
    if (!l_sys)image = draw2DLines(doProjection(figures, eyepoint), size, bg);
    // Output image
    std::ofstream fout("out.bmp", std::ios::binary);
    fout << image;
    fout.close();
    return image;
}


int main(int argc, char const *argv[]) {
    int retVal = 0;
    try {
        std::vector<std::string> args = std::vector<std::string>(argv + 1, argv + argc);
        if (args.empty()) {
            std::ifstream fileIn("filelist");
            std::string filelistName;
            while (std::getline(fileIn, filelistName)) {
                args.push_back(filelistName);
            }
        }
        for (std::string fileName: args) {
            ini::Configuration conf;
            try {
                std::ifstream fin(fileName);
                fin >> conf;
                fin.close();
            }
            catch (ini::ParseException &ex) {
                std::cerr << "Error parsing file: " << fileName << ": " << ex.what() << std::endl;
                retVal = 1;
                continue;
            }

            img::EasyImage image = generate_image(conf);
            if (image.get_height() > 0 && image.get_width() > 0) {
                std::string::size_type pos = fileName.rfind('.');
                if (pos == std::string::npos) {
                    //filename does not contain a '.' --> append a '.bmp' suffix
                    fileName += ".bmp";
                } else {
                    fileName = fileName.substr(0, pos) + ".bmp";
                }
                try {
                    std::ofstream f_out(fileName.c_str(), std::ios::trunc | std::ios::out | std::ios::binary);
                    f_out << image;

                }
                catch (std::exception &ex) {
                    std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                    retVal = 1;
                }
            } else {
                std::cout << "Could not generate image for " << fileName << std::endl;
            }
        }
    }
    catch (const std::bad_alloc &exception) {
        //When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
        //Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
        //(Unless of course you are already consuming the maximum allowed amount of memory)
        //If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
        //mark the test as failed while in reality it just needed a bit more memory
        std::cerr << "Error: insufficient memory" << std::endl;
        retVal = 100;
    }
    return retVal;
}
