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

    Brackets(double x_, double y_, double angle_) : x(x_), y(y_), angle(angle_){};
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
    for (Figure& fig: figs) {
        Matrix m = eyePointTrans(eyepoint) * rotateX(fig.rotateAngleX) * rotateY(fig.rotateAngleY) *
                   rotateZ(fig.rotateAngleZ) *
                   translate(fig.center) * scaleFigure(fig.scale);
        applyTransformation(fig, m);
        for (Face& face: fig.faces) {
            Vector3D p0 = fig.points[face.point_indexes[0]];
            Vector3D p1 = fig.points[face.point_indexes[1]];
            Point2D x = doProjection(p0, 1);
            Point2D y = doProjection(p1, 1);
            projection.push_back(Line2D(x, y, fig.color));
        }
    }
    return projection;
}


img::EasyImage generate_image(const ini::Configuration &configuration) {
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
//    LParser::LSystem2D l_system;
//
//    std::ifstream input_stream(configuration["2DLSystem"]["inputfile"]);
//    input_stream >> l_system;
//    input_stream.close();
//    std::vector<double> color = configuration["2DLSystem"]["color"];
//    std::vector<double> bg_col = configuration["General"]["backgroundcolor"];
//    Color c(color[0], color[1], color[2]);
//    img::Color bg(bg_col[0]*255, bg_col[1]*255, bg_col[2]*255);
//    img::EasyImage image = draw2DLines(drawLSystem(l_system, c) , configuration["General"]["size"], bg);

    std::vector<double> bg_col = configuration["General"]["backgroundcolor"];
    img::Color bg(bg_col[0] * 255, bg_col[1] * 255, bg_col[2] * 255);
    int size = configuration["General"]["size"];
    Figures3D figures;
    int nrFigures = configuration["General"]["nrFigures"];
    for (int i = 0; i < nrFigures; i++) {
        std::vector<Vector3D> points;
        int nrPoints = configuration["Figure" + std::to_string(i)]["nrPoints"];
        for (int indexp = 0; indexp < nrPoints; indexp++) {
            std::vector<double> point = configuration["Figure" + std::to_string(i)]["point" + std::to_string(indexp)];
            Vector3D p = Vector3D::point(point[0], point[1], point[2]);
            points.push_back(p);
        }

        std::vector<Face> faces;
        int nrLines = configuration["Figure" + std::to_string(i)]["nrLines"];
        for (int indexl = 0; indexl < nrLines; indexl++) {
            std::vector<int> indexes = configuration["Figure" + std::to_string(i)]["line" + std::to_string(indexl)];
            Face f = Face(indexes);
            faces.push_back(f);
        }
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
        Figure f(points, faces, color, center, scale, angleX, angleY, angleZ);
        figures.push_back(f);
    }
    std::vector<double> eyepoint_ = configuration["General"]["eye"];
    Vector3D eyepoint = Vector3D::point(eyepoint_[0], eyepoint_[1], eyepoint_[2]);
    img::EasyImage image = draw2DLines(doProjection(figures, eyepoint), size, bg);
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
