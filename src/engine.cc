#define _USE_MATH_DEFINES // Windows only

#include "../util/Easy_image/easy_image.h"
#include "../util/Ini_config/ini_configuration.h"

#include "lib/Line2D.h"
#include "lib/Face.h"

#include "include/Intro.h"
#include "include/LSystems2D.h"
#include "include/LSystems3D.h"
#include "include/Logic3D.h"
#include "include/Figures3D.h"
#include "include/ZBufferTriangles.h"
#include "include/Fractals3D.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <cmath>

img::EasyImage generate_image(const ini::Configuration &configuration) {
    if (((std::string) configuration["General"]["type"]).find("Intro") != std::string::npos) {
        //  ############################# INTRO ####################################
        img::EasyImage image((int) configuration["ImageProperties"]["width"],
                             (int) configuration["ImageProperties"]["height"]);
        if ((std::string) configuration["General"]["type"] == "IntroColorRectangle") {
            ColorRectangle(image);
        } else if (((std::string) configuration["General"]["type"] == "IntroBlocks")) {
            Blocks(image, configuration["ImageProperties"]["width"], configuration["ImageProperties"]["height"],
                   configuration["BlockProperties"]["nrXBlocks"], configuration["BlockProperties"]["nrYBlocks"],
                   configuration["BlockProperties"]["colorWhite"], configuration["BlockProperties"]["colorBlack"]);
        } else if (((std::string) configuration["General"]["type"] == "IntroLines")) {
            if ((std::string) configuration["LineProperties"]["figure"] == "QuarterCircle") {
                QuarterCircle(image, configuration["ImageProperties"]["height"],
                              configuration["ImageProperties"]["width"],
                              configuration["LineProperties"]["nrLines"], configuration["LineProperties"]["lineColor"],
                              configuration["LineProperties"]["backgroundcolor"]);
            } else if ((std::string) configuration["LineProperties"]["figure"] == "Eye") {
                Eye(image, configuration["ImageProperties"]["height"], configuration["ImageProperties"]["width"],
                    configuration["LineProperties"]["nrLines"], configuration["LineProperties"]["lineColor"],
                    configuration["LineProperties"]["backgroundcolor"]);
            } else if ((std::string) configuration["LineProperties"]["figure"] == "Diamond") {
                Diamond(image, configuration["ImageProperties"]["height"], configuration["ImageProperties"]["width"],
                        configuration["LineProperties"]["nrLines"], configuration["LineProperties"]["lineColor"],
                        configuration["LineProperties"]["backgroundcolor"]);
            }
        }
        std::ofstream fout("out.bmp", std::ios::binary);
        fout << image;
        fout.close();
        return image;
    } else if ((std::string) configuration["General"]["type"] == "2DLSystem") {
        //    ############################# 2D L-systems #############################
        LParser::LSystem2D l_system;
        std::ifstream input_stream(configuration["2DLSystem"]["inputfile"]);
        input_stream >> l_system;
        input_stream.close();
        std::vector<double> color = configuration["2DLSystem"]["color"];
        std::vector<double> bg_col = configuration["General"]["backgroundcolor"];
        Color c(color[0], color[1], color[2]);
        img::Color bg(bg_col[0] * 255, bg_col[1] * 255, bg_col[2] * 255);
        img::EasyImage image = draw2DLines(draw2DLSystem(l_system, c), configuration["General"]["size"], bg, false);
        std::ofstream fout("out.bmp", std::ios::binary);
        fout << image;
        fout.close();
        return image;
    } else if ((std::string) configuration["General"]["type"] == "Wireframe" ||
               (std::string) configuration["General"]["type"] == "ZBufferedWireframe" ||
               (std::string) configuration["General"]["type"] == "ZBuffering") {
        std::vector<double> bg_col = configuration["General"]["backgroundcolor"];
        img::Color bg(bg_col[0] * 255, bg_col[1] * 255, bg_col[2] * 255);
//        img::Color bg(0, 255, 0);
        int size = configuration["General"]["size"];
        bool zBuffer = ((std::string) configuration["General"]["type"] == "ZBufferedWireframe");
        bool toTriangulate = ((std::string) configuration["General"]["type"] == "ZBuffering");
        Figures3D figures;
        int nrFigures = configuration["General"]["nrFigures"];
        for (int i = 0; i < nrFigures; i++) {
            //    ############################# 3D Line drawings #############################
            if ((std::string) configuration["Figure" + std::to_string(i)]["type"] == "LineDrawing") {
                std::vector<Vector3D> points;
                int nrPoints = configuration["Figure" + std::to_string(i)]["nrPoints"];
                for (int indexp = 0; indexp < nrPoints; indexp++) {
                    std::vector<double> position = configuration["Figure" + std::to_string(i)]["point" +
                                                                                               std::to_string(indexp)];
                    Vector3D p = Vector3D::point(position[0], position[1], position[2]);
                    points.push_back(p);
                }

                std::vector<Face> faces;
                int nrLines = configuration["Figure" + std::to_string(i)]["nrLines"];
                for (int indexl = 0; indexl < nrLines; indexl++) {
                    std::vector<int> indexes = configuration["Figure" + std::to_string(i)]["line" +
                                                                                           std::to_string(indexl)];
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
            } else {
                //    ############################# 3D Figures #############################
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
                if (type == "Cube")
                    figures.push_back(createCube(color, center, scale, angleX, angleY, angleZ, toTriangulate));
                else if (type == "Tetrahedron")
                    figures.push_back(createTetrahedron(color, center, scale, angleX, angleY, angleZ, toTriangulate));
                else if (type == "Icosahedron")
                    figures.push_back(createIcosahedron(color, center, scale, angleX, angleY, angleZ, toTriangulate));
                else if (type == "Octahedron")
                    figures.push_back(createOctahedron(color, center, scale, angleX, angleY, angleZ, toTriangulate));
                else if (type == "Dodecahedron")
                    figures.push_back(createDodecahedron(color, center, scale, angleX, angleY, angleZ, toTriangulate));
                else if (type == "Cone") {
                    double height = configuration["Figure" + std::to_string(i)]["height"];
                    int n = configuration["Figure" + std::to_string(i)]["n"];
                    figures.push_back(
                            createCone(color, center, scale, angleX, angleY, angleZ, n, height, toTriangulate));
                } else if (type == "Cylinder") {
                    double height = configuration["Figure" + std::to_string(i)]["height"];
                    int n = configuration["Figure" + std::to_string(i)]["n"];
                    figures.push_back(
                            createCylinder(color, center, scale, angleX, angleY, angleZ, n, height, toTriangulate));
                } else if (type == "Sphere") {
                    int n = configuration["Figure" + std::to_string(i)]["n"];
                    figures.push_back(createSphere(color, center, scale, angleX, angleY, angleZ, n, toTriangulate));
                } else if (type == "Torus") {
                    double r = configuration["Figure" + std::to_string(i)]["r"];
                    double R = configuration["Figure" + std::to_string(i)]["R"];
                    int n = configuration["Figure" + std::to_string(i)]["m"];
                    int m = configuration["Figure" + std::to_string(i)]["n"];
                    figures.push_back(
                            createTorus(color, center, scale, angleX, angleY, angleZ, r, R, n, m, toTriangulate));
                } else if (type == "FractalTetrahedron") {
                    int nrIterations = configuration["Figure" + std::to_string(i)]["nrIterations"];
                    Figure tetrahedron = createTetrahedron(color, center, scale, angleX, angleY, angleZ, toTriangulate);
                    Figures3D fractal;
                    int fractalScale = configuration["Figure" + std::to_string(i)]["fractalScale"];
                    generateFractal(tetrahedron, fractal, nrIterations, fractalScale);
                    figures.insert(figures.end(), fractal.begin(), fractal.end());
                } else if (type == "FractalCube") {
                    int nrIterations = configuration["Figure" + std::to_string(i)]["nrIterations"];
                    Figure cube = createCube(color, center, scale, angleX, angleY, angleZ, toTriangulate);
                    Figures3D fractal;
                    int fractalScale = configuration["Figure" + std::to_string(i)]["fractalScale"];
                    generateFractal(cube, fractal, nrIterations, fractalScale);
                    figures.insert(figures.end(), fractal.begin(), fractal.end());
                } else if (type == "FractalIcosahedron") {
                    int nrIterations = configuration["Figure" + std::to_string(i)]["nrIterations"];
                    Figure icosahedron = createIcosahedron(color, center, scale, angleX, angleY, angleZ, toTriangulate);
                    Figures3D fractal;
                    int fractalScale = configuration["Figure" + std::to_string(i)]["fractalScale"];
                    generateFractal(icosahedron, fractal, nrIterations, fractalScale);
                    figures.insert(figures.end(), fractal.begin(), fractal.end());
                } else if (type == "FractalOctahedron") {
                    int nrIterations = configuration["Figure" + std::to_string(i)]["nrIterations"];
                    Figure octahedron = createOctahedron(color, center, scale, angleX, angleY, angleZ, toTriangulate);
                    Figures3D fractal;
                    int fractalScale = configuration["Figure" + std::to_string(i)]["fractalScale"];
                    generateFractal(octahedron, fractal, nrIterations, fractalScale);
                    figures.insert(figures.end(), fractal.begin(), fractal.end());
                } else if (type == "FractalDedocahedron") {
                    int nrIterations = configuration["Figure" + std::to_string(i)]["nrIterations"];
                    Figure dodecahedron = createDodecahedron(color, center, scale, angleX, angleY, angleZ, toTriangulate);
                    Figures3D fractal;
                    int fractalScale = configuration["Figure" + std::to_string(i)]["fractalScale"];
                    generateFractal(dodecahedron, fractal, nrIterations, fractalScale);
                    figures.insert(figures.end(), fractal.begin(), fractal.end());
                } else if (type == "3DLSystem") {
                    LParser::LSystem3D l_system;
                    std::ifstream input_stream(configuration["Figure" + std::to_string(i)]["inputfile"]);
                    input_stream >> l_system;
                    input_stream.close();
                    figures.push_back(draw3DLSystem(l_system, center, color, scale, angleX, angleY, angleZ));
                }
            }
        }
        std::vector<double> eyepoint_ = configuration["General"]["eye"];
        Vector3D eyepoint = Vector3D::point(eyepoint_[0], eyepoint_[1], eyepoint_[2]);
        img::EasyImage image;
        if (toTriangulate) {
            double d;
            double dx;
            double dy;
            double width;
            double height;
            getImageSpecs(doProjection(figures, eyepoint), size, d, dx, dy, width, height);

            ZBuffer zbuf(width, height);
            image = img::EasyImage(width, height, bg);
            for (Figure &figure: figures) {
                for (const Face &triangle: figure.faces) {
                    draw_zbuf_trag(zbuf, image, figure.points[triangle.point_indexes[0]],
                                   figure.points[triangle.point_indexes[1]],
                                   figure.points[triangle.point_indexes[2]],
                                   d, dx, dy,
                                   figure.getEzColor());
                }
            }
        } else image = draw2DLines(doProjection(figures, eyepoint), size, bg, zBuffer);
        std::ofstream fout("out.bmp", std::ios::binary);
        fout << image;
        fout.close();
        return image;
    }
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
