#ifndef ENGINE_INTRO_H
#define ENGINE_INTRO_H

#include "../../util/Easy_image/easy_image.h"

void ColorRectangle(img::EasyImage &img);

void
Blocks(img::EasyImage &img, int Wi, int Hi, int Nx, int Ny, std::vector<double> color1, std::vector<double> color2);

void QuarterCircle(img::EasyImage &img, int Hi, int Wi, int N, std::vector<double> lineColor,
                   const std::vector<double> &backgroundColor, const std::string &figure = "");

void
Eye(img::EasyImage &img, int Hi, int Wi, int N, std::vector<double> lineColor,
    const std::vector<double> &backgroundColor);

void Diamond(img::EasyImage &img, int Hi, int Wi, int N, std::vector<double> lineColor,
             const std::vector<double> &backgroundColor);

#endif //ENGINE_INTRO_H
