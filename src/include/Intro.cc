#include "Intro.h"

void ColorRectangle(img::EasyImage &img) {
    for (unsigned int i = 0; i < 500; i++) {
        for (unsigned int j = 0; j < 500; j++) {
            img(i, j).red = i * 256 / 500;
            img(i, j).green = j * 256 / 500;
            img(i, j).blue = (i * 256 / 500 + j * 256 / 500) % 255;
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
                   const std::vector<double> &backgroundColor, const std::string &figure) {
    // TODO: make sure the default parameter for figure is correct.
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
Eye(img::EasyImage &img, int Hi, int Wi, int N, std::vector<double> lineColor,
    const std::vector<double> &backgroundColor) {
    QuarterCircle(img, Hi, Wi, N, std::move(lineColor), backgroundColor, "eye");
}

void Diamond(img::EasyImage &img, int Hi, int Wi, int N, std::vector<double> lineColor,
             const std::vector<double> &backgroundColor) {
    QuarterCircle(img, Hi, Wi, N, lineColor, backgroundColor, "diamond");
}
