#define _USE_MATH_DEFINES // Windows only

#include "ZBufferWireframes.h"

#include <limits>
#include <assert.h>
#include <cmath>

ZBuffer::ZBuffer(const int width, const int height) {
    for (int i = 0; i <= width; i++) {
        std::vector<double> row;
        for (int j = 0; j <= height; j++) {
            row.push_back(std::numeric_limits<double>::infinity());
        }
        this->push_back(row);
    }
}

void
draw_zbuf_line(ZBuffer &zbuffer, img::EasyImage &image, unsigned int x0, unsigned int y0, const double z0,
               unsigned int x1, unsigned int y1, const double z1, img::Color color) {
    assert(x0 < image.get_width() && y0 < image.get_height());
    assert(x1 < image.get_width() && y1 < image.get_height());
    if (x0 == x1) {
        //special case for x0 == x1
        for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++) {
            // Define x and y values
            double x = x0;
            double y = i;
            // Fetch a
            double a = ((double) std::max(y0, y1));
            // Calculate p and Zp
            double p = i / a;
            double Zp = p / z0 + (1 - p) / z1;
            // Zbuffer logic
            if (Zp < zbuffer[x][y]) {
                (image)(x, y) = color;
                zbuffer[x][y] = Zp;
            }
        }
    } else if (y0 == y1) {
        //special case for y0 == y1
        for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++) {
            // Define x and y values
            double x = i;
            double y = y0;
            // Fetch a
            double a = ((double) std::max(x0, x1));
            // Calculate p and Zp
            double p = i / a;
            double Zp = p / z0 + (1 - p) / z1;
            // Zbuffer logic
            if (Zp < zbuffer[x][y]) {
                (image)(x, y) = color;
                zbuffer[x][y] = Zp;
            }
        }
    } else {
        if (x0 > x1) {
            //flip points if x1>x0: we want x0 to have the lowest value
            std::swap(x0, x1);
            std::swap(y0, y1);
        }
        double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
        if (-1.0 <= m && m <= 1.0) {
            for (unsigned int i = 0; i <= (x1 - x0); i++) {
                // Define x and y values
                double x = x0 + i;
                double y = round(y0 + m * i);
                // Fetch a
                double a = ((double) (x1 - x0));
                // Calculate p and Zp
                double p = i / a;
                double Zp = p / z0 + (1 - p) / z1;
                // Zbuffer logic
                if (Zp < zbuffer[x][y]) {
                    (image)(x, y) = color;
                    zbuffer[x][y] = Zp;
                }
            }
        } else if (m > 1.0) {
            for (unsigned int i = 0; i <= (y1 - y0); i++) {
                // Define x and y values
                double x = (unsigned int) round(x0 + (i / m));
                double y = y0 + i;
                // Fetch a
                double a = ((double) (y1 - y0));
                // Calculate p and Zp
                double p = i / a;
                double Zp = p / z0 + (1 - p) / z1;
                // Zbuffer logic
                if (Zp < zbuffer[x][y]) {
                    (image)(x, y) = color;
                    zbuffer[x][y] = Zp;
                }
            }
        } else if (m < -1.0) {
            for (unsigned int i = 0; i <= (y0 - y1); i++) {
                // Define x and y values
                double x = (unsigned int) round(x0 - (i / m));
                double y = y0 - i;
                // Fetch a
                double a = ((double) (y1 - y0));
                // Calculate p and Zp
                double p = i / a;
                double Zp = p / z0 + (1 - p) / z1;
                // Zbuffer logic
                if (Zp < zbuffer[x][y]) {
                    (image)(x, y) = color;
                    zbuffer[x][y] = Zp;
                }
            }
        }
    }
}

