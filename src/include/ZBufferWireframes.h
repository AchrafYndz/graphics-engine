#ifndef ENGINE_ZBUFFERWIREFRAMES_H
#define ENGINE_ZBUFFERWIREFRAMES_H

#include "../../util/Easy_image/easy_image.h"

class ZBuffer : public std::vector<std::vector<double>>{
public:
    ZBuffer(int width, int height);
};

void
draw_zbuf_line(ZBuffer &zbuffer, img::EasyImage &image, unsigned int x0, unsigned int y0, double z0,
               unsigned int x1, unsigned int y1, double z1, img::Color color);


#endif //ENGINE_ZBUFFERWIREFRAMES_H
