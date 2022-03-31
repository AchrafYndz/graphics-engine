#ifndef ENGINE_ZBUFFER_H
#define ENGINE_ZBUFFER_H

#include <vector>

class ZBuffer : public std::vector<std::vector<double>>{
public:
    ZBuffer(const int width, const int height);
};


#endif //ENGINE_ZBUFFER_H
