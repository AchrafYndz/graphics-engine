#ifndef ENGINE_LIGHT_H
#define ENGINE_LIGHT_H


#include <list>
#include "Color.h"

class Light {
public:
    Color ambientLight = Color(0, 0, 0);
    Color diffuseLight = Color(0, 0, 0);
    Color specularLight = Color(0, 0, 0);
};

typedef std::list<Light> Lights3D;

#endif //ENGINE_LIGHT_H
