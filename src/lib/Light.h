#ifndef ENGINE_LIGHT_H
#define ENGINE_LIGHT_H


#include <list>
#include "Color.h"

class Light {
public:
    Color ambientLight;
    Color diffuseLight;
    Color specularLight;
};

typedef std::list<Light> Lights3D;

#endif //ENGINE_LIGHT_H
