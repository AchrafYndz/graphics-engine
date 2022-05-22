#ifndef ENGINE_POINTLIGHT_H
#define ENGINE_POINTLIGHT_H


#include "../../util/Vector3D/vector3d.h"
#include "Light.h"

class PointLight: public Light {
public:
    Vector3D location;
    double spotAngle;
};


#endif //ENGINE_POINTLIGHT_H
