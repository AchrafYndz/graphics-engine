#ifndef ENGINE_INFLIGHT_H
#define ENGINE_INFLIGHT_H


#include "../../util/Vector3D/vector3d.h"
#include "Light.h"

class InfLight : public Light {
public:
    Vector3D ldVector;
};


#endif //ENGINE_INFLIGHT_H
