#ifndef ENGINE_FACE_H
#define ENGINE_FACE_H

#include <vector>

class Face {
public:
    std::vector<int> point_indexes;
    Face() = default;
    Face(std::vector<int> point_indexes_) : point_indexes(point_indexes_) {};
};


#endif //ENGINE_FACE_H
