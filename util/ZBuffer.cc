#include "ZBuffer.h"

#include <iostream>
#include <limits>

ZBuffer::ZBuffer(const int width, const int height) {
    for (int i=0; i<=width; i++) {
        for (int j=0; j<=height; j++) {
            this->at(i).push_back(std::numeric_limits<double>::infinity());
        }
    }
    std::cout << std::endl;
}
