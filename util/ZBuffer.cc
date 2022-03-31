#include "ZBuffer.h"

#include <iostream>
#include <limits>

ZBuffer::ZBuffer(const int width, const int height) {
    for (int i=0; i<=width; i++) {
        std::vector<double> row;
        for (int j=0; j<=height; j++) {
            row.push_back(std::numeric_limits<double>::infinity());
        }
        this->push_back(row);
    }
}
