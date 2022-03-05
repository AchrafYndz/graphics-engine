//
// Created by Achraf on 04/03/2022.
//

#include "Line2D.h"

img::Color Line2D::getEzColor() {
    return img::Color(color.red*255, color.green*255, color.blue*255);
}
