#include "Figure.h"


img::Color Figure::getEzColor() {
    return img::Color(color.red*255, color.green*255, color.blue*255);
}