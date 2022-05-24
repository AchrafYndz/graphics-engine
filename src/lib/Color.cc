#include "Color.h"

Color &Color::operator+=(Color &rhs) {
    this->blue += rhs.blue;
    this->green += rhs.green;
    this->red += rhs.red;

    if (this->blue > 1) this->blue = 1;
    if (this->green > 1) this->green = 1;
    if (this->red > 1) this->red = 1;

    return *this;
}
