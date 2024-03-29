#ifndef ENGINE_COLOR_H
#define ENGINE_COLOR_H


class Color {
public:
    double red;
    double green;
    double blue;
    Color(double r, double g, double b): red(r), green(g), blue(b) {};
    Color& operator+=(Color& rhs);
    Color& operator*=(Color& rhs);
};


#endif //ENGINE_COLOR_H
