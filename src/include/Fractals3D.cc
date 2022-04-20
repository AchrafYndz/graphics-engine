#include "Fractals3D.h"
#include "Logic3D.h"

void generateFractal(Figure &fig, Figures3D &fractal, const int nr_iterations, const double scale) {
    fractal.push_back(fig);
    for (int i = 0; i < nr_iterations; i++) {
        Figures3D newFractals;
        for (Figure fracFig: fractal) {
            for (int j = 0; j < fig.points.size(); j++) {
                Figure copyFig = fracFig;
                const Matrix scaleMatrix = scaleFigure(1/scale);
                applyTransformation(copyFig, scaleMatrix);
                Matrix translationMatrix = translate(fracFig.points[j] - copyFig.points[j]);
                applyTransformation(copyFig, translationMatrix);
                newFractals.push_back(copyFig);
            }
        }
        fractal = newFractals;
    }
}