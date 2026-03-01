//
// Created by 6anna on 15.02.2026.
//

#ifndef PMLCOEFFICIENTS_H
#define PMLCOEFFICIENTS_H

#include <vector>


class PMLCoefficients {
public:
    std::vector<double> gi1, gi2, gi3;  // E-поле curl
    std::vector<double> fi1, fi2, fi3;  // H-поле curl

    PMLCoefficients(int gridSize, int pmlThickness, double pmlCoeff = 0.33) {
        gi1.resize(gridSize, 0.0);
        gi2.resize(gridSize, 1.0);
        gi3.resize(gridSize, 1.0);
        fi1.resize(gridSize, 0.0);
        fi2.resize(gridSize, 1.0);
        fi3.resize(gridSize, 1.0);

        for (int n = 0; n < pmlThickness; ++n) {
            double xxn = (double)(pmlThickness - n) / pmlThickness;
            double xn = pmlCoeff * xxn * xxn * xxn;

            fi1[n] = xn;
            gi2[n] = 1.0 / (1.0 + xn);
            gi3[n] = (1.0 - xn) / (1.0 + xn);

            fi1[gridSize - n - 1] = xn;
            gi2[gridSize - n - 1] = 1.0 / (1.0 + xn);
            gi3[gridSize - n - 1] = (1.0 - xn) / (1.0 + xn);

            // Half-step versions
            xxn = (double)(pmlThickness - n - 0.5) / pmlThickness;
            xn = pmlCoeff * xxn * xxn * xxn;
            gi1[n] = xn;
            fi2[n] = 1.0 / (1.0 + xn);
            fi3[n] = (1.0 - xn) / (1.0 + xn);

            gi1[gridSize - n - 1] = xn;
            fi2[gridSize - n - 1] = 1.0 / (1.0 + xn);
            fi3[gridSize - n - 1] = (1.0 - xn) / (1.0 + xn);
        }
    }
};



#endif //PMLCOEFFICIENTS_H
