//
// Created by 6anna on 15.02.2026.
//

#ifndef PMLCOEFFICIENTS_H
#define PMLCOEFFICIENTS_H

#include <vector>


class PMLSigma {
public:
    std::vector<double> sigmaE; // size nx+1
    std::vector<double> sigmaM; // size nx

    PMLSigma(int nx, int pmlCells, double eps, double mu, double damping, int power, double dx);
};




#endif //PMLCOEFFICIENTS_H
