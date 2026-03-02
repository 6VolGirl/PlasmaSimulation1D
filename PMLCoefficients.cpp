//
// Created by 6anna on 15.02.2026.
//

#include <cmath>
#include "PMLCoefficients.h"

PMLSigma::PMLSigma(int nx, int pmlCells, double eps, double mu, double damping, int power, double dx)
        : sigmaE(nx+1, 0.0), sigmaM(nx, 0.0)
{
    double eta = std::sqrt(mu/eps);
    double L = pmlCells * dx;
    double s = -std::log(damping) / (eta * L);
    double s_m = s * mu / eps;

    auto prof = [&](int n, int N) {
        double x = double(n) / double(N); // 0..1
        double val = 1.0;
        for(int k=0;k<power;k++) val *= x;
        return val;
    };

    for (int i = 0; i < pmlCells; ++i) {
        double g = prof(i, pmlCells);
        sigmaE[i] = prof(pmlCells-1-i, pmlCells) * s;
        if(i < nx) sigmaM[i] = prof(pmlCells-1-i, pmlCells) * s_m;

        // Right side
        int ie = (nx+1) - 1 - i;
        int ih = nx - 1 - i;
        sigmaE[ie] = prof(pmlCells-1-i, pmlCells) * s;
        sigmaM[ih] = prof(pmlCells-1-i, pmlCells) * s_m;
    }
}
