

#include "DrudeADE.h"

void DrudeADE::resize(int nE) {
    Jn.assign(nE, 0.0);
    J_prev.assign(nE, 0.0);
    alpha.assign(nE, 0.0);
    xi.assign(nE, 0.0);
    gammaCoef.assign(nE, 0.0);
    active.assign(nE, 0);
}

void DrudeADE::configureRegion(int i0, int i1, double omegaP, double Gamma, double f, double dt) {
        i0 = std::max(i0, 0);
        i1 = std::min(i1, (int)Jn.size()-1);
        const double g = Gamma * dt * 0.5;

        // Drude => ωj = 0
        const double a = 2.0 / (1.0 + g);
        const double x = (g - 1.0) / (g + 1.0);
        const double c = (f * omegaP * omegaP * dt * dt) / (1.0 + g);

        for (int i = i0; i <= i1; ++i) {
            active[i] = 1;
            alpha[i] = a;
            xi[i] = x;
            gammaCoef[i] = c;
        }
    }

// inline double DrudeADE::Jhalf_noEterm(int i) const  {
//     if (!active[i]) return 0.0;
//     return 0.5 * ((1.0 + alpha[i]) * Jn[i] + xi[i] * J_prev[i]);
// }


void DrudeADE::updateJ(const std::vector<double>& Enp1, const std::vector<double>& Enm1, double dt)
{
    for (int i = 0; i < (int)Jn.size(); ++i) {
        if (!active[i]) continue;
        const double dE = (Enp1[i] - Enm1[i]) / (2.0 * dt);
        const double J_next = alpha[i] * Jn[i] + xi[i] * J_prev[i] + gammaCoef[i] * dE;
        J_prev[i] = Jn[i];
        Jn[i] = J_next;
    }
}