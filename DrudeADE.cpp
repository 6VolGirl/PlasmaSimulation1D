

#include "DrudeADE.h"

void DrudeADE::resize(int nE) {
    Jn.assign(nE, 0.0);
    Jnm1.assign(nE, 0.0);
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

inline double DrudeADE::Jhalf_noEterm(int i) const {
    if (!active[i]) return 0.0;
    return 0.5 * ((1.0 + alpha[i]) * Jn[i] + xi[i] * Jnm1[i]);
}

void DrudeADE::updateJ(const std::vector<double>& Enp1, const std::vector<double>& Enm1, double dt)
{
    const double inv2dt = 1.0 / (2.0 * dt);
    for (int i = 0; i < (int)Jn.size(); ++i) {
        if (!active[i]) continue;
        const double dE = (Enp1[i] - Enm1[i]) * inv2dt;      // (E^{n+1}-E^{n-1})/(2dt) [file:3]
        const double Jnp1 = alpha[i] * Jn[i] + xi[i] * Jnm1[i] + gammaCoef[i] * dE;
        Jnm1[i] = Jn[i];
        Jn[i] = Jnp1;
    }
}