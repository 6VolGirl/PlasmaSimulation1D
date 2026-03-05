 //
 // Created by 6anna on 15.02.2026.
 //

 #ifndef DRUDEADE_H
 #define DRUDEADE_H
 #include <vector>
#include <cstdint>
 #include "SimulationParameters.h"

class DrudeADE {

    std::vector<double> Jn, Jnm1;
    std::vector<double> alpha, xi, gammaCoef;
    std::vector<uint8_t> active;

public:
    explicit DrudeADE(int nE = 0) { resize(nE); }

    void resize(int nE);

    void configureRegion(int i0, int i1, double omegaP, double Gamma, double f, double dt);

    inline double gammaSum(int i) const { return active[i] ? gammaCoef[i] : 0.0; }

    inline double Jhalf_noEterm(int i) const;

    void updateJ(const std::vector<double>& Enp1, const std::vector<double>& Enm1, double dt);

};

 #endif //DRUDEADE_H
