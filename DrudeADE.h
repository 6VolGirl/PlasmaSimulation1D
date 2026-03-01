//
// Created by 6anna on 15.02.2026.
//

#ifndef DRUDEADE_H
#define DRUDEADE_H
#include <vector>
#include "SimulationParameters.h"


class DrudeADE {
private:
    int gridSize_;
    std::vector<double> alpha_, xi_, gamma_;  // ADE
    std::vector<double> J_curr_, J_prev_;     // Поляризация: n and n-1

    // E-поля update коэфф(зависят от материала + шага по времени)
    double C1_, C2_, C3_;

public:
    DrudeADE(int grid_size, const SimulationParameters& params, double dt);

    void updatePolarizationCurrent(const std::vector<double>& E_new, const std::vector<double>& E_old, double dt);
    void stepHistory();

    // Поправочный коэффициент
    double getPolarizationCorrection(int i) const;

    double getC1() const { return C1_; }
    double getC2() const { return C2_; }
    double getC3() const { return C3_; }
    const std::vector<double>& getJ() const { return J_curr_; }
    double getGamma(int i) const { return gamma_[i]; }
};



#endif //DRUDEADE_H
