//
// Created by 6anna on 15.02.2026.
//

#include "DrudeADE.h"


DrudeADE::DrudeADE(int grid_size, const SimulationParameters& params, double dt)
    : gridSize_(grid_size) {

    alpha_.resize(grid_size, 1.0);
    xi_.resize(grid_size, 0.0);
    gamma_.resize(grid_size, 0.0);
    J_curr_.resize(grid_size, 0.0);
    J_prev_.resize(grid_size, 0.0);

    double omega_p = params.omega_p;
    double gamma_coeff = params.gamma;
    double dt2 = dt * dt;
    double denom = 1.0 + gamma_coeff * dt / 2.0;

    int plasma_end = params.plasmaStart + params.plasmaWidth;
    for (int i = params.plasmaStart; i < plasma_end; ++i) {
        alpha_[i] = 2.0 / denom;  // ω² = 0 для Drude
        xi_[i] = (gamma_coeff * dt / 2.0 - 1.0) / denom;
        gamma_[i] = omega_p * omega_p * dt2 / denom;
    }
}

void DrudeADE::updatePolarizationCurrent(const std::vector<double>& E_new, const std::vector<double>& E_old, double dt) {
    for (int i = 0; i < gridSize_; ++i) {
        double dE_dt = (E_new[i] - E_old[i]) / (2.0 * dt);
        J_curr_[i] = alpha_[i] * J_curr_[i] + xi_[i] * J_prev_[i] + gamma_[i] * dE_dt;
    }
}


void DrudeADE::stepHistory() {
    std::swap(J_curr_, J_prev_);
}


double DrudeADE::getPolarizationCorrection(int i) const {
    return 0.5 * ((1.0 + alpha_[i]) * J_curr_[i] + xi_[i] * J_prev_[i]);
}

