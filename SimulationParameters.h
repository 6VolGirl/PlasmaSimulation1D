//
// Created by 6anna on 15.02.2026.
//

#ifndef SIMULATIONPARAMETERS_H
#define SIMULATIONPARAMETERS_H




struct SimulationParameters {
    int resolution;
    int nx;
    int numTimeSteps;
    double dx, dt, courantNumber;

    // Base
    double epsInf0;
    double mu0;

    // PML
    int pmlThickness;
    double pmlDamping = 1e-9;
    int pmlProfilePower = 3;

    // Source
    int source_pos;
    double sourceFreq;
    double sourceFWidth;

    // Drude ADE
    bool useDrude = false;
    int drudeStart = 0;
    int drudeEnd   = -1;
    double drudeOmegaP = 0.0;   // ωp
    double drudeGamma  = 0.0;   // Γ
    double drudeStrength = 1.0; // f
};



#endif //SIMULATIONPARAMETERS_H
