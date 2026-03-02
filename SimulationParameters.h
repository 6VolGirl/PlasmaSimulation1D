//
// Created by 6anna on 15.02.2026.
//

#ifndef SIMULATIONPARAMETERS_H
#define SIMULATIONPARAMETERS_H




struct SimulationParameters {
    int nx;
    int numTimeSteps;
    double dx;
    double dt;
    double courantNumber;

    // ДОБАВИТЬ:
    int resolution;   // число ячеек на единицу длины (как в python)

    // Material base
    double eps0 = 1.0;
    double mu0  = 1.0;

    // PML
    int pmlThickness;
    double pmlDamping = 1e-9;
    int pmlProfilePower = 3;

    // Source
    int source_pos;
    double sourceFreq;
    double sourceFWidth;

    // Monitors
    int monitorFront;
    int monitorBack;
};


#endif //SIMULATIONPARAMETERS_H
