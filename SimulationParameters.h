//
// Created by 6anna on 15.02.2026.
//

#ifndef SIMULATIONPARAMETERS_H
#define SIMULATIONPARAMETERS_H




struct SimulationParameters {
    int resolution;

    int nx;
    int numTimeSteps;
    double dx;
    double dt;
    double courantNumber;

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

    // Drude plasma (ADE)
    bool   useDrude = true;      // включить/выключить плазму
    int    plasmaStart = 150;    // индекс начала плазмы
    int    plasmaWidth = 100;    // ширина плазмы
    double omega_p = 1.0;        // плазменная частота
    double gamma   = 0.0;        // затухание
    double epsilon_inf = 1.0;
    double sigmaCond   = 0.0;    // доп. проводимость материала

};


#endif //SIMULATIONPARAMETERS_H
