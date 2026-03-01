//
// Created by 6anna on 15.02.2026.
//

#ifndef SIMULATIONPARAMETERS_H
#define SIMULATIONPARAMETERS_H



struct SimulationParameters {
    double a;          //Характерная длина [м]
    int resolution;   // количество ячеек в ячейке

    // Grid and time
    int nx, ny, nz;
    int numTimeSteps;
    double dx, dy, dz;
    double dt;
    double courantNumber;

    // Plasma (Drude)
    double omega_p;
    double gamma;
    double epsilon_inf;
    double sigmaCond;

    // Source
    int source_pos;
    double sourceFreqCenter;
    double sourceFreqWidth;
    double pulseWidth;

    // PML
    int pmlThickness;
    double pmlReflectCoeff;  // Теоретический
    double pmlGrading_m;

    // Monitors
    int monitorFront;
    int monitorBack;
    int plasmaStart;
    int plasmaWidth;
};



#endif //SIMULATIONPARAMETERS_H
