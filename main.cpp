#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "SimulationParameters.h"
#include "FDTD1D.h"



int main() {
    SimulationParameters params;

    // Сетка и шаги (нормированные)

    params.resolution = 20;
    params.nx = 400;
    params.dx = 1.0 / params.resolution;
    params.courantNumber = 0.5;
    params.dt = params.courantNumber * params.dx;
    params.numTimeSteps = 2000;


    // Базовый материал (вне плазмы)
    params.eps0 = 1.0;
    params.mu0  = 1.0;


    // PML
    params.pmlThickness = 20;     // в ячейках
    params.pmlDamping = 1e-9;
    params.pmlProfilePower = 3;


    // Источник GaussianSource
    double wavelength_min = 0.3;
    double wavelength_max = 0.8;
    params.sourceFreq = 2.0 / (wavelength_min + wavelength_max);
    params.sourceFWidth = 1.0 / wavelength_min - 1.0 / wavelength_max;

    params.source_pos = 50;          // не в PML

    // Плазма Друде (ADE)
    params.useDrude = true;
    params.plasmaStart = 100;
    params.plasmaWidth = 100;      // ширина
    params.epsilon_inf = 1.0;      // ε∞ в модели Друде
    params.omega_p = 50;          // пример! (в нормированных единицах)
    params.gamma = 5;            // пример! (в нормированных единицах)
    params.sigmaCond = 0.0;        // если хочешь доп. омические потери в плазме

    // // Мониторы
    // params.monitorFront = 80;
    // params.monitorBack  = 300;



    try {
        FDTD1D_PythonStyle sim(params);
        sim.run();

        sim.writeImpulsePlasmaCSV("ImpulsePlasma.csv");

        std::cout << "Simulation finished.\n";
        std::cout << "Written: ImpulsePlasma.csv\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
