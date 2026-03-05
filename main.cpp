#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "SimulationParameters.h"
#include "FDTD1D.h"


#include <iostream>
#include <exception>

int main() {
    SimulationParameters params;

    params.resolution = 20;

    // Сетка и шаги
    params.nx = 400;
    params.dx = 1.0 / params.resolution;
    params.courantNumber = 0.5;
    params.dt = params.courantNumber * params.dx;
    params.numTimeSteps = 1000;

    // Материал
    params.epsInf = 1.0;
    params.mu0 = 1.0;

    // PML
    params.pmlThickness = 20;
    params.pmlDamping = 1e-9;
    params.pmlProfilePower = 3;

    // Источник
    double wavelength_min = 0.3;
    double wavelength_max = 0.8;
    params.sourceFreq = 2.0 / (wavelength_min + wavelength_max);     // cycles/time
    params.sourceFWidth = 1.0 / wavelength_min - 1.0 / wavelength_max; // cycles/time
    params.source_pos = 50;

    // Drude ADE
    params.useDrude = true;
    params.plasmaStart = 100;
    params.plasmaEnd = 200;

    params.оmega_p = 8.0;
    params.gamma = 0.2;
    params.drudeStrength = 1.0;

    try {
        FDTD1D sim(params);
        sim.run();
        sim.writeImpulsePlasmaCSV("ImpulsePlasma.cvs");
        std::cout << "Simulation finished.\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

