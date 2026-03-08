#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "SimulationParameters.h"
#include "FDTD1D.h"


#include <iostream>
#include <exception>

inline void normalizeParamsOnPlasmaWavelength(SimulationParameters& p) {
    if (!p.useDrude) {
        throw std::runtime_error("normalizeParamsOnPlasmaWavelength: useDrude=false, no omega_p.");
    }
    if (p.оmega_p <= 0.0) {
        throw std::runtime_error("normalizeParamsOnPlasmaWavelength: drudeOmegaP must be > 0.");
    }

    const double L = 2.0 * M_PI / p.оmega_p; // λp в текущих единицах

    p.dx /= L;
    p.dt /= L; //по-другому надо сделать
    p.sourceFreq   *= L;
    p.sourceFWidth *= L;
    p.оmega_p  *= L;         // станет 2π
    p.gamma   *= L;         // Γ' = Γ/ωp * 2π


    // p.dt = p.courantNumber * p.dx;
}


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
    params.sourceFWidth = (1.0 / wavelength_min) - (1.0 / wavelength_max); // cycles/time
    params.source_pos = 50;

    // Drude ADE
    params.useDrude = true;
    params.plasmaStart = 100;
    params.plasmaEnd = 200;

    params.оmega_p = 10.0;
    params.gamma = 0.5;
    params.drudeStrength = 1.0;

    try {
        normalizeParamsOnPlasmaWavelength(params);

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

