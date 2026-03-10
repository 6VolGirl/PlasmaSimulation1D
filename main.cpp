#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "SimulationParameters.h"
#include "FDTD1D.h"
#include "SpectrumAnalyzer.h"


#include <iostream>
#include <exception>


void normalizeParamsOnPlasmaWavelength(SimulationParameters& p) {
    if (!p.useDrude) {
        throw std::runtime_error("normalizeParamsOnPlasmaWavelength: useDrude=false, no omega_p.");
    }
    if (p.оmega_p <= 0.0) {
        throw std::runtime_error("normalizeParamsOnPlasmaWavelength: drudeOmegaP must be > 0.");
    }

    const double L = 1.0 / p.оmega_p; // λp в текущих единицах

    p.dx /= L;
    p.dt /= L; //по-другому надо сделать?
    p.sourceFreq   *= L;
    p.sourceFWidth *= L;
    p.оmega_p  *= L;         // станет 2π
    p.gamma   *= L;         // Γ' = Γ/ωp * 2π

    // p.dt = p.courantNumber * p.dx;
}

int main() {
    // SimulationParameters params;
    //
    // params.resolution = 20;
    //
    // // Сетка и шаги
    // params.nx = 400;
    // params.dx = 1.0 / params.resolution;
    // params.courantNumber = 0.5;
    // params.dt = params.courantNumber * params.dx;
    // params.numTimeSteps = 1000;
    //
    // // Материал
    // params.epsInf = 1.0;
    // params.mu0 = 1.0;
    //
    // // PML
    // params.pmlThickness = 20;
    // params.pmlDamping = 1e-9;
    // params.pmlProfilePower = 3;
    //
    // // Источник
    // double wavelength_min = 0.3;
    // double wavelength_max = 0.8;
    // params.sourceFreq = 2.0 / (wavelength_min + wavelength_max);     // cycles/time
    // params.sourceFWidth = 1.0 / wavelength_min - 1.0 / wavelength_max; // cycles/time
    // params.source_pos = 50;
    //
    // // Drude ADE
    // params.useDrude = true;
    // params.plasmaStart = 100;
    // params.plasmaEnd = 200;
    //
    // params.оmega_p = 8.0;
    // params.gamma = 0.2;
    // params.drudeStrength = 1.0;

    SimulationParameters params;

    params.оmega_p = 8.0;
    const double L = 1.0 / params.оmega_p;
    params.plasmaWidth = 4.0 * L;
    params.sourceFreq = 0.5 / L;


    // Сетка и шаги
    params.nx = 400;
    params.dx = 0.06 * L;
    params.courantNumber = 0.5;
    params.dt = params.courantNumber * params.dx;
    params.numTimeSteps = 1000;

    params.resolution = static_cast<int>(std::lround(1.0 / params.dx));

    // Материал
    params.epsInf = 1.0;
    params.mu0 = 1.0;

    // PML
    params.pmlThickness = 20;
    params.pmlDamping = 1e-9;
    params.pmlProfilePower = 3;

    // Источник
    params.sourceFWidth = 0.5 / L;
    params.source_pos = static_cast<int>(std::lround(3.0 * L / params.dx));

    // Drude ADE
    params.useDrude = true;
    params.plasmaStart = static_cast<int>(std::lround(6.0 * L / params.dx));
    params.plasmaEnd = static_cast<int>(std::lround((6.0 * L + params.plasmaWidth) / params.dx));

    params.gamma = 0.2 / L;
    params.drudeStrength = 2.5;

    try {
        normalizeParamsOnPlasmaWavelength(params);

        SimulationParameters paramsVac = params;
        SimulationParameters paramsPlasma = params;

        paramsVac.useDrude = false;
        paramsPlasma.useDrude = true;

        FDTD1D simVac(paramsVac);
        simVac.addMonitor(params.plasmaStart);
        simVac.addMonitor(params.plasmaEnd);
        simVac.run();
        simVac.writeImpulsePlasmaCSV("ImpulseVac.cvs");
        simVac.writeAllMonitorsCSV("monitorsVac.cvs");


        FDTD1D simPlasma(paramsPlasma);
        simPlasma.addMonitor(params.plasmaStart);
        simPlasma.addMonitor(params.plasmaEnd);
        simPlasma.run();
        simPlasma.writeImpulsePlasmaCSV("ImpulsePlasma.cvs");
        simPlasma.writeAllMonitorsCSV("monitors.cvs");


        double t0Vac = simVac.getMonitor(0).centroidTime();
        double t1Vac = simVac.getMonitor(1).centroidTime();
        double t0Plasma = simPlasma.getMonitor(0).centroidTime();
        double t1Plasma = simPlasma.getMonitor(1).centroidTime();

        double tauTunPlasma = t1Plasma - t0Plasma;
        std::cout << "Total transit time through the plasma layer = " << tauTunPlasma << "\n";

        double tauVac = t1Vac - t0Vac;
        std::cout << "Total transit time through the vacuum layer = " << tauVac << "\n";

        double deltaTau = tauTunPlasma - tauVac;
        std::cout << "Delay due to plasma = " << deltaTau << "\n";


        int nSteps = 0.0;
        if (simVac.getMonitor(0).getEx().size() == simPlasma.getMonitor(1).getEx().size())
            nSteps = simPlasma.getMonitor(0).getEx().size();
        double rt_coeff = 0.0;
        double sumInc = 0.0;
        double sumPast = 0.0;

        for ( int k = 0; k < nSteps; k++) {
            double falling = simVac.getMonitor(0).getEx()[k];
            double pastPart = simPlasma.getMonitor(1).getEx()[k];

            sumInc += falling*falling;
            sumPast += pastPart* pastPart;
        }
        if (sumPast>1e-15) {
            rt_coeff = sumPast/sumInc;
        }
        else {
            rt_coeff = 0.0;
        }
        std::cout << "rt_coeff = " << rt_coeff << "\n";

        std::cout << "Simulation finished.\n";



        SpectrumAnalyzer specIn;
        specIn.buildFluxSpectrum(simPlasma.getMonitor(0));
        specIn.writeSpectrumCSV("spectrum_in.csv", params.sourceFreq);

        SpectrumAnalyzer specOut;
        specOut.buildFluxSpectrum(simPlasma.getMonitor(1));
        specOut.writeSpectrumCSV("spectrum_out.csv", params.sourceFreq);

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
};

