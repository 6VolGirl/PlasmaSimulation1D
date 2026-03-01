#include <iostream>
#include "SimulationParameters.h"
#include "FDTD1D.h"
#include "DrudeADE.h"

int main() {
    SimulationParameters params;

    params.a = 800e-9;
    params.resolution = 100;

    // Grid and time
    params.nx = 256;
    params.ny = 1;
    params.nz = 1;
    params.numTimeSteps = 20000;
    params.dx = 10e-9;
    params.dy = params.dz = params.dx;
    params.courantNumber = 0.5;

    // Plasma
    params.omega_p = 0.0; //1.4e16;
    params.gamma = 0.0; //1.0e14;
    params.epsilon_inf = 1.0;
    params.sigmaCond = 0.0;     // Drude


    //Тестим
    // const double c0 = 299792458.0;
    // double fL = c0 / params.a;          // Hz
    // double fp = 0.9 * fL;
    // params.sourceFreqCenter = 2.0 * M_PI * fp;
    // params.pulseWidth = 1.5 / fL;
    // Source
    params.source_pos = 40;
    params.sourceFreqCenter = 1e18;
    params.sourceFreqWidth = 0.5e16;
    params.pulseWidth = 20e-15;

    // PML
    params.pmlThickness = 20;
    params.pmlReflectCoeff = 1e-6;
    params.pmlGrading_m = 3.0;

    // Monitors
    params.monitorFront = 50;
    params.plasmaStart = 50;
    params.plasmaWidth = 50;
    params.monitorBack = 101;


    try {
        FDTD1D sim(params);

        std::vector<double> times_over_fL = {2.0, 10.0, 50.0, 100.0, 1000.0};
        sim.setSnapshotTimes_fl(times_over_fL);

        sim.run();
        sim.analyzeFourierSpectra("reflection_transmission1.txt");
        sim.writeImpulsePlasma("ImpulsePlasma.txt", times_over_fL);

        std::cout << "\nSimulation successful!\n";
        std::cout << "Output files:\n"
                  << "  - reflection_transmission.txt (R/T coefficients)\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
