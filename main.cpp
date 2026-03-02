#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "SimulationParameters.h"
#include "FDTD1D.h"

// здесь должны быть подключения заголовков/кода с SimulationParameters,
// GaussianSource, PMLSigma, FDTD1D_PythonStyle

int main() {
    SimulationParameters params;
    params.resolution= 20;

    // --- сетка и шаги ---
    params.nx             = 400;     // число ячеек H (E будет nx+1)
    params.dx             = 1.0 / params.resolution; // как в python: dx = 1/resolution
    params.courantNumber  = 0.5;
    params.dt             = params.courantNumber * params.dx;
    params.numTimeSteps   = 1000;

    // --- базовый материал ---
    params.eps0 = 1.0;
    params.mu0  = 1.0;

    // --- PML ---
    params.pmlThickness   = 20;      // в ячейках
    params.pmlDamping     = 1e-9;    // целевой уровень затухания (как python damping)
    params.pmlProfilePower= 3;       // cubic

    // --- источник (как GaussianSource в python) ---
    // пример: центральная длина волны ~ 0.5 (безразмерная),
    // тогда частота ~ 1/0.5 = 2 (cycles/единицу времени)
    double wavelength_min = 0.3;
    double wavelength_max = 0.8;
    params.sourceFreq   = 2.0 / (wavelength_min + wavelength_max);          // freq (cycles/time)
    params.sourceFWidth = 1.0 / wavelength_min - 1.0 / wavelength_max;      // bandwidth

    params.source_pos   = 50;       // ячейка источника (внутри области, не в PML)

    // --- мониторы (заготовка) ---
    params.monitorFront = 80;
    params.monitorBack  = 300;

    // --- запуск симуляции ---
    try {
        FDTD1D_PythonStyle sim(params);
        sim.run();
        std::vector<double> times_over_fL = {2.0, 10.0, 50.0, 100.0, 1000.0};
        sim.writeImpulsePlasmaCSV("ImpulsePlasma.cvs", times_over_fL);

        std::cout << "Simulation finished.\n";

        // здесь можно добавить:
        //  - сохранение Ex_/Hy_ (нужно сделать геттеры)
        //  - расчёт спектров в мониторах
        //  - вывод в файл

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
