//
// Created by 6anna on 15.02.2026.
//


#include <fstream>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "PMLCoefficients.h"
#include "SimulationParameters.h"
#include "DrudeADE.h"
#include "Sources.h"

class FDTD1D {

    const SimulationParameters& p_;
    GaussianSource src_;

    std::vector<double> Ex_prev, Ex_n_, Ex_next;
    std::vector<double> Hy_;
    std::vector<double> epsInf_, mu_;
    std::vector<double> sigmaE_, sigmaM_;
    PMLSigma pml_;
    DrudeADE drude_;

    std::vector<std::vector<double>> snapshotsEx_;

public:
    explicit FDTD1D(const SimulationParameters& p);

    void run();

    void writeImpulsePlasmaCSV(const std::string& filename) const ;


};
