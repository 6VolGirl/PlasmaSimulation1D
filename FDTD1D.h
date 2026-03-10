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
#include "Monitor.h"

class FDTD1D {

    const SimulationParameters& p_;
    GaussianSource src_;

    std::vector<double> Ex_nm1_, Ex_n_, Ex_np1_;
    std::vector<double> Hy_;
    std::vector<double> epsInf_, mu_;
    std::vector<double> sigmaE_, sigmaM_;
    PMLSigma pml_;
    DrudeADE drude_;

    std::vector<std::vector<double>> snapshotsEx_;
    std::vector<Monitor> monitors_;

public:
    explicit FDTD1D(const SimulationParameters& p);

    void run();

    void writeImpulsePlasmaCSV(const std::string& filename) const ;

    void addMonitor(int pos);
    void sampleMonitors(double t);
    double tunnelingTime(std::size_t inMonitor, std::size_t outMonitor) const;

    void writeAllMonitorsCSV(const std::string& filename) const;

    const Monitor& getMonitor(std::size_t index) const;
};


