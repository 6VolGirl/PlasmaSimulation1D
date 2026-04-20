#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include "SimulationParameters.h"
#include "FDTD1D.h"
#include "SpectrumAnalyzer.h"


#include <iostream>
#include <exception>

struct ScanResult {
    double freq;
    double fwidth;
    double plasmaWidth;
    double tauVac;
    double tauTunPlasma;
    double deltaTau;
    double T_coeff;
};


void normalizeParamsOnPlasmaWavelength(SimulationParameters& p) {
    if (!p.useDrude) {
        throw std::runtime_error("normalizeParamsOnPlasmaWavelength: useDrude=false, no omega_p.");
    }
    if (p.оmega_p <= 0.0) {
        throw std::runtime_error("normalizeParamsOnPlasmaWavelength: drudeOmegaP must be > 0.");
    }

    const double L = 2 * M_PI / p.оmega_p; // λp в текущих единицах

    p.dx /= L;
    p.dt /= L;
    p.sourceFreq   *= L;
    p.sourceFWidth *= L;
    p.оmega_p  *= L;         // станет 2π
    p.gamma   *= L;         // Γ' = Γ/ωp * 2π

    // p.dt = p.courantNumber * p.dx;
}

ScanResult runOnePoint(const SimulationParameters& baseParams) {
    SimulationParameters params = baseParams;
    normalizeParamsOnPlasmaWavelength(params);

    // Вакуум
    SimulationParameters paramsVac = params;
    paramsVac.useDrude = false;
    FDTD1D simVac(paramsVac);
    simVac.addMonitor(params.plasmaStart-5);
    simVac.addMonitor(params.plasmaEnd+5);
    simVac.run();

    // Плазма
    SimulationParameters paramsPlasma = params;
    paramsPlasma.useDrude = true;
    FDTD1D simPlasma(paramsPlasma);
    simPlasma.addMonitor(params.plasmaStart-5);
    simPlasma.addMonitor(params.plasmaEnd+5);
    simPlasma.run();

    // Времена
    double t0Vac    = simVac.getMonitor(0).centroidTime();
    double t1Vac    = simVac.getMonitor(1).centroidTime();
    double t1Plasma = simPlasma.getMonitor(1).centroidTime();

    double tauVac       = t1Vac - t0Vac;
    double tauTunPlasma = t1Plasma - t0Vac;
    double deltaTau     = tauTunPlasma - tauVac;

    // Прохождение
    double fluxIn  = simVac.getMonitor(0).integratedFlux();
    double fluxOut = simPlasma.getMonitor(1).integratedFlux();
    double T_coeff = (std::abs(fluxIn) > 1e-30) ? fluxOut / fluxIn : 0.0;

    return {
        baseParams.sourceFreq,
        baseParams.sourceFWidth,
        baseParams.plasmaWidth,
        tauVac,
        tauTunPlasma,
        deltaTau,
        T_coeff
    };
}




int main() {
    SimulationParameters params;

    params.оmega_p = 1.0;
    const double L = 2 * M_PI / params.оmega_p;
    params.plasmaWidth = 4.0 * L;
    params.sourceFreq = 0.60 / L;
    params.sourceFWidth = 1.0  / L;


    // Сетка и шаги
    params.nx = 1400;
    params.dx = 0.06 * L;
    params.courantNumber = 0.5;
    params.dt = params.courantNumber * params.dx;
    params.numTimeSteps = 3000;

    params.resolution = static_cast<int>(std::lround(1.0 / params.dx));

    // Материал
    params.epsInf = 1.0;
    params.mu0 = 1.0;

    // PML
    params.pmlThickness = 20;
    params.pmlDamping = 1e-9;
    params.pmlProfilePower = 3;

    // Источник
    params.source_pos = static_cast<int>(std::lround(3.0 * L / params.dx));

    // Drude ADE
    params.useDrude = true;
    params.plasmaStart = static_cast<int>(std::lround(6.0 * L / params.dx));
    params.plasmaEnd = static_cast<int>(std::lround((6.0 * L + params.plasmaWidth) / params.dx));

    params.gamma = 0.2 / L;
    params.drudeStrength = 1.0;

    params.chirpRate = 1.0;




    SimulationParameters baseParams = params;

    try {
        normalizeParamsOnPlasmaWavelength(params);

        // std::vector<double> fwidths = {1.0 / L, 1.5 / L, 2.0 / L};
        // std::vector<double> widths1 = {4.0 * L};
        // const int Nfreq = 20;
        // const double f_min = 0.2 / L;
        // const double f_max = 2.0 / L;
        //
        // // Анализ времени тунеллирования от частоты
        // std::ofstream csv1("tau_vs_freq.csv");
        // csv1 << std::scientific << std::setprecision(8);
        // csv1 << "f_over_omega_p,sourceFWidth,plasmaWidth_over_L,tauVac,tauTunPlasma,deltaTau,T_coeff\n";
        //
        // for (double pw : widths1) {
        //     for (double fw : fwidths) {
        //         for (int k = 0; k < Nfreq; ++k) {
        //             double f = f_min + (f_max - f_min) * k / (Nfreq - 1);
        //
        //             SimulationParameters p = baseParams;
        //
        //             p.sourceFreq   = f;
        //             p.sourceFWidth = fw;
        //             p.plasmaWidth  = pw;
        //
        //             p.plasmaStart = static_cast<int>(std::lround(8.0 * L / p.dx));
        //             p.plasmaEnd   = static_cast<int>(std::lround((8.0 * L + pw) / p.dx));
        //
        //             ScanResult r = runOnePoint(p);
        //
        //             double f_over_wp = f * L;
        //
        //             csv1 << f_over_wp << "," << (fw * L) << "," << pw / L << "," << r.tauVac << ","
        //                 << r.tauTunPlasma << "," << r.deltaTau << "," << r.T_coeff << "\n";
        //
        //             std::cout << "f/wp = " << f_over_wp << "  fw = " << (fw * L) << "  deltaTau = " << r.deltaTau
        //                       << "  T = " << r.T_coeff << "\n";
        //         }
        //     }
        // }
        //
        // csv1.close();
        // std::cout << "Scan written to tau_vs_freq.csv\n";
        //



        // std::vector<double> widths2 = {25.0 * L, 20.0 * L, 10.0 * L, 9.0 * L, 8.5 * L, 8.0 * L, 6.0 * L, 4.0 * L, 3.5 * L, 3.0 * L, 2.5 * L, 2.0 * L, 1.8 * L, 1.5 * L, 1.3 * L, 1.0 * L, 0.9 * L, 0.7 * L, 0.5 * L};
        //
        // // Анализ времени туннелирования от ширины плазмы
        // std::ofstream csv2("tau_vs_width.csv");
        // csv2 << std::scientific << std::setprecision(8);
        // csv2 << "PlasmaWidth,tauTunVacuum,tauTunPlasma,deltaTau,Tcoeff\n";
        //
        //
        // double fixedFreq = params.sourceFreq;
        // double fixedFWidth = params.sourceFWidth;
        //
        // for (double pw : widths2) {
        //     SimulationParameters p = baseParams;
        //
        //     p.sourceFreq   = fixedFreq;
        //     p.sourceFWidth = fixedFWidth;
        //     p.plasmaWidth  = pw;
        //
        //     p.plasmaStart = static_cast<int>(std::lround(6.0 * L / p.dx));
        //     p.plasmaEnd   = static_cast<int>(std::lround((6.0 * L + pw) / p.dx));
        //
        //     ScanResult r = runOnePoint(p);
        //
        //     csv2 << (pw / L) << ","<< r.tauVac << ","<< r.tauTunPlasma << ","
        //         << r.deltaTau << ","<< r.T_coeff << "\n";
        //
        //     std::cout << "PlasmaWidth = " << (pw / L) << "  tauVac = " << r.tauVac << "  tauPlasma = " << r.tauTunPlasma
        //               << "  deltaTau = " << r.deltaTau << "  T = " << r.T_coeff << "\n";
        // }
        //
        // csv2.close();
        // std::cout << "Scan written to tau_vs_width.csv\n";



        // std::ofstream csvChirp("tau_vs_chirp.csv");
        // csvChirp << std::scientific << std::setprecision(8);
        // csvChirp << "chirpRateNorm,tauTunVacuum,tauTunPlasma,deltaTau,Tcoeff\n";
        //
        // std::vector<double> chirpsNorm = {-10.0, -5.0, -2.0, -1.0, -0.5, -0.2, 0.0, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0};
        //
        // double freq = 0.60 / L;
        // double fWidth = 1.0 / L;
        // double width = 4.0 * L;
        //
        // for (double Cnorm : chirpsNorm) {
        //     SimulationParameters p = baseParams;
        //
        //     p.sourceFreq   = freq;
        //     p.sourceFWidth = fWidth;
        //     p.plasmaWidth  = width;
        //
        //     p.plasmaStart = static_cast<int>(std::lround(6.0 * L / p.dx));
        //     p.plasmaEnd   = static_cast<int>(std::lround((6.0 * L + p.plasmaWidth) / p.dx));
        //
        //     p.chirpRate  = Cnorm;
        //     ScanResult r = runOnePoint(p);
        //
        //     csvChirp << p.chirpRate << ","
        //              << r.tauVac << ","
        //              << r.tauTunPlasma << ","
        //              << r.deltaTau << ","
        //              << r.T_coeff << "\n";
        //
        //     std::cout << "C = " << Cnorm
        //               << "  tauVac = " << r.tauVac
        //               << "  tauPlasma = " << r.tauTunPlasma
        //               << "  deltaTau = " << r.deltaTau
        //               << "  T = " << r.T_coeff << "\n";
        // }
        //
        // csvChirp.close();
        // std::cout << "Scan written to tau_vs_chirp.csv\n";


        //τ_p = sqrt(4ln2) / fwidth
        //std::vector<double> fwidths = { 1.665 / L, 0.832 / L, 0.333 / L, 0.0832 / L};
         std::vector<double> fwidths = {1.0 / L, 0.5 / L, 0.2 / L, 0.05 / L};

        std::vector<double> widths2 = {
            0.5 * L, 0.6 * L, 0.7 * L, 0.8 * L, 0.9 * L,
            1.0 * L, 1.1 * L, 1.2 * L, 1.3 * L, 1.4 * L,
            1.5 * L, 1.6 * L, 1.7 * L, 1.8 * L, 1.9 * L, 2.0 * L,};

        // фиксируем частоту, как в твоей формулировке: f / fp = 1/4
        double fixedFreq = 0.25 / L;

        // файл для графика tau(d/L) при разных длительностях импульса
        std::ofstream csv_width_multi("tau_vs_width_multi.csv");
        csv_width_multi << std::scientific << std::setprecision(8);
        csv_width_multi << "PlasmaWidth_over_L,"
                          "sourceFWidth_over_1_over_L,"
                          "tauPulse_times_fL,"
                          "tauLight_times_fL,"
                          "tauVac,"
                          "tauTunPlasma,"
                          "deltaTau,"
                          "Tcoeff\n";

        for (double fw : fwidths) {
            for (double pw : widths2) {
                SimulationParameters p = baseParams;

                p.sourceFreq   = fixedFreq;
                p.sourceFWidth = fw;
                p.plasmaWidth  = pw;

                p.plasmaStart = static_cast<int>(std::lround(6.0 * L / p.dx));
                p.plasmaEnd   = static_cast<int>(std::lround((6.0 * L + pw) / p.dx));

                ScanResult r = runOnePoint(p);

                const double d_over_L = pw / L;
                const double fw_norm  = fw * L;

                const double tauPulse_times_fL = 1.0 / fw_norm;

                const double tauLight_times_fL = d_over_L;

                csv_width_multi
                    << d_over_L << ","
                    << fw_norm << ","
                    << tauPulse_times_fL << ","
                    << tauLight_times_fL << ","
                    << r.tauVac << ","
                    << r.tauTunPlasma << ","
                    << r.deltaTau << ","
                    << r.T_coeff << "\n";

                std::cout
                    << "d/L = " << d_over_L
                    << "  fw*L = " << fw_norm
                    << "  tauPulse*fL = " << tauPulse_times_fL
                    << "  tauPlasma = " << r.tauTunPlasma
                    << "  deltaTau = " << r.deltaTau
                    << "  T = " << r.T_coeff
                    << "\n";
            }
        }

        csv_width_multi.close();
        std::cout << "Scan written to tau_vs_width_multi.csv\n";





        //normalizeParamsOnPlasmaWavelength(params);

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


        double t0Vac    = simVac.getMonitor(0).centroidTimeByPoynting();
        double t1Vac    = simVac.getMonitor(1).centroidTimeByPoynting();
        double t1Plasma = simPlasma.getMonitor(1).centroidTimeByPoynting();

        double tauVac = t1Vac - t0Vac;
        std::cout << "Vacuum transit time between monitors = " << tauVac << "\n";


        double tauTunPlasma = t1Plasma - t0Vac;
        std::cout << "Plasma transit time (from vac input to plasma output) = " << tauTunPlasma << "\n";

        double deltaTau = tauTunPlasma - tauVac;
        std::cout << "Tunneling delay (plasma - vacuum) = " << deltaTau << "\n";


        double fluxInVac     = simVac.getMonitor(0).integratedFlux();
        double fluxOutPlasma = simPlasma.getMonitor(1).integratedFlux();

        double T_coeff = 0.0;
        if (std::abs(fluxInVac) > 1e-30) {
            T_coeff = fluxOutPlasma / fluxInVac;
        }

        std::cout << "Transmission coefficient T = " << T_coeff << "\n";
        std::cout << "|T|^2 = " << T_coeff * T_coeff << "\n";

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

