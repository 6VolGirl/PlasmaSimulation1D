//
// Created by 6anna on 15.02.2026.
//

#include <iostream>
#include "FDTD1D.h"

extern int n_test = 0;

FDTD1D::FDTD1D(const SimulationParameters& params_si)
        : pack_(makeNormalizedParams(params_si)),
          p_(pack_.nd),
          pml_(p_.nx, p_.pmlThickness),
          ade_(p_.nx, p_, calculateTimeStep())  {

    int size = p_.nx;
    Ex_.resize(size, 0.0);
    Hy_.resize(size, 0.0);
    Ex_prev_.resize(size, 0.0);

    monitorFront_.position = p_.monitorFront;
    monitorFront_.spectrum_e.resize(p_.numTimeSteps, 0.0);

    monitorBack_.position = p_.monitorBack;
    monitorBack_.spectrum_e.resize(p_.numTimeSteps, 0.0);

    // Setup source
    const double cutoff = 3.0;
    source_.start = 0.0;
    source_.width = p_.pulseWidth;
    source_.t0 = cutoff * source_.width;
    source_.finish = cutoff * source_.width;
    source_.freqCenter = p_.sourceFreqCenter;
    source_.freqSweep = 0.0;
}

double FDTD1D::calculateTimeStep() const {
    return p_.dt;
}

void FDTD1D::run() {
    const double dt = calculateTimeStep();
    const double factor = p_.courantNumber;

    std::cout << "Starting FDTD simulation...\n"
               << "  Time steps: " << p_.numTimeSteps << "\n"
               << "  Grid size: " << p_.nx << " cells\n"
               << "  dt = " << dt << ", dx = " << p_.dx << "\n"
               << "  ωₚ = " << p_.omega_p
               << ", Γ = " << p_.gamma << "\n";

    // if (monitorFront_.position < p_.pmlThickness) {
    //     std::cerr << "  WARNING: Front monitor is inside left PML!\n";
    // }
    // if (monitorBack_.position >= p_.nx - p_.pmlThickness) {
    //     std::cerr << "  WARNING: Back monitor is inside right PML!\n";
    // }

    snapshotsEx_.clear();
    snapshotsEx_.reserve(snapshotSteps_.size());
    snap_idx_ = 0;

    // // Чистый Yee
    // for (int step = 0; step < p_.numTimeSteps; ++step) {
    //
    //     for (int i = 0; i < p_.nx - 1; ++i) {
    //         double curl_e = Ex_[i + 1] - Ex_[i];
    //         Hy_[i] += factor * curl_e;
    //     }
    //
    //     for (int i = 1; i < p_.nx; ++i) {
    //         double curl_h = Hy_[i] - Hy_[i - 1];
    //         Ex_[i] += factor * curl_h;
    //     }
    //
    //     Ex_[0]      = 0.0;
    //     Ex_[p_.nx-1]= 0.0;

    // // Чистый Yee и PML
    // for (int step = 0; step < p_.numTimeSteps; ++step) {
    //
    //     for (int i = 0; i < p_.nx - 1; ++i) {
    //         double curl_e = Ex_[i + 1] - Ex_[i];
    //
    //         double C1 = factor;
    //         double fi1 = pml_.fi1[i];
    //         double fi2 = pml_.fi2[i];
    //         double fi3 = pml_.fi3[i];
    //
    //         double AH = fi3 + fi2 * fi1;
    //         double DH = fi3 * C1;
    //
    //         Hy_[i] = AH * Hy_[i] + DH * curl_e;
    //     }
    //
    //     for (int i = 1; i < p_.nx; ++i) {
    //         double curl_h = Hy_[i] - Hy_[i - 1];
    //
    //         double C2 = 1.0;       // вакуум
    //         double C3 = factor;    // стандартный Yee
    //
    //         double gi1 = pml_.gi1[i];
    //         double gi2 = pml_.gi2[i];
    //         double gi3 = pml_.gi3[i];
    //
    //         double A = gi3 * C2 + gi2 * gi1;
    //         double D = gi3 * C3;
    //
    //         Ex_[i] = A * Ex_[i] + D * curl_h;
    //     }
    //
    //     Ex_[0]     = 0.0;
    //     Ex_[p_.nx-1]  = 0.0;

    //Полная почти рабочая реализация с РML и ADE
    const int plasmaStart = p_.plasmaStart;
    const int plasmaEnd = p_.plasmaStart + p_.plasmaWidth;

    for (int step = 0; step < p_.numTimeSteps; ++step) {

        for (int i = 1; i < p_.nx - 1; ++i) {
            double curl_e = Ex_[i + 1] - Ex_[i];
            double C1 = factor;
            double fi1 = pml_.fi1[i];
            double fi2 = pml_.fi2[i];
            double fi3 = pml_.fi3[i];

            double AH = fi3 + fi2 * fi1;
            double DH = fi3 * C1;
            double H_new = AH * Hy_[i] + DH * curl_e;
            Hy_[i] = H_new;
        }

        ade_.updatePolarizationCurrent(Ex_, Ex_prev_, dt);

        for (int i = 1; i < p_.nx - 1; ++i) {
            double curl_h = Hy_[i] - Hy_[i - 1];

            double C1, C2, C3;
            double pol_correction = 0.0;

            if (i >= plasmaStart && i < plasmaEnd) {
                double gamma_sum = ade_.getGamma(i);
                double eps_inf = p_.epsilon_inf;
                double sigma = p_.sigmaCond;
                double denominator = 2.0 * eps_inf + sigma * dt + 0.5 * gamma_sum;

                C1 = (0.5 * gamma_sum) / denominator;
                C2 = (2.0 * eps_inf - sigma * dt) / denominator;
                C3 = (2.0 * dt) / denominator;

                pol_correction = ade_.getPolarizationCorrection(i);
            } else {//вакуум
                C1 = 0.0;
                C2 = 1.0;
                C3 = factor;
            }

            double gi1 = pml_.gi1[i];
            double gi2 = pml_.gi2[i];
            double gi3 = pml_.gi3[i];

            double A = gi3 * C2 + gi2 * gi1;
            double B = gi3 * C1;
            double D = gi3 * C3;

            double Ex_new = A * Ex_[i]
                          + B * Ex_prev_[i]
                          + D * (curl_h - pol_correction);

            Ex_prev_[i] = Ex_[i];
            Ex_[i]      = Ex_new;
        }



        // Чирпированный источник
        double t = step * dt;
        Ex_[p_.source_pos] += source_.evaluateChirp(t);

        monitorFront_.spectrum_e[step] = Ex_[monitorFront_.position];
        monitorBack_.spectrum_e[step] = Ex_[monitorBack_.position];

        // if (snap_idx_ < (int)snapshotSteps_.size() && step == snapshotSteps_[snap_idx_]) {
        //     snapshotsEx_.push_back(Ex_); // копируем весь профиль
        //     ++snap_idx_;
        // }
        if (step % 2 == 0) {
            snapshotsEx_.push_back(Ex_);
        }

        if ((step + 1) % (p_.numTimeSteps / 10) == 0) {
            std::cout << "  Progress: " << (100 * (step + 1) / p_.numTimeSteps)
                     << "%, max|E| = " << getMaxField() << "\n";
        }
        //ade_.stepHistory();
    }

    std::ofstream out_t;
    out_t.open("Ex.txt");
    for (int step = 0; step < p_.numTimeSteps; ++step) {
        out_t << monitorFront_.spectrum_e[step] <<  monitorBack_.spectrum_e[step] << std::endl;
    }
    out_t.close();

}




void FDTD1D::analyzeFourierSpectra(const std::string& output_file){
    std::ofstream out(output_file, std::ios_base::out);
    out << std::scientific << std::setprecision(6);


    const int nfreq = 256;
    std::vector<double> freqs(nfreq);
    std::vector<double> spec_front(nfreq, 0.0);
    std::vector<double> spec_back(nfreq, 0.0);

    double dt = calculateTimeStep();
    double df = 1.0 / (p_.numTimeSteps * dt);

    for (int k = 0; k < nfreq; ++k) {
        freqs[k] = k * df;
        std::complex<double> sum_f(0, 0), sum_b(0, 0);

        for (int n = 0; n < p_.numTimeSteps; ++n) {
            double t = (n + 0.5) * dt; // E на полушагах времени
            double phase = -2.0 * M_PI * freqs[k] * t;
            std::complex<double> phasor(std::cos(phase), std::sin(phase));

            sum_f += monitorFront_.spectrum_e[n] * phasor * dt;
            sum_b += monitorBack_.spectrum_e[n]  * phasor * dt;

        }

        spec_front[k] = std::abs(sum_f);
        spec_back[k] = std::abs(sum_b);
    }


    out << "# Frequency  | Spectrum (Front) | Spectrum (Back) | R/T Coeff\n";
    for (int k = 0; k < nfreq; ++k) {
        double rt_coeff = (spec_front[k] > 1e-15) ? spec_back[k] / spec_front[k] : 0.0;
        out << freqs[k] << "\t" << spec_front[k] << "\t"
            << spec_back[k] << "\t" << rt_coeff << "\n";
    }
    out.close();

    std::cout << "Spectra written to " << output_file << "\n";
}

double FDTD1D::getMaxField() const {
    return *std::max_element(Ex_.begin(), Ex_.end(),
        [](double a, double b) { return std::abs(a) < std::abs(b); });
}

