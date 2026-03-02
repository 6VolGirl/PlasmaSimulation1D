//
// Created by 6anna on 15.02.2026.
//


#include <fstream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <complex>
#include "PMLCoefficients.h"
#include "SimulationParameters.h"
#include "DrudeADE.h"
#include "Sources.h"


class FDTD1D_PythonStyle {
public:
    explicit FDTD1D_PythonStyle(const SimulationParameters& p)
        : p_(p),
          src_(p.sourceFreq, p.sourceFWidth),
          Ex_(p.nx + 1, 0.0),
          Ex_nm1_(p.nx + 1, 0.0),
          Hy_(p.nx, 0.0),
          eps_(p.nx + 1, p.eps0),
          mu_(p.nx, p.mu0),
          pml_(p.nx, p.pmlThickness, p.eps0, p.mu0,
               p.pmlDamping, p.pmlProfilePower, p.dx),
          ade_(p)
    {
        sigmaE_ = pml_.sigmaE;
        sigmaM_ = pml_.sigmaM;

        if (p_.useDrude) {
            int i0 = std::max(0, p_.plasmaStart);
            int i1 = std::min(p_.nx, p_.plasmaStart + p_.plasmaWidth);
            for (int i = i0; i <= i1; ++i) eps_[i] = p_.epsilon_inf;
        }
    }

    void run() {
        std::vector<double> ca(p_.nx + 1), cb(p_.nx + 1), db(p_.nx);

        // коэффициенты для НЕдисперсионных узлов
        for (int i = 0; i < p_.nx + 1; ++i) {
            const double sigma_total = sigmaE_[i]; // PML conductivity
            const double denom = (2.0 * eps_[i] + sigma_total * p_.dt);
            ca[i] = (2.0 * eps_[i] - sigma_total * p_.dt) / denom;
            cb[i] = 2.0 * p_.dt / (p_.dx * denom);
        }
        for (int i = 0; i < p_.nx; ++i) {
            const double denom = (2.0 * mu_[i] + sigmaM_[i] * p_.dt);
            db[i] = 2.0 * p_.dt / (p_.dx * denom);
        }

        snapshotsEx_.clear();

        for (int n = 0; n < p_.numTimeSteps; ++n) {

            // E update: вне плазмы по python, в плазме по ADE
            for (int i = 1; i < p_.nx; ++i) {
                const double curl_h = (Hy_[i - 1] - Hy_[i]);

                bool inPlasma = p_.useDrude &&
                                (i >= p_.plasmaStart) &&
                                (i <= p_.plasmaStart + p_.plasmaWidth);

                if (!inPlasma) {
                    Ex_[i] = ca[i] * Ex_[i] + cb[i] * curl_h;
                } else {
                    const double sigma_total = sigmaE_[i] + p_.sigmaCond;
                    const double eps_inf = p_.epsilon_inf;

                    const double denom = (2.0 * eps_inf + sigma_total * p_.dt);
                    const double C2 = (2.0 * eps_inf - sigma_total * p_.dt) / denom;
                    const double C3 = (2.0 * p_.dt) / (p_.dx * denom);

                    const double pol = ade_.polCorrection(i);
                    const double Ex_new = C2 * Ex_[i] + C3 * (curl_h - pol);

                    Ex_[i] = Ex_new;
                }
            }

            // Источник
            Ex_[p_.source_pos] += src_(n * p_.dt);

            if (p_.useDrude) {
                ade_.updateJ(Ex_, Ex_nm1_);
            }

            for (int i = 0; i < p_.nx; ++i) {
                Hy_[i] = Hy_[i] + db[i] * (Ex_[i] - Ex_[i + 1]);
            }

            Ex_nm1_ = Ex_;

            if (n % 2 == 0) snapshotsEx_.push_back(Ex_);
        }
    }

    void writeImpulsePlasmaCSV(const std::string& filename) const {
        std::ofstream out(filename);
        if (!out.is_open()) {
            std::cerr << "Cannot open " << filename << " for writing\n";
            return;
        }

        out << std::scientific << std::setprecision(9);

        const size_t Nt = snapshotsEx_.size();
        const int Nx = p_.nx;
        const double invRes = 1.0 / static_cast<double>(p_.resolution);

        out << "time_over_fL,x_tilde,Ez\n";
        for (size_t k = 0; k < Nt; ++k) {
            const double t_over_fL = (2.0 * k) * p_.dt;
            for (int i = 0; i < Nx; ++i) {
                const double x_tilde = (static_cast<double>(i) - p_.source_pos) * invRes;
                const double Ez = snapshotsEx_[k][i];
                out << t_over_fL << "," << x_tilde << "," << Ez << "\n";
            }
        }
    }

private:
    const SimulationParameters& p_;
    GaussianSource src_;

    std::vector<double> Ex_, Ex_nm1_, Hy_;
    std::vector<double> eps_, mu_;
    std::vector<double> sigmaE_, sigmaM_;
    PMLSigma pml_;

    DrudeADE ade_;

    std::vector<std::vector<double>> snapshotsEx_;
};
