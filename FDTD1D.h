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
          Ex_nm1_(p.nx + 1, 0.0),  // E^{n-1}
          Ex_n_(p.nx + 1, 0.0),    // E^{n}
          Ex_np1_(p.nx + 1, 0.0),  // E^{n+1}
          Hy_(p.nx, 0.0),          // H^{n+1/2} (храним полушаг)
          epsInf_(p.nx + 1, p.epsInf0),
          mu_(p.nx, p.mu0),
          pml_(p.nx, p.pmlThickness, p.epsInf0, p.mu0,
               p.pmlDamping, p.pmlProfilePower, p.dx),
          drude_(p.nx + 1)
    {
        sigmaE_ = pml_.sigmaE;
        sigmaM_ = pml_.sigmaM;

        if (p_.useDrude && p_.drudeEnd >= p_.drudeStart) {
            drude_.configureRegion(p_.drudeStart, p_.drudeEnd,
                                   p_.drudeOmegaP, p_.drudeGamma,
                                   p_.drudeStrength, p_.dt);
        }
    }

    void run() {
        // Магнитные коэффициенты (как в python PML-версии: da, db) [file:41]
        std::vector<double> da(p_.nx), db(p_.nx);
        for (int i = 0; i < p_.nx; ++i) {
            const double denom = (2.0 * mu_[i] + sigmaM_[i] * p_.dt);
            da[i] = (2.0 * mu_[i] - sigmaM_[i] * p_.dt) / denom;
            db[i] = (2.0 * p_.dt) / (p_.dx * denom);
        }

        snapshotsEx_.clear();


        for (int n = 0; n < p_.numTimeSteps; ++n) {
            // 1) E^{n+1} из H^{n+1/2} и (E^n, E^{n-1}, J^n, J^{n-1})  [file:3]
            Ex_np1_ = Ex_n_; // можно как старт, потом перезапишем внутренние узлы

            for (int i = 1; i < p_.nx; ++i) {
                const double gsum = drude_.gammaSum(i); // Σγj [file:3]
                const double denom = (2.0 * epsInf_[i] + sigmaE_[i] * p_.dt + 0.5 * gsum); // [file:3]

                const double C1 = (0.5 * gsum) / denom;                                    // [file:3]
                const double C2 = (2.0 * epsInf_[i] - sigmaE_[i] * p_.dt) / denom;         // [file:3]
                const double C3 = (2.0 * p_.dt) / denom;                                   // [file:3]

                const double curlH = (Hy_[i - 1] - Hy_[i]) / p_.dx;                        // ∇×H 1D
                const double Jterm = drude_.Jhalf_noEterm(i);                              // [file:3]

                Ex_np1_[i] = C1 * Ex_nm1_[i] + C2 * Ex_n_[i] + C3 * (curlH - Jterm);       // [file:3]
            }

            // Границы (можно 0; при PML обычно и так ок)
            Ex_np1_[0] = 0.0;
            Ex_np1_[p_.nx] = 0.0;

            // Источник (как в python: добавка в E после апдейта E) [file:41]
            const double t = n * p_.dt;
            Ex_np1_[p_.source_pos] += src_(t);

            // 2) J^{n+1} (Drude/Lorentz ADE) [file:3]
            if (p_.useDrude) {
                drude_.updateJ(Ex_np1_, Ex_nm1_, p_.dt);
            }

            // 3) H^{n+3/2} из H^{n+1/2} и E^{n+1} (da/db как в python PML) [file:41]
            for (int i = 0; i < p_.nx; ++i) {
                Hy_[i] = da[i] * Hy_[i] + db[i] * (Ex_np1_[i] - Ex_np1_[i + 1]);
            }

            // Сдвиг времен: (n-1,n,n+1) -> (n,n+1,...) [file:3]
            Ex_nm1_.swap(Ex_n_);
            Ex_n_.swap(Ex_np1_);

            if (n % 2 == 0) snapshotsEx_.push_back(Ex_n_);
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

                out << t_over_fL << ","
                    << x_tilde   << ","
                    << Ez        << "\n";
            }
        }

        out.close();
        std::cout << "Impulse snapshots written to " << filename << " (CSV)\n";
    }

    // остальное без изменений...

private:
    const SimulationParameters& p_;
    GaussianSource src_;

    std::vector<double> Ex_nm1_, Ex_n_, Ex_np1_;
    std::vector<double> Hy_;
    std::vector<double> epsInf_, mu_;
    std::vector<double> sigmaE_, sigmaM_;
    PMLSigma pml_;
    DrudeADE drude_;

    std::vector<std::vector<double>> snapshotsEx_;
};
