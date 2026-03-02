//
// Created by 6anna on 15.02.2026.
//


#include <fstream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <complex>
#include "PMLCoefficients.h"
//#include "DrudeADE.h"
//#include "Normalization.h"
#include "SimulationParameters.h"

    struct GaussianSource {
        double freq;
        double fwidth;
        double start_time = 0.0;
        double cutoff = 3.0;

        double w, t0, finish_time;

        GaussianSource(double freq_, double fwidth_, double start=0.0, double cutoff_=3.0)
            : freq(freq_), fwidth(fwidth_), start_time(start), cutoff(cutoff_) {
            w = 1.0 / fwidth;
            t0 = start_time + cutoff * w;
            finish_time = start_time + 2.0 * cutoff * w;
        }

        double operator()(double t) const {
            if (t < start_time || t > finish_time) return 0.0;
            double tau = t - t0;
            return std::exp(-0.5 * (tau*tau) / (w*w)) * std::sin(2.0*M_PI*freq*t);
        }
    };

class FDTD1D_PythonStyle {
public:
    explicit FDTD1D_PythonStyle(const SimulationParameters& p)
        : p_(p),
          src_(p.sourceFreq, p.sourceFWidth),
          Ex_(p.nx + 1, 0.0),
          Hy_(p.nx, 0.0),
          eps_(p.nx + 1, p.eps0),
          mu_(p.nx, p.mu0),
          pml_(p.nx, p.pmlThickness, p.eps0, p.mu0,
               p.pmlDamping, p.pmlProfilePower, p.dx)
    {
        sigmaE_ = pml_.sigmaE;
        sigmaM_ = pml_.sigmaM;
    }

    void run() {
        std::vector<double> ca(p_.nx + 1), cb(p_.nx + 1), db(p_.nx);

        for (int i = 0; i < p_.nx + 1; ++i) {
            double denom = (2.0 * eps_[i] + sigmaE_[i] * p_.dt);
            ca[i] = (2.0 * eps_[i] - sigmaE_[i] * p_.dt) / denom;
            cb[i] = 2.0 * p_.dt / (p_.dx * denom);
        }
        for (int i = 0; i < p_.nx; ++i) {
            double denom = (2.0 * mu_[i] + sigmaM_[i] * p_.dt);
            db[i] = 2.0 * p_.dt / (p_.dx * denom);
        }

        snapshotsEx_.clear();

        for (int n = 0; n < p_.numTimeSteps; ++n) {
            for (int i = 1; i < p_.nx; ++i) {
                Ex_[i] = ca[i] * Ex_[i] + cb[i] * (Hy_[i - 1] - Hy_[i]);
            }

            double t = n * p_.dt;
            Ex_[p_.source_pos] += src_(t);

            for (int i = 0; i < p_.nx; ++i) {
                Hy_[i] = Hy_[i] + db[i] * (Ex_[i] - Ex_[i + 1]);
            }

            if (n % 2 == 0) {
                snapshotsEx_.push_back(Ex_);
            }
        }
    }


    void writeImpulsePlasmaCSV(const std::string& filename) const
    {
        std::ofstream out(filename);
        if (!out.is_open()) {
            std::cerr << "Cannot open " << filename << " for writing\n";
            return;
        }

        out << std::scientific << std::setprecision(9);

        const size_t Nt = snapshotsEx_.size();
        const int    Nx = p_.nx;

        const double invRes = 1.0 / static_cast<double>(p_.resolution);

        out << "time_over_fL,x_tilde,Ez\n";

        for (size_t k = 0; k < Nt; ++k) {
            double t_over_fL = (2.0 * k) * p_.dt;

            for (int i = 0; i < Nx; ++i) {
                const double x_tilde =
                    (static_cast<double>(i) - p_.source_pos) * invRes;
                const double Ez = snapshotsEx_[k][i];

                out << t_over_fL << ","
                    << x_tilde      << ","
                    << Ez           << "\n";
            }
        }

        out.close();
        std::cout << "Impulse snapshots written to " << filename << " (CSV)\n";
    }

private:
    const SimulationParameters& p_;
    GaussianSource src_;

    std::vector<double> Ex_, Hy_;
    std::vector<double> eps_, mu_;
    std::vector<double> sigmaE_, sigmaM_;
    PMLSigma pml_;

    std::vector<int> snapshotSteps_;
    std::vector<std::vector<double>> snapshotsEx_;
};


