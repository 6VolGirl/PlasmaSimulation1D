
#include <fstream>
#include <iomanip>
#include <iostream>
#include "FDTD1D.h"

FDTD1D::FDTD1D(const SimulationParameters& p)
    : p_(p),
      src_(p.sourceFreq, p.sourceFWidth),
      Ex_nm1_(p.nx + 1, 0.0),  // E^{n-1}
      Ex_n_(p.nx + 1, 0.0),    // E^{n}
      Ex_np1_(p.nx + 1, 0.0),  // E^{n+1}
      Hy_(p.nx, 0.0),          // H^{n+1/2} (храним полушаг)
      epsInf_(p.nx + 1, p.epsInf),
      mu_(p.nx, p.mu0),
      pml_(p.nx, p.pmlThickness, p.epsInf, p.mu0,
           p.pmlDamping, p.pmlProfilePower, p.dx),
      drude_(p.nx + 1)
{
    sigmaE_ = pml_.sigmaE;
    sigmaM_ = pml_.sigmaM;

    if (p_.useDrude && p_.plasmaEnd >= p_.plasmaStart) {
        drude_.configureRegion(p_.plasmaStart, p_.plasmaEnd,
                               p_.оmega_p, p_.gamma,
                               p_.drudeStrength, p_.dt);
    }
}

void FDTD1D::run() {
        std::vector<double> da(p_.nx), db(p_.nx);
        for (int i = 0; i < p_.nx; ++i) {
            const double denom = (2.0 * mu_[i] + sigmaM_[i] * p_.dt);
            da[i] = (2.0 * mu_[i] - sigmaM_[i] * p_.dt) / denom;
            db[i] = (2.0 * p_.dt) / (p_.dx * denom);
        }

        snapshotsEx_.clear();


        for (int n = 0; n < p_.numTimeSteps; ++n) {
            Ex_np1_ = Ex_n_;

            for (int i = 1; i < p_.nx; ++i) {
                const double gsum = drude_.gammaSum(i);
                const double denom = (2.0 * epsInf_[i] + sigmaE_[i] * p_.dt + 0.5 * gsum);

                const double C1 = (0.5 * gsum) / denom;
                const double C2 = (2.0 * epsInf_[i] - sigmaE_[i] * p_.dt) / denom;
                const double C3 = (2.0 * p_.dt) / denom;

                const double curlH = (Hy_[i - 1] - Hy_[i]) / p_.dx;
                const double Jterm = drude_.Jhalf_noEterm(i);

                Ex_np1_[i] = C1 * Ex_nm1_[i] + C2 * Ex_n_[i] + C3 * (curlH - Jterm);
            }

            // Границы
            Ex_np1_[0] = 0.0;
            Ex_np1_[p_.nx] = 0.0;

            // Источник
            const double t = n * p_.dt;
            Ex_np1_[p_.source_pos] += src_(t);

            // 2) J^{n+1} (Drude/Lorentz ADE)
            if (p_.useDrude) {
                drude_.updateJ(Ex_np1_, Ex_nm1_, p_.dt);
            }

            // 3) H^{n+3/2} из H^{n+1/2} и E^{n+1}
            for (int i = 0; i < p_.nx; ++i) {
                Hy_[i] = da[i] * Hy_[i] + db[i] * (Ex_np1_[i] - Ex_np1_[i + 1]);
            }

            Ex_nm1_.swap(Ex_n_);
            Ex_n_.swap(Ex_np1_);

            sampleMonitors(t);

            if (n % 2 == 0) snapshotsEx_.push_back(Ex_n_);
        }
    }

void FDTD1D::writeImpulsePlasmaCSV(const std::string& filename) const {
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


void FDTD1D::addMonitor(int pos) {
    if (pos < p_.pmlThickness || pos >= static_cast<int>(Ex_n_.size())-p_.pmlThickness) {
        throw std::out_of_range("addMonitor: position is out of Ex range or in PML");
    }

    Monitor m;
    m.position = pos;
    m.reserve(p_.numTimeSteps);

    monitors_.push_back(m);

}


void FDTD1D::sampleMonitors(double t) {
    for (auto& m : monitors_) {
        m.sample(t, Ex_n_[m.position]);
    }
}


double FDTD1D::tunnelingTime(std::size_t inMonitor, std::size_t outMonitor) const {
    if (inMonitor >= monitors_.size() || outMonitor >= monitors_.size()) {
        throw std::out_of_range("tunnelingTime: monitor index is out of range");
    }

    const Monitor& in  = monitors_[inMonitor];
    const Monitor& out = monitors_[outMonitor];

    if (in.time.empty() || out.time.empty()) {
        throw std::runtime_error("tunnelingTime: monitor data is empty");
    }

    return out.centroidTime() - in.centroidTime();
    }


void FDTD1D::writeAllMonitorsCSV(const std::string& filename) const {
    if (monitors_.empty()) {
        std::cerr << "writeAllMonitorsCSV: no monitors to write\n";
        return;
    }

    const std::size_t nSteps = monitors_[0].time.size();

    for (std::size_t i = 0; i < monitors_.size(); ++i) {
        if (monitors_[i].time.size() != nSteps ||
            monitors_[i].fieldEx.size() != nSteps ||
            monitors_[i].intensity.size() != nSteps) {
            throw std::runtime_error("writeAllMonitorsCSV: monitor arrays have different sizes");
            }
    }

    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Cannot open " << filename << " for writing\n";
        return;
    }

    out << std::scientific << std::setprecision(9);

    out << "time";
    for (const auto& m : monitors_) {
        out << ",Ex" << m.position
            << ",I" << m.position;
    }
    out << "\n";

    for (int k = 0; k < nSteps; k++) {
        out << monitors_[0].time[k];
        double rt_coeff = (monitors_[0].fieldEx[k] > 1e-15) ? monitors_[1].fieldEx[k] / monitors_[0].fieldEx[k] : 0.0;
        out << "," << monitors_[0].fieldEx[k] << ", " << monitors_[0].intensity[k]
            << "," << monitors_[1].fieldEx[k] << ", " << monitors_[1].intensity[k]
            << rt_coeff << "\n";
    }

    out.close();
    std::cout << "All monitors in file" << filename << "\n";
}
