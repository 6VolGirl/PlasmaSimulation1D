//
// Created by 6anna on 15.02.2026.
//

#ifndef FDTD1D_H
#define FDTD1D_H

#include <fstream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <complex>
#include "PMLCoefficients.h"
#include "DrudeADE.h"
#include "Normalization.h"


class FDTD1D {
private:
    ParamsPack pack_;
    const SimulationParameters& p_;
    PMLCoefficients pml_;
    DrudeADE ade_;

    std::vector<double> Ex_, Hy_;           // Electric and magnetic fields
    std::vector<double> Ex_prev_;           // E^{n-1} for history

    struct Monitor {
        std::vector<std::complex<double>> spectrum_e;
        int position;
    };
    Monitor monitorFront_, monitorBack_;

    struct ChirpSource {
        double start;   // starttime (безразмерное)
        double finish;  // finishtime (безразмерное)

        double t0;
        double width;
        double freqCenter, freqSweep;

        double evaluateChirp(double t) const {
            if (t < start || t > finish ) return 0.0;

            const double tau = t - t0;
            //const double envelope = std::exp(-(tau * tau) / (2.0 * width * width));
            const double envelope  =1.0;
            //const double phase = freqCenter * tau + 0.5 * freqSweep * tau * tau;
            const double phase = freqCenter * tau;
            return envelope * std::sin(phase);
        }
    } source_;


    std::vector<int> snapshotSteps_;
    std::vector<std::vector<double>> snapshotsEx_;   // snapshotsEx_[k][i] = Ex(i)
    int snap_idx_ = 0;

public:
    FDTD1D(const SimulationParameters& params);

    double calculateTimeStep() const ;
    void run();

    //Фурье анализ
    void analyzeFourierSpectra(const std::string& output_file);

    double getMaxField() const;

    const std::vector<double>& getEx() const { return Ex_; }
    const std::vector<double>& getHy() const { return Hy_; }

    void setSnapshotTimes_fl(const std::vector<double>& times_over_fL) {
        snapshotSteps_.clear();
        snapshotsEx_.clear();
        snap_idx_ = 0;

        // times_over_fL заданы как {12, 30, 100} в единицах 1/fL
        // В вашей нормировке: t_unit = a/c0, а fL = c0/a => fL * t_unit = 1
        // Значит безразмерное время t~ численно равно (t * fL).
        // Следовательно step = round(t~ / dt~) = round((t/fL) / dt~) = round(times_over_fL / dt)
        for (double tau : times_over_fL) {
            int step = (int)std::llround(tau / p_.dt);
            if (step < 0) step = 0;
            if (step >= p_.numTimeSteps) step = p_.numTimeSteps - 1;
            snapshotSteps_.push_back(step);
        }

        std::sort(snapshotSteps_.begin(), snapshotSteps_.end());
        snapshotSteps_.erase(std::unique(snapshotSteps_.begin(), snapshotSteps_.end()), snapshotSteps_.end());
    }

    void writeImpulsePlasma(const std::string& filename, const std::vector<double>& times_over_fL) const
    {
        std::ofstream out(filename);
        out << std::scientific << std::setprecision(9);

        if (snapshotsEx_.size() < times_over_fL.size()) {
            out << "# ERROR: not enough snapshots saved. Saved = " << snapshotsEx_.size()
                << ", requested = " << times_over_fL.size() << "\n";
            return;
        }

        out << "# x_tilde = (i - source_pos)/resolution  (this equals x*fL/c)\n";
        out << "# plasma region in x_tilde: [plasmaStart, plasmaStart+plasmaWidth] / resolution\n";
        //out << "# Columns: x_tilde  Ez(t=" << times_over_fL[0] << "/fL)"
        //    << "  Ez(t=" << times_over_fL[1] << "/fL)"
        //    << "  Ez(t=" << times_over_fL[2] << "/fL)\n";

        // const double invRes = 1.0 / (double)p_.resolution;
        //
        // for (int i = 0; i < p_.nx; ++i) {
        //     double x_tilde = (i - p_.source_pos) * invRes;
        //
        //     out << x_tilde;
        //     for (size_t k = 0; k < snapshotsEx_[k].size(); ++k) {
        //         out << "\t" << snapshotsEx_[k][i];
        //     }
        //     out << "\n";
        // }


        const double invRes = 1.0 / (double)p_.resolution;

        for (int i = 0; i < p_.nx; ++i) {
            double x_tilde = (i - p_.source_pos) * invRes;
            out << x_tilde;

            for (size_t k = 0; k < snapshotsEx_.size(); ++k) {
                out << "\t" << snapshotsEx_[k][i];
            }
            out << "\n";
        }


        out.close();
        std::cout << "Impulse snapshots written to " << filename << "\n";
    }
};



#endif //FDTD1D_H
