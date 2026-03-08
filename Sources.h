//
// Created by 6anna on 03.03.2026.
//

#ifndef SOURCES_H
#define SOURCES_H

#include <cmath>

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

    // double operator()(double t) const {
    //     if (t < start_time || t > finish_time) return 0.0;
    //     double tau = t - t0;
    //     double env = std::exp(-0.5 * tau * tau / (w * w));
    //     double phase = 2.0 * M_PI * (freq * tau + 0.5 * chirpRate * tau * tau);
    //
    //     return env * std::sin(phase);
    // }


    double operator()(double t) const {
        if (t < start_time || t > finish_time) return 0.0;
        double tau = t - t0;
        return std::exp(-0.5 * (tau*tau) / (w*w)) * std::sin(2.0*M_PI*freq*t);
    }
};

#endif //SOURCES_H
