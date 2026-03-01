//
// Created by 6anna on 16.02.2026.
//

#ifndef NORMALIZATION_H
#define NORMALIZATION_H

#include "SimulationParameters.h"

struct Normalization {
    double a;       // [m]
    double c0;      // [m/s]
    double t_unit;  // [s] = a/c0

    explicit Normalization(double a_, double c0_ = 299792458.0)
        : a(a_), c0(c0_), t_unit(a_/c0_) {}

    // SI -> dimensionless
    double timeToDimless(double t_si) const { return t_si / t_unit; }
    double freqHzToDimless(double f_hz) const { return f_hz * t_unit; }
    double omegaToDimless(double omega_rad_s) const { return omega_rad_s * t_unit; }

    // dimensionless -> SI
    double freqDimlessToHz(double f_tilde) const { return f_tilde / t_unit; }
};

struct ParamsPack {
    SimulationParameters si;
    SimulationParameters nd;
    Normalization norm;
};


inline ParamsPack makeNormalizedParams(const SimulationParameters& p_si) {
    ParamsPack pack{p_si, p_si, Normalization(p_si.a)};

    // Сетка нормирована: dx~ = 1/resolution
    pack.nd.dx = 1.0 / p_si.resolution;
    pack.nd.dy = 1.0 / p_si.resolution;
    pack.nd.dz = 1.0 / p_si.resolution;

    // Время: dt~ = Q * dx~
    pack.nd.dt = p_si.courantNumber * pack.nd.dx;
    pack.nd.courantNumber = p_si.courantNumber;

    // Частоты: выбираем, что p_si хранит ω (rad/s)
    pack.nd.omega_p = pack.norm.omegaToDimless(p_si.omega_p);
    pack.nd.gamma   = pack.norm.omegaToDimless(p_si.gamma);

    // Источник ω (rad/s), фаза: omega*t
    pack.nd.sourceFreqCenter = pack.norm.omegaToDimless(p_si.sourceFreqCenter);
    pack.nd.sourceFreqWidth  = pack.norm.omegaToDimless(p_si.sourceFreqWidth);

    // pulseWidth: в секундах
    pack.nd.pulseWidth = pack.norm.timeToDimless(p_si.pulseWidth);

    pack.nd.nx = p_si.nx;
    pack.nd.ny = p_si.ny;
    pack.nd.nz = p_si.nz;
    pack.nd.numTimeSteps = p_si.numTimeSteps;

    pack.nd.pmlThickness = p_si.pmlThickness;
    pack.nd.pmlReflectCoeff = p_si.pmlReflectCoeff;
    pack.nd.pmlGrading_m = p_si.pmlGrading_m;

    pack.nd.source_pos = p_si.source_pos;
    pack.nd.monitorFront = p_si.monitorFront;
    pack.nd.monitorBack = p_si.monitorBack;
    pack.nd.plasmaStart = p_si.plasmaStart;
    pack.nd.plasmaWidth = p_si.plasmaWidth;

    pack.nd.epsilon_inf = p_si.epsilon_inf;
    pack.nd.sigmaCond   = p_si.sigmaCond;

    return pack;
}

//     ParamsPack pack{p_si, p_si, Normalization(p_si.dx)};
//
//     // Сетка: a=dx
//     pack.nd.dx = 1.0;
//     pack.nd.dy = 1.0;
//     pack.nd.dz = 1.0;
//
//     // Время: t~ = t / (a/c0), поэтому dt~ = Q*dx~ = Q
//     pack.nd.dt = pack.nd.courantNumber * pack.nd.dx; // = Q
//
//     pack.nd.omega_p = pack.norm.omegaToDimless(p_si.omega_p);
//     pack.nd.gamma   = pack.norm.omegaToDimless(p_si.gamma);
//     pack.nd.sourceFreqCenter = pack.norm.omegaToDimless(p_si.sourceFreqCenter);
//     pack.nd.sourceFreqWidth  = pack.norm.omegaToDimless(p_si.sourceFreqWidth);
//
//     pack.nd.pulseWidth = pack.norm.timeToDimless(p_si.pulseWidth);
//
//     return pack;
// }




#endif //NORMALIZATION_H
