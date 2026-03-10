//
// Created by 6anna on 10.03.2026.
//

#ifndef SPECTRUMANALYZER_H
#define SPECTRUMANALYZER_H

#include <vector>
#include <cmath>
#include <complex>
#include <iomanip>
#include <fstream>

#include "Monitor.h"

class SpectrumAnalyzer {
    std::vector<double> freq;   // ось частот
    std::vector<std::complex<double>>  poynting_f;  // FluxSpectrum S(f)


    std::vector<std::complex<double>> fftRecursive(const std::vector<std::complex<double>>& a);

    std::vector<std::complex<double>> computeFFT(const std::vector<double>& signal, double dt);

public:


    void buildFluxSpectrum(const Monitor& m){
        if (m.fieldEx.size() != m.fieldHy.size() || m.fieldEx.size() != m.time.size()) {
            throw std::runtime_error("buildFluxSpectrum: monitor arrays have different sizes");
        }

        if (m.time.size() < 2) {
            throw std::runtime_error("buildFluxSpectrum: not enough time samples");
        }


        int nStep = m.fieldEx.size();
        const double dt = m.time[1] - m.time[0];
        if (dt <= 0.0) {
            throw std::runtime_error("buildFluxSpectrum: invalid dt reconstructed from monitor time");
        }

        std::vector<double> poynting_t(nStep-1, 0);   // S(t) = E(t)*H(t)

        // Так как H хранить на полшага сдвинутое относительно E(t)
        // то H(t)≈ 0.5 * ( H_n−1/2 + H_n+1/2)
        for (int i = 0; i < m.fieldEx.size()-1; i++){
            double Hy_i = 0.5 * (m.fieldHy[i] + m.fieldHy[i+1]);
            poynting_t[i] = m.fieldEx[i] * Hy_i;
        }
        poynting_f = computeFFT(poynting_t, dt);

    }

    void writeSpectrumCSV(const std::string& filename, double fL) const {
        if (freq.empty() || poynting_f.empty()) {
            throw std::runtime_error("writeSpectrumCSV: spectrum is empty");
        }
        if (freq.size() != poynting_f.size()) {
            throw std::runtime_error("writeSpectrumCSV: freq and spectrum sizes differ");
        }
        if (fL <= 0.0) {
            throw std::runtime_error("writeSpectrumCSV: fL must be > 0");
        }

        std::ofstream out(filename);
        if (!out.is_open()) {
            throw std::runtime_error("writeSpectrumCSV: cannot open output file");
        }

        out << std::scientific << std::setprecision(10);
        out << "f,f_over_fL,ReS,ImS,absS2\n";

        for (std::size_t k = 0; k < freq.size(); ++k) {
            const double f = freq[k];
            const double absS2 = std::norm(poynting_f[k]);

            out << f << ','
                << f / fL << ','
                << poynting_f[k].real() << ','
                << poynting_f[k].imag() << ','
                << absS2 << '\n';
        }
    }


};



#endif //SPECTRUMANALYZER_H
