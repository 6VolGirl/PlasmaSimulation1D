//
// Created by 6anna on 10.03.2026.
//

#include "SpectrumAnalyzer.h"


std::vector<std::complex<double>> SpectrumAnalyzer::fftRecursive(const std::vector<std::complex<double>>& a) {
    /* Функция принимает комплексный вектор и возвращает его FFT */

    int n = (int)a.size();
    // для рекурсии
    if (n == 1) {
        return a;
    }

    std::vector<std::complex<double>> even(n / 2);  // чётные
    std::vector<std::complex<double>> odd(n / 2);   // нечётные

    for (int i = 0; i < n / 2; ++i) {
        even[i] = a[2 * i];
        odd[i]  = a[2 * i + 1];
    }

    // рекурсия
    auto Fe = fftRecursive(even);
    auto Fo = fftRecursive(odd);

    std::vector<std::complex<double>> y(n); // итоговый результат

    // комплексное превращения
    for (int k = 0; k < n / 2; ++k) {
        double ang = -2.0 * M_PI * k / n;
        std::complex<double> w(std::cos(ang), std::sin(ang));

        y[k] = Fe[k] + w * Fo[k];
        y[k + n / 2] = Fe[k] - w * Fo[k];
    }

    return y;
}

std::vector<std::complex<double>> SpectrumAnalyzer::computeFFT(const std::vector<double>& signal, double dt) {
    if (signal.empty()) {
        throw std::runtime_error("computeFFT: signal is empty");
    }
    if (dt <= 0.0) {
        throw std::runtime_error("computeFFT: dt must be > 0");
    }

    int N0 = (int)signal.size();

    int N = 1;
    while (N < N0) {
        N = N * 2;
    }

    std::vector<std::complex<double>> a(N, std::complex<double>(0.0, 0.0));
    // реальный сигнал в комплексный
    for (int i = 0; i < N0; ++i) {
        a[i] = std::complex<double>(signal[i], 0.0);
    }

    // рекурсия
    auto fullSpectrum = fftRecursive(a);

    // Вычисляем число уникальных точек спектра для вещественного сигнала
    int Nh = N / 2 + 1;
    freq.resize(Nh);

    std::vector<std::complex<double>> spectrum(Nh);
    // укороченный спектр(неотрицательные частоты)
    for (int k = 0; k < Nh; ++k) {
        freq[k] = k / (N * dt);
        spectrum[k] = fullSpectrum[k];
    }

    return spectrum;
}



//void SpectrumAnalyzer::writeSpectrumCSV(const std::string& filename, double fL) const {
//    if (freq.empty() || poynting_f.empty()) {
//        throw std::runtime_error("writeSpectrumCSV: spectrum is empty");
//    }
//    if (freq.size() != poynting_f.size()) {
//        throw std::runtime_error("writeSpectrumCSV: freq and spectrum sizes differ");
//    }
//    if (fL <= 0.0) {
//        throw std::runtime_error("writeSpectrumCSV: fL must be > 0");
//    }
//
//    std::ofstream out(filename);
//    if (!out.is_open()) {
//        throw std::runtime_error("writeSpectrumCSV: cannot open output file");
//    }
//
//    out << std::scientific << std::setprecision(10);
//    out << "f,f_over_fL,ReS,ImS,absS2\n";
//
//    for (std::size_t k = 0; k < freq.size(); ++k) {
//        const double f = freq[k];
//        const double absS2 = std::norm(poynting_f[k]);
//
//        out << f << ','
//            << f / fL << ','
//            << poynting_f[k].real() << ','
//            << poynting_f[k].imag() << ','
//            << absS2 << '\n';
//    }
//}

