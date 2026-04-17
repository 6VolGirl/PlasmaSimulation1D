//
// Created by 6anna on 10.03.2026.
//

#include "SpectrumAnalyzer.h"

#include <stdexcept>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

namespace {
    constexpr double PI = 3.14159265358979323846;
}

std::vector<std::complex<double>> SpectrumAnalyzer::fftRecursive(
    const std::vector<std::complex<double>>& a
) {
    const int n = static_cast<int>(a.size());

    if (n == 1) {
        return a;
    }

    std::vector<std::complex<double>> even(n / 2);
    std::vector<std::complex<double>> odd(n / 2);

    for (int i = 0; i < n / 2; ++i) {
        even[i] = a[2 * i];
        odd[i]  = a[2 * i + 1];
    }

    auto Fe = fftRecursive(even);
    auto Fo = fftRecursive(odd);

    std::vector<std::complex<double>> y(n);

    for (int k = 0; k < n / 2; ++k) {
        const double ang = -2.0 * PI * static_cast<double>(k) / static_cast<double>(n);
        const std::complex<double> w(std::cos(ang), std::sin(ang));

        y[k]         = Fe[k] + w * Fo[k];
        y[k + n / 2] = Fe[k] - w * Fo[k];
    }

    return y;
}

std::vector<std::complex<double>> SpectrumAnalyzer::computeFFT(
    const std::vector<double>& signal,
    double dt
) {
    if (signal.size() < 2) {
        throw std::runtime_error("computeFFT: signal must contain at least 2 samples");
    }
    if (dt <= 0.0) {
        throw std::runtime_error("computeFFT: dt must be > 0");
    }

    const int Nsig = static_cast<int>(signal.size());

    int N = 1;
    while (N < Nsig) {
        N <<= 1;
    }

    std::vector<std::complex<double>> a(N, std::complex<double>(0.0, 0.0));

    // FIX: вычитаем среднее, чтобы убрать паразитный пик около f = 0
    double mean = 0.0;
    for (int i = 0; i < Nsig; ++i) {
        mean += signal[i];
    }
    mean /= static_cast<double>(Nsig);

    // FIX: записываем сигнал с вычитанием среднего + zero padding
    for (int i = 0; i < N; ++i) {
        if (i < Nsig) {
            a[i] = std::complex<double>(signal[i] - mean, 0.0);
        } else {
            a[i] = std::complex<double>(0.0, 0.0);
        }
    }

    auto fullSpectrum = fftRecursive(a);

    // FIX: берём только неотрицательные частоты
    const int Nh = N / 2 + 1;
    freq.resize(Nh);

    std::vector<std::complex<double>> spectrum(Nh);

    const double df = 1.0 / (static_cast<double>(N) * dt);

    // FIX: масштаб dt полезен, если хочешь, чтобы амплитуда была ближе к непрерывному Фурье-интегралу
    const double scale = dt;

    for (int k = 0; k < Nh; ++k) {
        freq[k] = static_cast<double>(k) * df;
        spectrum[k] = fullSpectrum[k] * scale;
    }

    return spectrum;
}

void SpectrumAnalyzer::buildFluxSpectrum(const Monitor& m) {
    const size_t N = std::min(m.time.size(), m.poynting.size());

    if (N < 2) {
        throw std::runtime_error("buildFluxSpectrum: not enough monitor samples");
    }

    std::vector<double> signal(N);
    for (size_t i = 0; i < N; ++i) {
        signal[i] = m.poynting[i];
    }

    // Надёжнее брать dt по времени монитора
    const double dt = (m.time[N - 1] - m.time[0]) / static_cast<double>(N - 1);
    if (dt <= 0.0) {
        throw std::runtime_error("buildFluxSpectrum: invalid dt from monitor time");
    }

    // FIX: спектр строится прямо из уже сохранённого S(t),
    // а не пересчитывается заново из fieldEx/fieldHy/fieldHy2
    poynting_f = computeFFT(signal, dt);
}

void SpectrumAnalyzer::writeSpectrumCSV(const std::string& filename, double fL) const {
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

        out << f << ","
            << f / fL << ","
            << poynting_f[k].real() << ","
            << poynting_f[k].imag() << ","
            << absS2 << "\n";
    }
}
// std::vector<std::complex<double>> SpectrumAnalyzer::fftRecursive(const std::vector<std::complex<double>>& a) {
//     /* Функция принимает комплексный вектор и возвращает его FFT */
//
//     int n = (int)a.size();
//     // для рекурсии
//     if (n == 1) {
//         return a;
//     }
//
//     std::vector<std::complex<double>> even(n / 2);  // чётные
//     std::vector<std::complex<double>> odd(n / 2);   // нечётные
//
//     for (int i = 0; i < n / 2; ++i) {
//         even[i] = a[2 * i];
//         odd[i]  = a[2 * i + 1];
//     }
//
//     // рекурсия
//     auto Fe = fftRecursive(even);
//     auto Fo = fftRecursive(odd);
//
//     std::vector<std::complex<double>> y(n); // итоговый результат
//
//     // комплексное превращения
//     for (int k = 0; k < n / 2; ++k) {
//         double ang = -2.0 * M_PI * k / n;
//         std::complex<double> w(std::cos(ang), std::sin(ang));
//
//         y[k] = Fe[k] + w * Fo[k];
//         y[k + n / 2] = Fe[k] - w * Fo[k];
//     }
//
//     return y;
// }
//
// std::vector<std::complex<double>> SpectrumAnalyzer::computeFFT(const std::vector<double>& signal, double dt) {
//     if (signal.empty()) {
//         throw std::runtime_error("computeFFT: signal is empty");
//     }
//     if (dt <= 0.0) {
//         throw std::runtime_error("computeFFT: dt must be > 0");
//     }
//
//      const int N0 = static_cast<int>(signal.size());
//
//     int N = 1;
//     while (N < N0) {
//         N = N * 2;
//     }
//
//     std::vector<std::complex<double>> a(N, std::complex<double>(0.0, 0.0));
//     // реальный сигнал в комплексный
//     for (int i = 0; i < N0; ++i) {
//         a[i] = std::complex<double>(signal[i], 0.0);
//     }
//
//     // рекурсия
//     auto fullSpectrum = fftRecursive(a);
//
//     // Вычисляем число уникальных точек спектра для вещественного сигнала
//     const int Nh = N / 2 + 1;
//     freq.resize(Nh);
//
//     std::vector<std::complex<double>> spectrum(Nh);
//     // укороченный спектр(неотрицательные частоты)
//     for (int k = 0; k < Nh; ++k) {
//         freq[k] = k / (N * dt);
//         spectrum[k] = fullSpectrum[k];
//     }
//
//     return spectrum;
// }
//


