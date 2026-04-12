//
// Created by 6anna on 08.03.2026.
//

#ifndef MONITOR_H
#define MONITOR_H

#include <vector>
#include <string>


struct Monitor {
  int position;
  std::vector<double> time;
  std::vector<double> fieldEx;
  std::vector<double> fieldHy;    // Hy (pos - 1/2)
  std::vector<double> fieldHy2;   // Hy (pos + 1/2)
  std::vector<double> intensity;
  std::vector<double> poynting;   // S = E * H_interp

  void reserve(size_t n) {
    time.reserve(n);
    fieldEx.reserve(n);
    fieldHy.reserve(n);
    fieldHy2.reserve(n);
    intensity.reserve(n);
    poynting.reserve(n);
  }

  void sample(double t, double e, double h, double h2) {
    time.push_back(t);
    fieldEx.push_back(e);
    fieldHy.push_back(h);
    fieldHy2.push_back(h2);
    intensity.push_back(e * e);

    double h_interp = 0.5 * (h + h2);
    poynting.push_back(e * h_interp);
  }

  // 𝑡_𝑐 = ∑ (𝑡_𝑘 * 𝐼_𝑘) / ∑𝐼_𝑘    - центр тяжести
  double centroidTime() const {
    double num = 0.0, den = 0.0;
    for (size_t k = 0; k < time.size(); ++k) {
      num += time[k] * intensity[k];
      den += intensity[k];
    }
    return (den > 0.0) ? num / den : 0.0;
  }

  double centroidTimeByPoynting() const {
    double num = 0.0, den = 0.0;
    for (size_t k = 0; k < time.size(); ++k) {
      double w = std::abs(poynting[k]); // FIX
      num += time[k] * w;
      den += w;
    }
    return (den > 0.0) ? num / den : 0.0;
  }

  double integratedFlux() const {
    if (time.size() < 2) return 0.0;

    double dt = time[1] - time[0];
    double sum = 0.0;

    for (size_t k = 0; k < poynting.size(); ++k) {
      sum += poynting[k];
    }

    return sum * dt;
  }

  void writeCSV(const std::string& filename) const {
    std::ofstream out(filename);
    out << "t,Ex,Hy,Ie2,S\n";
    for (size_t k = 0; k < time.size(); ++k) {
      out << time[k] << ","
          << fieldEx[k] << ","
          << fieldHy[k] << ","
          << intensity[k] << ","
          << poynting[k] << "\n";
    }
  }

  std::vector<double> getEx() const {return fieldEx;};
};



#endif //MONITOR_H
