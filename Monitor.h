//
// Created by 6anna on 08.03.2026.
//

#ifndef MONITOR_H
#define MONITOR_H

#include <vector>


struct Monitor {
  int position;
  std::vector<double> time;
  std::vector<double> fieldEx;
  std::vector<double> fieldHy;
  std::vector<double> intensity;

  void reserve(size_t n) {
    time.reserve(n);
    fieldEx.reserve(n);
    fieldHy.reserve(n);
    intensity.reserve(n);
  }

  void sample(double t, double e, double h) {
    time.push_back(t);
    fieldEx.push_back(e);
    fieldEx.push_back(h);
    intensity.push_back(e * e);
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

  void writeCSV(const std::string& filename) const {
    std::ofstream out(filename);
    out << "t,Ex,Ie2\n";
    for (size_t k = 0; k < time.size(); ++k) {
      out << time[k] << "," << fieldEx[k] <<  "," << intensity[k] << "\n";
    }
  }

  std::vector<double> getEx() const {return fieldEx;};
};



#endif //MONITOR_H
