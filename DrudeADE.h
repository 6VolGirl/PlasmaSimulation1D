 //
 // Created by 6anna on 15.02.2026.
 //

 #ifndef DRUDEADE_H
 #define DRUDEADE_H
 #include <vector>
 #include "SimulationParameters.h"


 class DrudeADE {
 private:
     const SimulationParameters& p_;
     std::vector<double> alpha_, xi_, gcoef_;  // ADE
     std::vector<double> Jn_, Jnm1_;     // Поляризация: n and n-1

 public:
     explicit DrudeADE(const SimulationParameters& p)
     : p_(p),
       alpha_(p.nx + 1, 0.0),
       xi_(p.nx + 1, 0.0),
       gcoef_(p.nx + 1, 0.0),
       Jn_(p.nx + 1, 0.0),
       Jnm1_(p.nx + 1, 0.0)
     {
         // коэффициенты в области плазмы
         const int i0 = std::max(0, p_.plasmaStart);
         const int i1 = std::min(p_.nx, p_.plasmaStart + p_.plasmaWidth);

         const double dt = p_.dt;
         const double denom = 1.0 + p_.gamma * dt / 2.0;

         for (int i = i0; i <= i1; ++i) {
             alpha_[i] = 2.0 / denom;
             xi_[i]    = (p_.gamma * dt / 2.0 - 1.0) / denom;
             gcoef_[i] = (p_.omega_p * p_.omega_p) * (dt * dt) / denom;
         }
     }


     // 0.5*((1+alpha)*J^n + xi*J^{n-1})
     inline double polCorrection(int i) const {
         return 0.5 * ((1.0 + alpha_[i]) * Jn_[i] + xi_[i] * Jnm1_[i]);
     }

     void updateJ(const std::vector<double>& E_np1,
                  const std::vector<double>& E_nm1)
     {
         const double dt = p_.dt;
         const int i0 = std::max(0, p_.plasmaStart);
         const int i1 = std::min(p_.nx, p_.plasmaStart + p_.plasmaWidth);

         for (int i = i0; i <= i1; ++i) {
             const double dE_dt = (E_np1[i] - E_nm1[i]) / (2.0 * dt);
             const double J_next =
                 alpha_[i] * Jn_[i] + xi_[i] * Jnm1_[i] + gcoef_[i] * dE_dt;

             Jnm1_[i] = Jn_[i];
             Jn_[i]   = J_next;
         }
     }
 };



 #endif //DRUDEADE_H
