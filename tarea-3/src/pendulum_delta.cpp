#include "../include/physim/core/state.h"
#include <cmath>
#include <matplot/matplot.h>
#include <string>
#include <vector>

using namespace matplot;
using simulacra::State;

// Sistema físico: péndulo no lineal amortiguado y forzado
State pendulumSystem(const State &s, double t, double F_D) {
  double theta = s[0];
  double omega = s[1];

  // Parámetros físicos
  double g = 9.8, l = 9.8, q = 0.5, Omega_D = 2.0 / 3.0;
  double dtheta = omega;
  double domega =
      -(g / l) * std::sin(theta) - q * omega + F_D * std::sin(Omega_D * t);

  return State({dtheta, domega});
}

void simulate_delta_theta(double F_D, const std::string &filename) {
  double dt = 0.04;
  double tmax = 200.0;

  // Condiciones iniciales
  State s1({0.2, 0.0});
  State s2({0.2001, 0.0}); // ligeramente distinto

  std::vector<double> time, delta_theta;

  double t = 0.0;
  while (t <= tmax) {
    // Calcular derivadas
    State ds1 = pendulumSystem(s1, t, F_D);
    State ds2 = pendulumSystem(s2, t, F_D);

    // Integrar (Euler-Cromer)
    s1[1] += ds1[1] * dt;
    s1[0] += s1[1] * dt;

    s2[1] += ds2[1] * dt;
    s2[0] += s2[1] * dt;

    // Mantener θ en [-π, π]
    if (s1[0] > M_PI)
      s1[0] -= 2 * M_PI;
    if (s1[0] < -M_PI)
      s1[0] += 2 * M_PI;
    if (s2[0] > M_PI)
      s2[0] -= 2 * M_PI;
    if (s2[0] < -M_PI)
      s2[0] += 2 * M_PI;

    double dtheta = std::fabs(s1[0] - s2[0]);
    if (dtheta > M_PI)
      dtheta = 2 * M_PI - dtheta; // continuidad angular

    time.push_back(t);
    delta_theta.push_back(dtheta);

    t += dt;
  }

  // Graficar Δθ(t)
  matplot::figure();
  matplot::semilogy(time, delta_theta, "-b")->line_width(1.5);
  matplot::title("Evolución de Δθ(t)  |  F_D = " + std::to_string(F_D));
  matplot::xlabel("Tiempo (s)");
  matplot::ylabel("Δθ (rad)");
  matplot::grid(true);
  matplot::save(filename);
}

int main() {
  simulate_delta_theta(0.5, "delta_theta_FD05.png"); // régimen periódico
  simulate_delta_theta(1.2, "delta_theta_FD12.png"); // régimen caótico
  return 0;
}
