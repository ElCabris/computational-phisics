#include "../include/physim/core/state.h"
#include <cmath>
#include <matplot/matplot.h>
#include <string>
#include <vector>

using namespace matplot;
using simulacra::State;

State pendulumSystem(const State &s, double t, double F_D) {
  double theta = s[0];
  double omega = s[1];
  double g = 9.8, l = 9.8, q = 0.5, Omega_D = 2.0 / 3.0;
  double dtheta = omega;
  double domega =
      -(g / l) * std::sin(theta) - q * omega + F_D * std::sin(Omega_D * t);
  return State({dtheta, domega});
}

void simulate_phase_space(double F_D, const std::string &filename) {
  double dt = 0.04;
  double tmax = 300.0;
  State s({0.2, 0.0});

  std::vector<double> thetas, omegas;

  double t = 0.0;
  while (t <= tmax) {
    State ds = pendulumSystem(s, t, F_D);

    // Integración Euler-Cromer
    s[1] += ds[1] * dt;
    s[0] += s[1] * dt;

    if (s[0] > M_PI)
      s[0] -= 2 * M_PI;
    if (s[0] < -M_PI)
      s[0] += 2 * M_PI;

    thetas.push_back(s[0]);
    omegas.push_back(s[1]);

    t += dt;
  }

  matplot::figure();
  matplot::plot(thetas, omegas, "-k")->line_width(1.0);
  matplot::title("Espacio de fase ω vs θ | F_D = " + std::to_string(F_D));
  matplot::xlabel("θ (radianes)");
  matplot::ylabel("ω (rad/s)");
  matplot::grid(true);
  matplot::save(filename);
}

int main() {
  simulate_phase_space(0.5, "phase_space_FD05.png");
  simulate_phase_space(1.2, "phase_space_FD12.png");
  return 0;
}
