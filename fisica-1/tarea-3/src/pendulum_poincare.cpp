#include "../include/physim/core/euler_cromer.h"
#include "../include/physim/core/pendulum_system.h"
#include "../include/physim/core/state.h"
#include <iostream>
#include <matplot/matplot.h>
#include <vector>

using namespace matplot;

int main() {
  using namespace simulacra;

  const double dt = 0.04;
  const double Omega_D = 2.0 / 3.0;
  const double two_pi = 2 * M_PI;
  const double t_max = 1000.0;
  const int n_steps = static_cast<int>(t_max / dt);

  // Crear sistema físico
  PendulumSystem system(9.8, 9.8, 0.5, 1.2, 2.0 / 3.0);

  // Integrador
  EulerCromerIntegrator integrator;

  // Estado inicial
  State state({0.2, 0.0});
  double t = 0.0;

  std::vector<double> theta_points;
  std::vector<double> omega_points;

  // Variables de fase
  double next_phase_time = two_pi / Omega_D;
  int phase_index = 1;

  for (int i = 0; i < n_steps; ++i) {
    integrator.step(system, state, t, dt);

    // Mantener theta en rango [-pi, pi]
    if (state[0] > M_PI)
      state[0] -= 2 * M_PI;
    if (state[0] < -M_PI)
      state[0] += 2 * M_PI;

    // Guardar puntos en fase con la fuerza impulsora
    if (t >= phase_index * next_phase_time) {
      theta_points.push_back(state[0]);
      omega_points.push_back(state[1]);
      phase_index++;
    }
  }

  // Graficar y guardar figura
  figure(true);
  scatter(theta_points, omega_points, 5.0);
  title("Poincaré Map for Driven Pendulum (F_D = 1.2)");
  xlabel("\\theta (radians)");
  ylabel("\\omega (radians/s)");
  xlim({-4, 4});
  ylim({-2, 1});
  grid(true);
  save("poincare_map.png");

  std::cout << "✅ Figura guardada como poincare_map.png\n";
  return 0;
}
