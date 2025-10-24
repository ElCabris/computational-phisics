#include "../include/physim/core/state.h"
#include <cmath>
#include <matplot/matplot.h>
#include <vector>

using namespace matplot;

namespace simulacra {

// Modelo físico del péndulo no lineal, amortiguado y accionado
State pendulumSystem(const State &s, double t, double F_D) {
  double theta = s[0];
  double omega = s[1];

  // Parámetros físicos
  double g = 9.8, l = 9.8, q = 0.5, Omega_D = 2.0 / 3.0;

  double dtheta = omega;
  double domega =
      -(g / l) * std::sin(theta) - q * omega + F_D * std::sin(Omega_D * t);

  return State(std::vector<double>{dtheta, domega});
}

} // namespace simulacra

int main() {
  using namespace simulacra;

  double dt = 0.04;
  double tmax = 200.0;

  std::vector<double> F_D_values = {0.0, 0.5, 1.2};

  // === FIGURA 1: θ(t) con reinicio para tres F_D ===
  auto f1 = figure(true);
  hold(on);
  for (double F_D : F_D_values) {
    std::vector<double> time, theta_values;
    State state(std::vector<double>{0.2, 0.0});
    double t = 0.0;

    while (t <= tmax) {
      State deriv = pendulumSystem(state, t, F_D);
      state[1] += deriv[1] * dt;
      state[0] += state[1] * dt;

      // Reiniciar θ al rango [-π, π]
      if (state[0] > M_PI)
        state[0] -= 2 * M_PI;
      if (state[0] < -M_PI)
        state[0] += 2 * M_PI;

      time.push_back(t);
      theta_values.push_back(state[0]);
      t += dt;
    }
    plot(time, theta_values)->display_name("F_D = " + std::to_string(F_D));
  }
  title("θ(t) con reinicio en [-π, π] para varios F_D");
  xlabel("Tiempo (s)");
  ylabel("θ (rad)");
  legend();
  grid(on);

  // === FIGURA 2: θ(t) con y sin reinicio para F_D = 1.2 ===
  double F_D = 1.2;
  auto f2 = figure(true);
  hold(on);

  // --- Con reinicio ---
  std::vector<double> time_r, theta_r;
  State s1(std::vector<double>{0.2, 0.0});
  double t = 0.0;
  while (t <= tmax) {
    State deriv = pendulumSystem(s1, t, F_D);
    s1[1] += deriv[1] * dt;
    s1[0] += s1[1] * dt;

    if (s1[0] > M_PI)
      s1[0] -= 2 * M_PI;
    if (s1[0] < -M_PI)
      s1[0] += 2 * M_PI;

    time_r.push_back(t);
    theta_r.push_back(s1[0]);
    t += dt;
  }

  // --- Sin reinicio ---
  std::vector<double> time_nr, theta_nr;
  State s2(std::vector<double>{0.2, 0.0});
  t = 0.0;
  while (t <= tmax) {
    State deriv = pendulumSystem(s2, t, F_D);
    s2[1] += deriv[1] * dt;
    s2[0] += s2[1] * dt;
    time_nr.push_back(t);
    theta_nr.push_back(s2[0]);
    t += dt;
  }

  plot(time_r, theta_r)->display_name("Con reinicio");
  plot(time_nr, theta_nr)->display_name("Sin reinicio");
  title("θ(t) para F_D = 1.2 con y sin reinicio");
  xlabel("Tiempo (s)");
  ylabel("θ (rad)");
  legend();
  grid(on);

  // === FIGURA 3: ω(t) para F_D = 1.2 ===
  auto f3 = figure(true);
  hold(on);

  std::vector<double> time_w, omega_values;
  State s3(std::vector<double>{0.2, 0.0});
  t = 0.0;
  while (t <= tmax) {
    State deriv = pendulumSystem(s3, t, F_D);
    s3[1] += deriv[1] * dt;
    s3[0] += s3[1] * dt;

    if (s3[0] > M_PI)
      s3[0] -= 2 * M_PI;
    if (s3[0] < -M_PI)
      s3[0] += 2 * M_PI;

    time_w.push_back(t);
    omega_values.push_back(s3[1]);
    t += dt;
  }

  plot(time_w, omega_values, "r-")->display_name("ω(t)");
  title("Velocidad angular ω(t) para F_D = 1.2");
  xlabel("Tiempo (s)");
  ylabel("ω (rad/s)");
  legend();
  grid(on);

  show();
  return 0;
}
