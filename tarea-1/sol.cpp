#include <math.h>
#include <matplot/matplot.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#define ELEMENTS 3
#define MAX_TIME 10000 // dias

struct Simulation {
  double h;
  double lambda;
  int iterations; // number of times the approximation will be made
  double *results;
  double *time_points;
};

void save_simulation_plot(const Simulation &sim, const std::string &filename);
// df(n + 1) = - \lambda  f(n)
double df(double lambda, double f_n);
// f(n + 1) = f(n) + h * df(n)
double euler_aprox(double f_n, double df_n, double h);
// f(t) = f0 e^{-\lambda t}
double analytical_solution(double f0, double lambda, double t);
void free_simulation(struct Simulation *simulation);

int main(void) {

  // number of atoms in 2g of:
  const double atoms[] = {
      5.059542800944304e+21, // U-238
      8.601066473464549e+22, // C-14
      2.0095969372253427e+22 // Co-60
  };

  // steps for Euler's approximation
  const double steps[] = {0.5, 0.3, 0.1};
  const int num_steps = sizeof(steps) / sizeof(steps[0]);

  // half live in seconds
  const double half_life[] = {
      4.468e9, // U-238 en dias
      5730,    // C-14
      5.27     // Co-60
  };

  double lambda[ELEMENTS];
  for (int i = 0; i < ELEMENTS; ++i)
    lambda[i] = log(2) / half_life[i];

  struct Simulation sims[ELEMENTS][num_steps];

  // initializacitons of simulations
  for (int i = 0; i < ELEMENTS; ++i) {
    for (int j = 0; j < num_steps; ++j) {
      sims[i][j].h = steps[j];
      sims[i][j].iterations = (int)ceil(MAX_TIME / steps[j]) + 1;
      sims[i][j].lambda = lambda[i];

      sims[i][j].results =
          (double *)malloc(sims[i][j].iterations * sizeof(double));
      sims[i][j].time_points =
          (double *)malloc(sims[i][j].iterations * sizeof(double));

      if (!sims[i][j].results || !sims[i][j].time_points) {
        fprintf(stderr, "error allocating memory\n");
        exit(EXIT_FAILURE);
      }

      sims[i][j].results[0] = atoms[i];
      sims[i][j].time_points[0] = 0.0;
    }
  }

  for (int element = 0; element < ELEMENTS; ++element) {
    for (int step = 0; step < num_steps; ++step) {
      for (int k = 1; k < sims[element][step].iterations; ++k) {
        double t_prev = sims[element][step].time_points[k - 1],
               f_prev = sims[element][step].results[k - 1],

               deriv = df(sims[element][step].lambda, f_prev),
               f_next = euler_aprox(f_prev, deriv, sims[element][step].h);

        sims[element][step].time_points[k] = t_prev + sims[element][step].h;
        sims[element][step].results[k] = f_next;
      }
    }
  }

  for (int i = 0; i < ELEMENTS; ++i) {
    for (int j = 0; j < num_steps; ++j) {
      save_simulation_plot(sims[i][j], "simulation-" + std::to_string(i) + "-" +
                                           std::to_string(j));
      free_simulation(sims[i] + j);
    }
  }
  return 0;
}

double df(double lambda, double f_n) { return -1 * lambda * f_n; }

double euler_aprox(double f_n, double df_n, double h) {
  return f_n + (df_n * h);
}

double analytical_solution(double f0, double lambda, double t) {
  return f0 * exp(-1 * lambda * t);
}

void save_simulation_plot(const Simulation &sim, const std::string &filename) {
  // Configurar el backend (no siempre necesario en versiones recientes)
  // matplot::backend::gnuplot();

  // Crear figura sin mostrar
  auto fig = matplot::figure(true);
  fig->quiet_mode(true);

  // Convertir arrays a vectores
  std::vector<double> t(sim.time_points, sim.time_points + sim.iterations);
  std::vector<double> euler(sim.results, sim.results + sim.iterations);

  // Calcular solución analítica
  std::vector<double> analytical;
  const double f0 = sim.results[0];
  for (const auto &time : t) {
    analytical.push_back(f0 * std::exp(-sim.lambda * time));
  }

  // Configurar gráfica
  auto ax = fig->current_axes();
  matplot::hold(ax, true);
  matplot::plot(ax, t, euler, "b-")->line_width(1.5).display_name("Euler");
  matplot::plot(ax, t, analytical, "r--")
      ->line_width(1.5)
      .display_name("Analítico");

  // Añadir etiquetas
  matplot::xlabel(ax, "Tiempo (t)");
  matplot::ylabel(ax, "Átomos restantes");
  matplot::title(ax, "Simulación (h = " + matplot::num2str(sim.h) + ")");
  matplot::legend(ax, "show");

  // Guardar imagen
  matplot::save(filename + ".png");

  // Liberar recursos (alternativa sin clf)
  ax->clear();
}
void free_simulation(struct Simulation *simulation) {
  free(simulation->results);
  free(simulation->time_points);
  simulation->results = NULL;
  simulation->time_points = NULL;
}
