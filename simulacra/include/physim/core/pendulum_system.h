#ifndef PENDULUM_SYSTEM_H
#define PENDULUM_SYSTEM_H

#include <cmath>
#include "system.h"
#include "state.h"

namespace simulacra {

class PendulumSystem : public IPhysicalSystem {
private:
    double g;
    double l;
    double q;
    double F_D;
    double Omega_D;

public:
    PendulumSystem(double g_ = 9.8, double l_ = 9.8, double q_ = 0.5,
                   double F_D_ = 1.2, double Omega_D_ = 2.0 / 3.0)
        : g(g_), l(l_), q(q_), F_D(F_D_), Omega_D(Omega_D_) {}

    State derivatives(const State& state, double t) const override {
        double theta = state[0];
        double omega = state[1];

        double dtheta = omega;
        double domega = - (g / l) * std::sin(theta) - q * omega + F_D * std::sin(Omega_D * t);

        return State({dtheta, domega});
    }
};

}

#endif

