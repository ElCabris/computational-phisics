#ifndef INTEGRATOR
#define INTEGRATOR

#include "integrator.h"
#include "state.h"
#include "system.h"
namespace simulacra {
class EulerCromerIntegrator : public IIntegrator {
	public:
	void step(IPhysicalSystem& system, State& state, double& t, double dt) const override {
		State dstate = system.derivatives(state, t);

		for (size_t i = 1; i< state.size(); i+=2) {
			state[i] += dstate[i] * dt;
		}

		for (size_t i = 0; i < state.size(); i += 2) {
			state[i] += state[i + 1] * dt;
		}

		t +=dt;
	}
};
}
#endif
