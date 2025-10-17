#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "integrator.h"
#include "state.h"
#include "system.h"
#include <cstdlib>
#include <vector>
namespace simulacra {
class Simulator {
private:
	IIntegrator& integrator;
	IPhysicalSystem& system;
	double dt;
public:
	Simulator(IIntegrator& integrator_, IPhysicalSystem& system_, double dt_) : integrator(integrator_), system(system_), dt(dt_) {}

	std::vector<State> run(State initialState, double t0, double tf) {
		std::vector<State> results;

		State state = initialState;
		double t = t0;
		while(t < tf) {
			integrator.step(system, state, t, dt);
			results.push_back(state);
		}
		return results;
	}

};
}

#endif
