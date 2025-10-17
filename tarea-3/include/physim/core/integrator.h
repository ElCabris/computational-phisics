#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "system.h"
namespace simulacra {

	class IIntegrator {
	public:
	virtual ~IIntegrator() = default;

	virtual void step(IPhysicalSystem& system, State& state, double& t, double dt) const = 0;

	};
}

#endif
