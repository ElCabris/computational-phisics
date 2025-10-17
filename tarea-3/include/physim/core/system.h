#ifndef SYSTEM_H
#define SYSTEM_H


#include "state.h"
namespace simulacra {
class IPhysicalSystem {
public:
	virtual ~IPhysicalSystem() = default;
	virtual State derivatives(const State& state, double t) const = 0;
};

}

#endif
