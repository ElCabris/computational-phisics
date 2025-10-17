#ifndef STATE_H 
#define STATE_H

#include <cstddef>
#include <vector>
namespace simulacra {
class State {
private:
	std::vector<double> data;
public:
	State() = default;
	explicit State(std::vector<double> values) : data(std::move(values)) {}

	double& operator[](size_t i) { return data[i]; } 
	const double& operator[](size_t i) const { return data[i]; }
	size_t size() const { return data.size(); }
	std::vector<double>& values() { return data; }
	const std::vector<double>& values() const { return data; }
};
}

#endif
