#pragma once

#include "RK4.hpp"
#include "Star_Constants.hpp"
#include "Graphable.hpp"

class Star;

class Tau : public Graphable {
public:

	RK4 solver;
	Star* star;

	long double dtau_r(long double r, long double tau);

	Tau(Star* star) : star(star), solver(bind(&Tau::dtau_r, this, placeholders::_1, placeholders::_2), tau_0) {

	}

	void iterate() {
		solver.iterate();
	}

	long double get() {
		return solver.get();
	}

	virtual void pushValues() {
		arr[0].push_back(solver.t);
		arr[1].push_back(solver.y);
	}
};