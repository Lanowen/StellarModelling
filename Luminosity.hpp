#pragma once

#include "RK4.hpp"
#include "Star_Constants.hpp"
#include "Graphable.hpp"

class Star;

class Luminosity : public Graphable {
public:

	RK4 solver;
	Star* star;
	

	long double dL_r(long double r, long double L);

	Luminosity(Star* star) : star(star), solver(bind(&Luminosity::dL_r, this, placeholders::_1, placeholders::_2), L_0) {

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