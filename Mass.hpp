#pragma once

#include "RK4.hpp"
#include "Star_Constants.hpp"
#include "Graphable.hpp"

class Star;

class Mass : public Graphable {
public:

	RK4 solver;
	Star* star;
	

	long double dM_dr(long double r, long double M);

	Mass(Star* star) : star(star), solver(bind(&Mass::dM_dr, this, placeholders::_1, placeholders::_2), M_0) {

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