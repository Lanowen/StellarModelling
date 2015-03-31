#pragma once

#include "RK4.hpp"
#include "Star_Constants.hpp"
#include "Graphable.hpp"

class Star;

class Density : public Graphable  {
public:

	RK4 solver;
	Star* star;

	long double dRho_dr(long double r, long double rho);

	Density(Star* star, long double rho_c) : star(star), solver(bind(&Density::dRho_dr, this, placeholders::_1, placeholders::_2), rho_c, DBL_EPSILON) {

	}

	inline void iterate() {
		solver.iterate();
	}

	inline long double get() {
		return solver.get();
	}

	inline virtual void pushValues() {
		arr[0].push_back(solver.t);
		arr[1].push_back(solver.y);
	}

	inline virtual void popValues() {
		arr[0].pop_back();
		arr[1].pop_back();
	}
};