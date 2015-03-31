#pragma once

#include "RK4.hpp"
#include "Star_Constants.hpp"
#include "Graphable.hpp"

class Star;

class Tau {
public:

	RK4 solver;
	Star* star;
	vector<vector<long double>> arr;

	long double dtau_r(long double r, long double tau);

	Tau(Star* star) : star(star), solver(bind(&Tau::dtau_r, this, placeholders::_1, placeholders::_2), tau_0) {
		arr.resize(2);
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