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

	Luminosity(Star* star) : star(star), solver(bind(&Luminosity::dL_r, this, placeholders::_1, placeholders::_2), L_0), Graphable(6) {

	}

	inline void iterate() {
		solver.iterate();
	}

	inline long double get() {
		return solver.get();
	}

	virtual void pushValues();

	inline virtual void popValues() {
		arr[0].pop_back();
		arr[1].pop_back();
		arr[2].pop_back();
		arr[3].pop_back();
		arr[4].pop_back();
		arr[5].pop_back();
	}
};