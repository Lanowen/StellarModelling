#pragma once
#include "Star_Constants.hpp"

class Star;

class Mu{
public:
	long double X;
	long double Y;
	long double Z;

	long double He_cutoff;

	Star* star;

	Mu(Star* star, long double x, long double y, long double z, long double He_cutoff = 0.0L) : star(star), X(x), Y(y), Z(z), He_cutoff(He_cutoff) {
	//Mu() : X(0.1), Y(0.9), Z(0.01), mu(1.0 / (2.0*X + 0.75*Y + 0.5*Z)) {
	}

	inline long double get_mu() {
		return 1.0 / (2.0*get_X() + 0.75*get_Y() + 0.5*get_Z());
	}

	long double get_X();

	long double get_Y();

	long double get_Z();
};

