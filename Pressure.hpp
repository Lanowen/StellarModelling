#pragma once

#include "Graphable.hpp"
#include "Star_Constants.hpp"
#include <vector>

class Star;

class Pressure : public Graphable{
public:

	vector<double> *deg, *rad, *therm;

	Star* star;

	long double get_P();

	long double dP_T(long double rho, long double T);

	long double dP_rho(long double rho, long double T);

	Pressure(Star* star) : star(star), Graphable(4) {
		deg = &arr[1];
		rad = &arr[2];
		therm = &arr[3];
	}

	long double get_dP_T();

	long double get_dP_rho();

	virtual void pushValues();
};