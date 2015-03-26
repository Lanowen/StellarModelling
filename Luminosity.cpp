#include "Luminosity.hpp"
#include "Star.hpp"
#include <iostream>

long double Luminosity::dL_r(long double r, long double L) {
	return 4.0 * pi * pow(r,2) * star->density.get() * star->energy.get();
}