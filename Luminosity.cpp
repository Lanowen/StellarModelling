#include "Luminosity.hpp"
#include "Star.hpp"
#include <iostream>
#include "Energy.hpp"

long double Luminosity::dL_r(long double r, long double L) {
	return 4.0 * pi * pow(r,2) * star->density.get() * star->energy.get();
}

void Luminosity::pushValues() {
	arr[0].push_back(solver.t);
	arr[1].push_back(solver.y);

	long double r = solver.t;

	long double dLpp_dr = 4.0 * pi * pow(r, 2) * star->density.get() * _e_pp;
	long double dLCNO_dr = 4.0 * pi * pow(r, 2) * star->density.get() * _e_CNO;
	long double dL3a_dr = 4.0 * pi * pow(r, 2) * star->density.get() * _e_3alpha;

	arr[2].push_back(dLpp_dr);
	arr[3].push_back(dLCNO_dr);
	arr[4].push_back(dL3a_dr);
	arr[5].push_back(dLpp_dr + dLCNO_dr + dL3a_dr);
}