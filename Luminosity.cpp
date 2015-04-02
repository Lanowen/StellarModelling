#include "Luminosity.hpp"
#include "Star.hpp"
#include <iostream>
#include "Energy.hpp"

long double Luminosity::dL_r(long double r, long double L) {
	return 4.0L * pi * pow(r,2) * star->density.get() * star->energy.get();
}

void Luminosity::pushValues() {
	arr[0].push_back(solver.t);
	arr[1].push_back(solver.y);

	long double r = solver.t;

	long double dLpp_dr = star-> ePP ? (4.0L * pi * pow(r, 2) * star->density.get() * _e_pp) : 0.0L;
	long double dLCNO_dr = star->eCNO ? (4.0L * pi * pow(r, 2) * star->density.get() * _e_CNO) : 0.0L;
	long double dL3a_dr = star->e3a ? (4.0L * pi * pow(r, 2) * star->density.get() * _e_3alpha) : 0.0L;

	arr[2].push_back(dLpp_dr);
	arr[3].push_back(dLCNO_dr);
	arr[4].push_back(dL3a_dr);
	arr[5].push_back(dLpp_dr + dLCNO_dr + dL3a_dr);
}