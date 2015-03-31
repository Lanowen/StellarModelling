#pragma once

#include "RK4.hpp"
#include <algorithm>
#include <functional>
#include "Star_Constants.hpp"
#include <vector>
#include "Graphable.hpp"

class Star;

class Temperature : public Graphable {
public:
	RK4 solver;
	vector<double> *dlnPdlnT;// , *dlnPdlnT_opal;
	Star* star;
	

	int test_switch = 0;

	long double dT_dr_conv(long double r, long double T);

	long double dT_dr_rad(long double r, long double T);

	long double dT_dr(long double r, long double T);

	long double criterion();

	//long double get_gamma();
 
	Temperature(Star* star, long double T_c) : star(star), solver(bind(&Temperature::dT_dr, this, placeholders::_1, placeholders::_2), T_c, DBL_EPSILON), Graphable(3) {
		dlnPdlnT = &arr[2];
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
		dlnPdlnT->push_back(1.0 / criterion());
		//dlnPdlnT_opal->push_back(1.0 / ((3.0*star->opal.get_kappa(get(), star->density.get())*star->luminosity.get()*star->pressure.get_P()) / (pow(get(), 4) * 16.0 * a_rad * c_0 * pi *G * star->mass.get())));
	}

	inline virtual void popValues() {
		arr[0].pop_back();
		arr[1].pop_back();
		dlnPdlnT->pop_back();
		//dlnPdlnT_opal->push_back(1.0 / ((3.0*star->opal.get_kappa(get(), star->density.get())*star->luminosity.get()*star->pressure.get_P()) / (pow(get(), 4) * 16.0 * a_rad * c_0 * pi *G * star->mass.get())));
	}
};