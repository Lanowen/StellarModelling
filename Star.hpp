#pragma once

#include "Density.hpp"
#include "Energy.hpp"
#include "Kappa.hpp"
#include "Luminosity.hpp"
#include "Mu.hpp"
#include "Mass.hpp"
#include "Tau.hpp"
#include "Pressure.hpp"
#include "Temperature.hpp"
#include "Opal.hpp"

#include <discpp.h>
#include <sstream>

using namespace std;

static Dislin dlin;

static char legendText[255] = { '/0' };
static int legendIndex = 1;

static void dislinInit(char* xName, char* yName, char* title, int i) {
	stringstream ss;
	ss << "star" << i << ".pdf";
	dlin.setfil(ss.str().c_str());
	dlin.disini();
	//dlin.nochek();
	dlin.pagera();
	dlin.name(xName, "X");
	dlin.name(yName, "Y");
	dlin.titlin(title, 1);

	legendIndex = 1;
	dlin.complx();
	dlin.texmod("ON");
};

class Star {
public:

	Density density;
	Energy energy;
	Kappa kappa;
	Luminosity luminosity;
	Mu mu;
	Mass mass;
	Tau tau;
	Pressure pressure;
	Temperature temperature;
	Opal opal;

	long double step_min;
	long double step_max;
	long double int_R_stop;

	long double R_star, T_star;
	long double T_c, rho_c;
	long double T_max, L_max, M_max, Rho_max,P_max;
	bool ePP, eCNO, e3a;
	bool use_opal;
	long double err_sensitivity;

	Star(long double T_c, long double rho_c, long double x, long double y, long double z, long double step_min, long double step_max, long double int_R_stop,
		 bool ePP, bool eCNO, bool e3a, bool use_opal, long double err_sensitivity, string opal_table, long double He_cutoff = 0.0L) :
		 use_opal(use_opal), opal(this, opal_table, use_opal), density(this, rho_c), energy(this), kappa(this), luminosity(this), mu(this, x, y, z, He_cutoff), mass(this), tau(this), pressure(this), temperature(this, T_c),
		R_star(0), T_star(0), T_c(T_c), rho_c(rho_c), step_min(step_min), step_max(step_max), int_R_stop(int_R_stop), ePP(ePP), eCNO(eCNO), e3a(e3a), err_sensitivity(err_sensitivity)
	{
		//add plot points for values
		pushValues();
	}

	void solve();

	void graph(int starnum, bool makepdf = true);

	inline void push() {
		temperature.solver.push();
		density.solver.push();
		luminosity.solver.push();
		mass.solver.push();
	}

	inline void pop() {
		temperature.solver.pop();
		density.solver.pop();
		luminosity.solver.pop();
		mass.solver.pop();
	}

	inline void clear() {
		temperature.solver.clear();
		density.solver.clear();
		luminosity.solver.clear();
		mass.solver.clear();
	}

	inline void pushValues() {
		temperature.pushValues();
		density.pushValues();
		luminosity.pushValues();
		mass.pushValues();
		tau.pushValues();
		pressure.pushValues();
		kappa.pushValues();
		opal.pushValues();
	}

	inline void popValues() {
		temperature.popValues();
		density.popValues();
		luminosity.popValues();
		mass.popValues();
		tau.popValues();
		pressure.popValues();
		kappa.popValues();
		opal.popValues();
	}

	inline long double frac_diff() {
		return (L_max - 4.0*pi*sigma_B*R_star*R_star*pow(T_star, 4)) / sqrt(4.0*pi*sigma_B*R_star*R_star*pow(T_star, 4)*L_max);
	}

	inline void step() {
		//calculate k values before iteration
		do {
			temperature.solver.updateK();
			density.solver.updateK();
			luminosity.solver.updateK();
			mass.solver.updateK();
		} while (tau.solver.updateK());

		//iterate rk4
		temperature.iterate();
		density.iterate();
		luminosity.iterate();
		mass.iterate();
		tau.iterate();
	}

	inline void iterate();

	inline void iterate(long double rk_step) {
		RK4::step = rk_step;
		step();

		pushValues();
	}
};