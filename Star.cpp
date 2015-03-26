#include "Star.hpp"
#include <sstream>

#define lineWidth 4

void Star::step() {
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

void Star::iterate() {
	vector<double long> T1, T2;

	do {
		push();
		step();
		T1.push_back(temperature.solver.y);
		T1.push_back(density.solver.y);
		T1.push_back(luminosity.solver.y);
		T1.push_back(mass.solver.y);
		pop();

		//push();
		RK4::step /= 2;
		step();
		step();
		T2.push_back(temperature.solver.y);
		T2.push_back(density.solver.y);
		T2.push_back(luminosity.solver.y);
		T2.push_back(mass.solver.y);

		long double err, terr = LDBL_MAX;

		for (int i = 0; i < T1.size(); i++) {
			err = 1E-14L / abs((T2[i] - T1[i]) / T2[i]);
			terr = (err < terr) ? err : terr;
		}
		

		long double lastStep = RK4::step;
		RK4::step = max(1000.0L, min(RK4::step*1.8L*min(terr, 2.0L), 20000.0L));

		//clear();

		if (lastStep != 1000.0 && RK4::step != 1000.0 && terr < 1E-3L) {
			T1.clear();
			T2.clear();
			continue;
		}
		break;

	}while(true);

	//cout << RK4::step << " " << abs((T2[0] - T1[0]) / T2[0]) << " " << terr << endl;

	pushValues();
	
}

void Star::solve() {
	dlin.metafl("PDF");
	dlin.imgfmt("RGB");
	dlin.filmod("COUNT");
	dlin.scrmod("REVERS");
#ifdef _DEBUG
	dlin.metafl("XWIN");
#endif

	try {
		//while (this->density.get() > 0.0001 && this->temperature.get() > 0 && this->temperature.get() < 1E10) {
		while (this->temperature.arr[0].back() < 1.2 * Rsun) {
			this->iterate();
		}

		//RK4::step = 1000;
		for (int i = 0; i < 1000; i++) {
			this->iterate();
		}
	}
	catch (OpalException oe) {
		cout << oe.what() << endl;
	}

	
	long double Tau_inf = this->tau.arr[1].back();
	for (int i = this->tau.arr[1].size() - 1; i >= 0; i--) {
		//cout << Tau_inf - this->tau.arr[1][i] << endl;
		if (Tau_inf - this->tau.arr[1][i] > 2.0 / 3.0) {
			T_star = this->temperature.arr[1][i + 1];
			R_star = this->temperature.arr[0][i + 1];
			break;
		}
	}

	T_max = *max_element(this->temperature.arr[1].begin(), this->temperature.arr[1].end());
	L_max = *max_element(this->luminosity.arr[1].begin(), this->luminosity.arr[1].end());
	M_max = *max_element(this->mass.arr[1].begin(), this->mass.arr[1].end());
	Rho_max = *max_element(this->density.arr[1].begin(), this->density.arr[1].end());
	P_max = *max_element(this->pressure.arr[0].begin(), this->pressure.arr[0].end());

	//R_max = *max_element(this->temperature.arr[0].begin(), this->temperature.arr[0].end());

	cout << "Temperature = " << T_star << " @r= " << R_star / Rsun << " Rsun" << endl;
	cout << "Mass: " << M_max / Msun << " Msun" << endl;
	cout << "Luminosity: " << L_max / Lsun << " Lsun" << endl;

	cout << L_max << " - " << 4.0*pi*sigma_B*R_star*R_star*pow(T_star, 4) << " : " << L_max - 4.0*pi*sigma_B*R_star*R_star*pow(T_star, 4) << ". fractional: " << (L_max - 4.0*pi*sigma_B*R_star*R_star*pow(T_star, 4)) / sqrt(4.0*pi*sigma_B*R_star*R_star*pow(T_star, 4)) << endl;

}

void Star::graph() {
	for (int i = 0; i < this->kappa.k->size(); i++) {
		this->kappa.k->at(i) = log10(this->kappa.k->at(i));
		this->kappa.kes->at(i) = log10(this->kappa.kes->at(i));
		this->kappa.kff->at(i) = log10(this->kappa.kff->at(i));
		this->kappa.kh->at(i) = log10(this->kappa.kh->at(i));
	}

	dislinInit("$\\log(T)$$ (K)", "$\\log(\\kappa)$ ($\\frac{m^2}{kg}$) ", "Kappa");
	dlin.graf(3, 8, 3, 1, -2, 10, -2, 1);
	dlin.nochek();
	dlin.title();
	dlin.linwid(lineWidth);

	dlin.color("FORE");
	dlin.curve(this->kappa.logT->data(), this->kappa.k->data(), this->kappa.k->size());

	dlin.lintyp(LINE_DOT);
	dlin.color("BLUE");
	dlin.curve(this->kappa.logT->data(), this->kappa.kes->data(), this->kappa.kes->size());

	dlin.lintyp(LINE_DASHL);
	dlin.color("GREEN");
	dlin.curve(this->kappa.logT->data(), this->kappa.kff->data(), this->kappa.kff->size());

	dlin.lintyp(LINE_DASH);
	dlin.color("RED");
	dlin.curve(this->kappa.logT->data(), this->kappa.kh->data(), this->kappa.kh->size());

	//dlin.lintyp(LINE_SOLID);
	//dlin.color("ORANGE");
	//dlin.curve(this->opal.v_t.data(), this->opal.v_y.data(), this->opal.v_t.size());

	dlin.disfin();

	for (int i = 0; i < this->kappa.logT->size(); i++) {
		this->kappa.logT->at(i) = this->temperature.arr[0][i] / Rsun;
	}

	dislinInit("$\\frac{ R }{R_\\odot}$", "$\\log(\\kappa)$ ($\\frac{m^2}{kg}$) ", "Kappa");

	dlin.graf(0, 1.2, 0, 0.1, -2, 10, -2, 1);
	dlin.nochek();
	dlin.title();
	dlin.linwid(lineWidth);

	dlin.color("FORE");
	dlin.curve(this->kappa.logT->data(), this->kappa.k->data(), this->kappa.k->size());

	dlin.lintyp(LINE_DOT);
	dlin.color("BLUE");
	dlin.curve(this->kappa.logT->data(), this->kappa.kes->data(), this->kappa.kes->size());

	dlin.lintyp(LINE_DASHL);
	dlin.color("GREEN");
	dlin.curve(this->kappa.logT->data(), this->kappa.kff->data(), this->kappa.kff->size());

	dlin.lintyp(LINE_DASH);
	dlin.color("RED");
	dlin.curve(this->kappa.logT->data(), this->kappa.kh->data(), this->kappa.kh->size());

	//dlin.lintyp(LINE_SOLID);
	//dlin.color("ORANGE");
	//dlin.curve(this->opal.v_t.data(), this->opal.v_y.data(), this->opal.v_t.size());

	dlin.disfin();




	dislinInit("$\\frac{R}{R_\\odot}$", "$\\rho/\\rho_{max}$ $T/T_{max}$ $M/M_{max}$ $L/L_{max}$", "Star");

	dlin.graf(0, 1.2, 0, 0.1, 0, 1, 0, 0.1);
	dlin.nochek();
	dlin.title();
	dlin.linwid(lineWidth);
	dlin.legini(legendText, 4, 50);

	if (this->temperature.arr[1].size() > this->luminosity.arr[1].size()) {
		this->temperature.arr[1].pop_back();
		this->temperature.arr[0].pop_back();
	}
	

	for (int i = 0; i < this->temperature.arr[0].size(); i++) {
		this->temperature.arr[1][i] /= T_max;
		this->luminosity.arr[1][i] /= L_max;
		this->mass.arr[1][i] /= M_max;
		this->density.arr[1][i] /= Rho_max;
		this->pressure.arr[0][i] /= P_max;
		this->pressure.deg->at(i) /= P_max;
		this->pressure.therm->at(i) /= P_max;
		this->pressure.rad->at(i) /= P_max;
		this->temperature.dlnPdlnT[i] = (this->temperature.dlnPdlnT[i]);
		this->temperature.arr[0][i] /= Rsun;
	}

	dlin.lintyp(LINE_DASH);
	dlin.color("RED");
	dlin.leglin(legendText, "Temperature", legendIndex++);
	dlin.curve(this->temperature.arr[0].data(), this->temperature.arr[1].data(), this->temperature.arr[0].size());

	dlin.lintyp(LINE_DOT);
	dlin.color("BLUE");
	dlin.leglin(legendText, "Luminosity", legendIndex++);
	dlin.curve(this->temperature.arr[0].data(), this->luminosity.arr[1].data(), this->temperature.arr[0].size());

	dlin.lintyp(LINE_DASHL);
	dlin.color("GREEN");
	dlin.leglin(legendText, "Mass", legendIndex++);
	dlin.curve(this->temperature.arr[0].data(), this->mass.arr[1].data(), this->temperature.arr[0].size());

	//dlin.color("ORANGE");
	//dlin.curve(this->temperature.arr[0].data(), this->tau.arr[1].data(), this->temperature.arr[0].size());

	dlin.lintyp(LINE_SOLID);
	dlin.color("FORE");
	dlin.leglin(legendText, "Density", legendIndex++);
	dlin.curve(this->temperature.arr[0].data(), this->density.arr[1].data(), this->temperature.arr[0].size());

	dlin.color("FORE");
	dlin.legtit("Legend");
	dlin.legpos(1835, 420);
	dlin.legend(legendText, legendIndex - 1);

	stringstream ss;
	ss << "$T_{*}$ = " << T_star << " K";
	dlin.messag(ss.str().c_str(), 1500, 800);
	ss.str("");
	streamsize ssize = ss.precision();
	ss.precision(19);
	ss << "$\\rho_c$ = " << rho_c << " $kg/m^3$";
	dlin.messag(ss.str().c_str(), 1500, 860);
	ss.precision(ssize);
	ss.str("");
	ss << "$R_{*}$ = " << R_star / Rsun << " $R_\\odot$";
	dlin.messag(ss.str().c_str(), 1500, 920);
	ss.str("");
	ss << "$M_{*}$ = " << M_max / Msun << " $M_\\odot$";
	dlin.messag(ss.str().c_str(), 1500, 980);
	ss.str("");
	ss << "$L_{*}$ = " << L_max / Lsun << " $L_\\odot$";
	dlin.messag(ss.str().c_str(), 1500, 1040);

	dlin.disfin();

	



	dislinInit("$\\frac{ R }{R_\\odot}$", "$P/P_{Tot_{max}}$", "Pressure");

	//dlin.labdig(-1, "X");
	dlin.graf(0, 1.2, 0, 0.1, 0, 1, 0, 0.1);
	dlin.nochek();
	dlin.title();
	dlin.linwid(lineWidth);
	dlin.legini(legendText, 4, 50);

	dlin.color("FORE");
	dlin.leglin(legendText, "$P_{Tot}$", legendIndex++);
	dlin.curve(this->temperature.arr[0].data(), this->pressure.arr[0].data(), this->temperature.arr[0].size());
	dlin.lintyp(LINE_DASH);
	dlin.color("RED");
	dlin.leglin(legendText, "$P_{Deg}$", legendIndex++);
	dlin.curve(this->temperature.arr[0].data(), this->pressure.deg->data(), this->temperature.arr[0].size());
	dlin.lintyp(LINE_DASHL);
	dlin.color("GREEN");
	dlin.leglin(legendText, "$P_{Thm}$", legendIndex++);
	dlin.curve(this->temperature.arr[0].data(), this->pressure.therm->data(), this->temperature.arr[0].size());
	dlin.lintyp(LINE_DOT);
	dlin.color("BLUE");
	dlin.leglin(legendText, "$P_{Rad}$", legendIndex++);
	dlin.curve(this->temperature.arr[0].data(), this->pressure.rad->data(), this->temperature.arr[0].size());

	dlin.color("FORE");
	dlin.legtit("Legend");
	dlin.legpos(2100, 420);
	dlin.legend(legendText, legendIndex - 1);
	dlin.disfin();

	dislinInit("$\\frac{ R }{R_\\odot}$", "dlnP / dlnT", "Schwarzschild Criterion, kappa fitted, kappa opal");

	dlin.graf(0, 1.2, 0, 0.1, 0, 10, 0, 1);
	dlin.nochek();
	dlin.title();
	dlin.linwid(lineWidth);
	dlin.lintyp(LINE_SOLID);
	dlin.color("FORE");
	dlin.curve(this->temperature.arr[0].data(), this->temperature.dlnPdlnT->data(), this->temperature.arr[0].size());

	//dlin.lintyp(LINE_SOLID);
	//dlin.color("orange");
	//dlin.curve(this->temperature.arr[0].data(), this->temperature.dlnPdlnT_opal.data(), this->temperature.arr[0].size());

	dlin.lintyp(LINE_DOT);
	dlin.color("RED");
	dlin.rline(0, 1.0 / (1.0 - 1.0 / gamma), 1.2, 1.0 / (1.0 - 1.0 / gamma));

	dlin.messag("1/(1-$\\frac{1}{\\gamma})$", 2500, 1450);

	dlin.disfin();
}