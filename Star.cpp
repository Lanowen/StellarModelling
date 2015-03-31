#include "Star.hpp"
#include <sstream>
#include <iostream>

#define lineWidth 4

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
		RK4::step = max(step_min, min(RK4::step*1.8L*min(terr, 2.0L), step_max));

		//clear();

		if (lastStep != step_min && RK4::step != step_min && terr < 1E-14L) {
			T1.clear();
			T2.clear();
			continue;
		}
		break;

	}while(true);

	//cout << RK4::step << " " << abs((T2[0] - T1[0]) / T2[0]) << " " << terr << endl;

	//cout << this->temperature.arr[1].back() << endl;
	
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
	//cout << kappa.get() << " " << kappa.get()*pow(density.get(), 2) / abs(density.dRho_dr(density.arr[0].back(), density.get())) << endl;
	while (kappa.get()*pow(density.get(), 2) / abs(density.dRho_dr(density.arr[0].back(), density.get())) > 1E-30 && temperature.arr[0].back() < int_R_stop*Rsun && mass.get() < 1E3*Msun) {
		this->iterate();
		//cout << this->temperature.arr[1].back() << endl;
	}

	//for (int i = 0; i < 200000*int_R_stop; i++) {
	//	this->iterate(250);
	//}

	while (isnan(temperature.arr[1].back())) {
		popValues();
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

	if (R_star < LDBL_EPSILON) {
		R_star = this->temperature.arr[1].back();
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

void Star::graph(int starnum, bool makepdf) {

	if (!makepdf) {
		dlin.metafl("XWIN");
	}

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
		this->temperature.arr[0][i] /= R_star;
	}

	//find convective regeions;
	vector<double> conv;
	bool below = false;
	vector<double> *crit = temperature.dlnPdlnT;
	for (int i = 0; i < crit->size(); i++) {
		if (crit->at(i) < 1.0 / (1.0 - 1.0 / gamma) && !below) {
			conv.push_back(i);
			below = true;
		}
		else if (crit->at(i) > 1.0 / (1.0 - 1.0 / gamma) && below) {
			conv.push_back(i);
			below = false;
		}
	}

	vector<vector<double>> conv_x;
	vector<double> conv_y;
	vector<double> conv_end;
	conv_y.push_back(-1000.0);
	conv_y.push_back(1000.0);

	conv_end.push_back(1.0);
	conv_end.push_back(1.0);

	//conv_x.resize(conv.size() / 2);
	for (int i = 0; i < conv.size(); i += 2) {
		conv_x.resize(conv_x.size() + 1);
		conv_x.back().push_back(temperature.arr[0][conv[i]]);
		conv_x.back().push_back(temperature.arr[0][conv[i]]);

		if (i + 1 < conv.size()) { //othwise, end of star reached
			conv_x.back().push_back(temperature.arr[0][conv[i + 1]]);
			conv_x.back().push_back(temperature.arr[0][conv[i + 1]]);
		}
	}

	auto plotconv = [&]() {
		for (int i = 0; i < conv_x.size(); i++) {
			dlin.shdpat(SHADING_FILLED);
			dlin.setrgb(0.75294, 0.75294, 0.85098);

			if (conv_x[i].size() <= 2) //end of star reached
				dlin.shdcrv(&conv_x[i][0], conv_y.data(), 2, conv_end.data(), conv_y.data(), 2);
			else
				dlin.shdcrv(&conv_x[i][0], conv_y.data(), 2, &conv_x[i][2], conv_y.data(), 2);
		}

		//for (int i = 0; i < conv_x.size(); i += 2) {
		//	if (i + 1 >= conv_x.size()) //end of star reached
		//		dlin.shdcrv(&conv_x[i][0], conv_y.data(), 2, conv_end.data(), conv_y.data(), 2);
		//	else
		//		dlin.shdcrv(&conv_x[i][0], conv_y.data(), 2, &conv_x[i + 1][0], conv_y.data(), 2);
		//}
	};

	for (int i = 0; i < this->kappa.k->size(); i++) {
		this->kappa.k->at(i) = log10(this->kappa.k->at(i));
		this->kappa.kes->at(i) = log10(this->kappa.kes->at(i));
		this->kappa.kff->at(i) = log10(this->kappa.kff->at(i));
		this->kappa.kh->at(i) = log10(this->kappa.kh->at(i));
	}

	//dislinInit("$\\log(T)$$ (K)", "$\\log(\\kappa)$ ($\\frac{m^2}{kg}$) ", "Kappa");
	//dlin.graf(3, 8, 3, 1, -2, 10, -2, 1);
	//dlin.nochek();
	//dlin.title();
	//dlin.linwid(lineWidth);

	//dlin.color("FORE");
	//dlin.curve(this->kappa.logT->data(), this->kappa.k->data(), this->kappa.k->size());

	//dlin.lintyp(LINE_DOT);
	//dlin.color("BLUE");
	//dlin.curve(this->kappa.logT->data(), this->kappa.kes->data(), this->kappa.kes->size());

	//dlin.lintyp(LINE_DASHL);
	//dlin.color("GREEN");
	//dlin.curve(this->kappa.logT->data(), this->kappa.kff->data(), this->kappa.kff->size());

	//dlin.lintyp(LINE_DASH);
	//dlin.color("RED");
	//dlin.curve(this->kappa.logT->data(), this->kappa.kh->data(), this->kappa.kh->size());

	////dlin.lintyp(LINE_SOLID);
	////dlin.color("ORANGE");
	////dlin.curve(this->opal.v_t.data(), this->opal.v_y.data(), this->opal.v_t.size());

	//dlin.disfin();

	for (int i = 0; i < this->kappa.logT->size(); i++) {
		this->kappa.logT->at(i) = this->temperature.arr[0][i];
	}

	dislinInit("$\\frac{ R }{R_*}$", "$\\log(\\kappa)$ ($\\frac{m^2}{kg}$) ", "Kappa", starnum);

	dlin.graf(0, 1, 0, 0.1, -2, 10, -2, 1);
	dlin.nochek();
	dlin.title();
	dlin.linwid(lineWidth);

	dlin.legini(legendText, 4, 50);

	

	plotconv();

	dlin.color("FORE");
	dlin.leglin(legendText, "$\\kappa_{Tot}$", legendIndex++);
	dlin.curve(this->kappa.logT->data(), this->kappa.k->data(), this->kappa.k->size());

	dlin.lintyp(LINE_DOT);
	dlin.color("BLUE");
	dlin.leglin(legendText, "$\\kappa_{es}$", legendIndex++);
	dlin.curve(this->kappa.logT->data(), this->kappa.kes->data(), this->kappa.kes->size());

	dlin.lintyp(LINE_DASHL);
	dlin.color("GREEN");
	dlin.leglin(legendText, "$\\kappa_{ff}$", legendIndex++);
	dlin.curve(this->kappa.logT->data(), this->kappa.kff->data(), this->kappa.kff->size());

	dlin.lintyp(LINE_DASH);
	dlin.color("RED");
	dlin.leglin(legendText, "$\\kappa_{H^{-}}$", legendIndex++);
	dlin.curve(this->kappa.logT->data(), this->kappa.kh->data(), this->kappa.kh->size());

	dlin.color("FORE");
	dlin.legtit("Legend");
	dlin.legpos(1835, 420);
	dlin.legend(legendText, legendIndex - 1);

	//dlin.lintyp(LINE_SOLID);
	//dlin.color("ORANGE");
	//dlin.curve(this->opal.v_t.data(), this->opal.v_y.data(), this->opal.v_t.size());

	dlin.disfin();

	dislinInit("$\\frac{R}{R_*}$", "$\\rho/\\rho_{c}$ $T/T_{c}$ $M/M_{*}$ $L/L_{*}$", "Star", starnum);

	dlin.graf(0, 1, 0, 0.1, 0, 1, 0, 0.1);
	dlin.nochek();
	dlin.title();
	dlin.linwid(lineWidth);
	dlin.legini(legendText, 4, 50);

	plotconv();




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
	ss.str("");
	ss.precision(19);
	ss << "$T_{c}$ = " << T_c << " K";
	dlin.messag(ss.str().c_str(), 1500, 1100);
	ss.str("");
	ss.precision(ssize);
	ss << "$P_{c}$ = " << P_max << " Pa";
	dlin.messag(ss.str().c_str(), 1500, 1160);

	dlin.disfin();





	dislinInit("$\\frac{ R }{R_*}$", "$P/P_{c}$", "Pressure", starnum);

	//dlin.labdig(-1, "X");
	dlin.graf(0, 1, 0, 0.1, 0, 1, 0, 0.1);
	dlin.nochek();
	dlin.title();
	dlin.linwid(lineWidth);
	dlin.legini(legendText, 4, 50);

	plotconv();

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


	dislinInit("$\\frac{ R }{R_*}$", "dlnP / dlnT", "Schwarzschild Criterion", starnum);

	dlin.graf(0, 1, 0, 0.1, 0, 10, 0, 1);
	dlin.nochek();
	dlin.title();
	dlin.linwid(lineWidth);

	plotconv();



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



	dislinInit("$\\frac{ R }{R_*}$", "$dL/dr$ $(L_*/R_*)$", "Luminosity", starnum);

	//dlin.labdig(-1, "X");
	dlin.graf(0, 1, 0, 0.1, 0, 15, 0, 1);
	dlin.nochek();
	dlin.title();
	dlin.linwid(lineWidth);
	dlin.legini(legendText, 4, 50);

	plotconv();

	//double L_3a_max = *max_element(this->luminosity.arr[4].begin(), this->luminosity.arr[4].end());

	for (int i = 0; i < luminosity.arr[0].size(); i++) {
		luminosity.arr[5][i] /= L_max / R_star;
		luminosity.arr[2][i] /= L_max / R_star;
		luminosity.arr[3][i] /= L_max / R_star;
		luminosity.arr[4][i] /= L_max / R_star;
	}

	dlin.color("FORE");
	dlin.leglin(legendText, "$dL/dr$", legendIndex++);
	dlin.curve(this->temperature.arr[0].data(), this->luminosity.arr[5].data(), this->temperature.arr[0].size());
	dlin.lintyp(LINE_DASH);
	dlin.color("RED");
	dlin.leglin(legendText, "$dL_{pp}/dr$", legendIndex++);
	dlin.curve(this->temperature.arr[0].data(), this->luminosity.arr[2].data(), this->temperature.arr[0].size());
	dlin.lintyp(LINE_DASHM);
	dlin.color("BLUE");
	dlin.leglin(legendText, "$dL_{CNO}/dr$", legendIndex++);
	dlin.curve(this->temperature.arr[0].data(), this->luminosity.arr[3].data(), this->temperature.arr[0].size());
	dlin.lintyp(LINE_DASHL);
	dlin.color("GREEN");
	dlin.leglin(legendText, "$dL_{3\\alpha}/dr$", legendIndex++);
	dlin.curve(this->temperature.arr[0].data(), this->luminosity.arr[4].data(), this->temperature.arr[0].size());

	dlin.color("FORE");
	dlin.legtit("Legend");
	dlin.legpos(1930, 420);
	dlin.legend(legendText, legendIndex - 1);
	dlin.disfin();
}