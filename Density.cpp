#include "Density.hpp"
#include "Star.hpp"

long double Density::dRho_dr(long double r, long double rho) {
	return -(G*star->mass.get()*rho/pow(r, 2) + star->pressure.dP_T(rho, star->temperature.get())*star->temperature.dT_dr(r, star->temperature.get())) / star->pressure.dP_rho(rho, star->temperature.get());
}