#include "Mu.hpp"
#include "Star.hpp"


long double Mu::get_X() {
	return star->mass.get() / Msun < He_cutoff ? 0.0L : X;
}

long double Mu::get_Y() {
	return star->mass.get() / Msun < He_cutoff ? 1.0L : Y;
}

long double Mu::get_Z() {
	return star->mass.get() / Msun < He_cutoff ? 0.0L : Z;
}