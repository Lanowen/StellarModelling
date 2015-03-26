#include "Mass.hpp"
#include "Star.hpp"

long double Mass::dM_dr(long double r, long double M) {
	return 4.0*pi*pow(r,2)*star->density.get();
}