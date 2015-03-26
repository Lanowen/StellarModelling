#pragma once

class Mu{
public:
	long double X;
	long double Y;
	long double Z;

	long double mu;

	Mu(long double x, long double y, long double z) : X(x), Y(y), Z(z), mu(1.0 / (2.0*X + 0.75*Y + 0.5*Z)) {
	//Mu() : X(0.1), Y(0.9), Z(0.01), mu(1.0 / (2.0*X + 0.75*Y + 0.5*Z)) {
	}

	void iterate() {
		mu = 1.0 / (2.0*X + 0.75*Y + 0.5*Z);
	}

	long double get_mu() {
		return mu;
	}

	long double get_X() {
		return X;
	}

	long double get_Y() {
		return Y;
	}

	long double get_Z() {
		return Z;
	}
};

