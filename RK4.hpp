#pragma once

#include <functional>
#include "Star_Constants.hpp"
#include <vector>

class RK4 {
public:
	vector<long double> saved_y, saved_t;
	long double y, intermed_y;
	long double k1, k2, k3, k4, k5, k6;
	long double t;
	static long double step;
	function<long double(long double, long double)> f;

	unsigned int currK;
	bool updatingK;
	
	inline bool updateK() {
		switch (currK) {
		case 1:
			updatingK = true;
			intermed_y = y;
			k1 = f(t, y)*step;
			break;
		case 2:
			intermed_y = y + k1 / 4.0;
			break;
		case 3:
			k2 = f(t + step / 4.0, intermed_y)*step;
			break;
		case 4:
			intermed_y = y + k1 * 3.0 / 32.0 + k2 * 9.0 / 32.0;
			break;
		case 5:
			k3 = f(t + step * 3.0 / 8.0, intermed_y) * step;
			break;
		case 6:
			intermed_y = y  + k1*1932.0/2197.0 - k2*7200.0/2197.0 + k3*7296.0/2197.0;
			break;
		case 7:
			k4 = f(t + step*12.0/13.0, intermed_y)*step;
			break;
		case 8:
			intermed_y = y + k1*439.0 / 216.0 - k2*8.0 + k3*3680.0 / 513.0 - k4*845.0 / 4104.0;
			break;
		case 9:
			k5 = f(t + step, intermed_y)*step;
			break;
		case 10:
			intermed_y = y - k1*8.0 / 27.0 + k2*2.0 - k3*3544.0 / 2565.0 + k4*1859.0 / 4104.0 - k5*11.0 / 40.0;
			break;
		case 11:
			k6 = f(t + step / 2.0, intermed_y)*step;
			break;
		default:
			updatingK = false;
			currK = 1;
			return false;
		}

		currK++;

		return true;
	}

	inline void push() {
		saved_y.push_back(y);
		saved_t.push_back(t);
	}

	inline void pop() {
		y = saved_y.back();
		t = saved_t.back();
		saved_y.pop_back();
		saved_t.pop_back();
	}

	inline void clear() {
		saved_y.clear();
		saved_t.clear();
	}

	RK4(function<long double(long double, long double)> func_derivative, long double initial_value, long double initial_t = 0) : t(initial_t), f(func_derivative), y(initial_value), currK(1), updatingK(false) {
	}

	inline void iterate() {
		y += k1*25.0 / 216.0 + k3*1408.0 / 2565.0 + k4*2197.0 / 4104.0 - k5 / 5.0;
		t += step;
		//v_t.push_back(t);
		//v_y.push_back(y);
	}

	inline long double getError() {
		return abs((k1*16.0/135.0 + k3*6656.0/12825.0 + k4*28561.0/56430.0 - k5*9.0/50.0 + k6*2.0/55.0) - (k1*25.0 / 216.0 + k3*1408.0 / 2565.0 + k4*2197.0 / 4104.0 - k5 / 5.0));
	}

	inline long double getRelError() {
		long double yi = k1*25.0 / 216.0 + k3*1408.0 / 2565.0 + k4*2197.0 / 4104.0 - k5 / 5.0;
		long double zi = k1*16.0 / 135.0 + k3*6656.0 / 12825.0 + k4*28561.0 / 56430.0 - k5*9.0 / 50.0 + k6*2.0 / 55.0;

		return abs((zi - yi) / yi);
	}

	inline long double get() {
		if (updatingK)
			return intermed_y;

		return y;
	}

};