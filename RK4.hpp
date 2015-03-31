#pragma once

#include <functional>
#include "Star_Constants.hpp"
#include <vector>

class RK4 {
public:
	vector<long double> saved_y, saved_t;
	long double y, intermed_y;
	long double k1, k2, k3, k4;
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
			intermed_y = y + k1 / 2.0;
			break;
		case 3:
			k2 = f(t + step / 2.0, intermed_y)*step;
			break;
		case 4:
			intermed_y = y + k2 / 2.0;
			break;
		case 5:
			k3 = f(t + step / 2.0, intermed_y) * step;
			break;
		case 6:
			intermed_y = y + k3;
			break;
		case 7:
			k4 = f(t + step, intermed_y)*step;
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
		y += (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
		t += step;
		//v_t.push_back(t);
		//v_y.push_back(y);
	}

	inline long double get() {
		if (updatingK)
			return intermed_y;

		return y;
	}

};