#pragma once

#include <vector>

using namespace std;

class Graphable{
public:
	vector<vector<double>> arr;

	Graphable(int size = 2) {
		arr.resize(size);
	}

	virtual void pushValues() = 0;
};