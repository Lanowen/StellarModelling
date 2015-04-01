#include "Energy.hpp"
#include "Star.hpp"
#include <iostream>


long double Energy::get() {
	long double pp = star->ePP ? _e_pp : 0.0L;
	long double CNO = star->eCNO ? _e_CNO : 0.0L;
	long double ee3a = star->e3a ?  _e_3alpha : 0.0L;

	return pp + CNO + ee3a;
}

 void Energy::pushValues(){
	 long double pp = star->ePP ? _e_pp : 0.0L;
	 long double CNO = star->eCNO ? _e_CNO : 0.0L;
	 long double ee3a = star->e3a ? _e_3alpha : 0.0L;

	 arr[0].push_back(get());
	 arr[1].push_back( pp );
	 arr[2].push_back(CNO);
	 arr[3].push_back(ee3a);	 
 }