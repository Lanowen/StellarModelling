#pragma once


#define Msun 1.989E30  //Mass of Sun
#define Rsun  6.958E8     //Radius of sun
#define Lsun  3.846E26  //Luminostity of the sun
#define Tsun  5.778E3     //Temperature at the surface of the sun
#define kb  1.3806488E-23 //Boltzmann Constant
#define sigma_B  5.670373E-8
#define me  9.10938291E-31 //mass of electron kg
#define mp  1.67262178E-27 //mass of proton kg
#define u_atomic  1.66053873E-27 //kg
#define m_H  (1.00794*u_atomic) //mass of hydrogen kg
#define eV  1.60217657E-19 //joules
#define h_  6.62606957E-34
#define G  6.67384E-11 //Gravitational Constant
#define hbar  1.054571596E-34
#define pi  3.141592654
#define c_0  299792458.0
#define q_e  1.602176462E-19
#define epsilon_0 8.854187817E-12
#define H_0  (70.0*1000.0 / 3.08567758E22)
#define GYR  (365.25*24.0*60.0*60.0*1E9)
#define AU  149597870700.0 //meters
#define pc  (206264.81 * AU)
//#define a_rad  (4.0*sigma_B / c_0)
#define a_rad 7.5657E-16

//Initial conditions
#define M_0  0
#define dM_0  LDBL_EPSILON
#define L_0  0
#define dL_0  0.0
//#define rho_0  1.622E5 //rho_c from section 3
//#define rho_0 58560
#define drho_0  0.0
//#define T_0  1.571E7  //T_c from section 3
//#define T_0 8.23E6
#define dT_0  0.0
#define tau_0  0.0
//#define kappa  1.0
#define dtau_0  (star->kappa.get()*rho_0)
#define gamma (5.0/3.0)//(5.0 / 3.0)


//Surface Boundary Conditions
#define tau_out  (2.0 / 3.0) //This is the opacity at the outside of the star,
//comes from tau(infinity) - tau(Radius of star)

// Sample Values at suface from section 3
#define Tstar  3.056E3
#define Rstar  (0.865*Rsun)
#define Mstar  (0.673*Msun)
#define Lstar  ((5.86 * 10E-2)*Lsun)

//steping interval
#define iteration_step  10000


using namespace std;