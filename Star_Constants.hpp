#pragma once


#define Msun 1.989E30L  //Mass of Sun
#define Rsun  6.958E8L     //Radius of sun
#define Lsun  3.846E26L  //Luminostity of the sun
#define Tsun  5.778E3L     //Temperature at the surface of the sun
#define kb  1.3806488E-23L //Boltzmann Constant
#define sigma_B  5.670373E-8L
#define me  9.10938291E-31L //mass of electron kg
#define mp  1.67262178E-27L //mass of proton kg
#define u_atomic  1.66053873E-27L //kg
#define m_H  (1.00794L*u_atomic) //mass of hydrogen kg
#define eV  1.60217657E-19L //joules
#define h_  6.62606957E-34L
#define G  6.67384E-11L //Gravitational Constant
#define hbar  1.054571596E-34L
#define pi  3.141592654L
#define c_0  299792458.0L
#define q_e  1.602176462E-19L
#define epsilon_0 8.854187817E-12L
#define H_0  (70.0L*1000.0L / 3.08567758E22L)
#define GYR  (365.25l*24.0L*60.0L*60.0L*1E9L)
#define AU  149597870700.0L //meters
#define pc  (206264.81L * AU)
//#define a_rad  (4.0*sigma_B / c_0)
#define a_rad 7.5657E-16L

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

//steping interval
#define iteration_step  10000


using namespace std;