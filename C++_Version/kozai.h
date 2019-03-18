#include <cstdio>
#include <cstdlib>      // stdlib.h with namespace std
#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <cmath>        // math.h   with namespace std
#include <cstring>
#include <string>
//#include <fstream>

using namespace std;

#include "vector.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>

#define sqr(x) ((x)*(x))
#define DIMENSION 19

#define PRI(x) {for (int __pri__ = 0; __pri__ < x; __pri__++) cerr << " ";}
#define PR(x)  cerr << #x << " = " << x << " "
#define PRC(x) cerr << #x << " = " << x << ",  "
#define PRL(x) cerr << #x << " = " << x << endl

#define isnan(x) ((x) != (x))

//G's and c's (with some additional powers)
//All in CGS, just to annoy the relativists
const double G=6.67407e-8; // cm^3 / g / s^2
const double c=29979245800.; // cm / s
const double Gc2=G/c/c; // cm / g
const double G3=G*G*G;
const double c5=c*c*c*c*c;
const double G3c5 = G*G*G / (c*c*c*c*c);

//Conversion Factors
const double YEAR=3.15569e7; // s
const double AU=14959787070000.0; // cm
const double RSUN=69570000000.0; // cm
const double MSUN=1.9884754153381438e+33; // g
const double DEG=0.017453292519943295;

//Some pN constants
const double c304o15=20.266666666666666;
const double c64o5=12.8;
const double c121o304=0.3980263157894737;
const double c73o24=3.0416666666666665;
const double c37o96=0.3854166666666667;

const double c29o19 = 29./19.;
const double c1181o2299 = 1181./2299.;
const double c12o19 = 12./19.;
const double c870o2299= 870./2299.;

const double SMALLNUMBER=1e-10;
const double TWOPI=6.28318530718;
const double PI=3.14159265359;

class kozai_struct;

int rhs(double t, const double y[], double f[], void *kozai_ptr);
void set_parameters(int argc, char **argv, kozai_struct *kozai, double& t_min, double& t_step, bool& IGNORE_GSL_ERROR);
void print_state(double t_next, kozai_struct *kozai, ostream& stream);
void print_header_and_initial_state(kozai_struct* kozai);
double peters_t(kozai_struct *kozai);
double peters_integral(double e, void *params);

class kozai_struct
{
	private:

		//triple configuration (initial)
		//
		//these are *not* returned by the getter
		//functions below; instead, they are
		//computed from the evolved y vector
		double a1;
		double a2;
		double e1;
		double e2;
		double inc; 
		double g1;
		double g2;
		double m1;
		double m2;
		double m3;
		double Omega1;
		double Omega2;
		double r1;
		double r2;

		//Spin information for inner binary
		double chi1;
		double chi2;
		double theta1;
		double theta2;
		double phi1;
		double phi2;

		//Things to compute ahead of time...
		double L1_no_a;
		double L2;
		double mu1;
		double mu2;
		double m1m2;
		double m2m1;
		double eta;

		//flags determining what terms to include
		bool quadrupole;
		bool octupole;
		bool pericenter;
		bool spinorbit;
		bool spinspin;
		bool radiation;

		//The initial configuration of the system in j,e,s space
		//and the initial 18 dimensional vector (allocated outside class)
		vec j1_init;
		vec j2_init;
		vec e1_init;
		vec e2_init;
		vec s1_init;
		vec s2_init;

		double *y;


	public:

		//Sloppy af constructor
		kozai_struct(double ia1=365.11437, double ia2=2865.7776, double ie1=0.32229033, double ie2=0.35653592, double iinc=93.920405, 
			  double ig1=199.31831, double ig2=237.15944, double im1=24.22645, double im2=14.999986, double im3=21.509239,
			  double iOmega1=321.97666, double iOmega2=141.97666, double ir1=0, double ir2=0,  double ichi1=1, double ichi2=1,
			  double itheta1=0, double itheta2=0, double iphi1=0, double iphi2=0,
			  bool iquad=false, bool ioct=false, bool iperi=false, bool iso=false, bool iss=false, bool rad=false)
			{	
				a1=ia1*AU;
				a2=ia2*AU;
				e1=ie1;
				e2=ie2;
				inc=iinc*DEG;
				g1=ig1*DEG;
				g2=ig2*DEG;
				m1=im1*MSUN;
				m2=im2*MSUN;
				m3=im3*MSUN;
				Omega1=iOmega1*DEG;
				Omega2=iOmega2*DEG;
				r1=ir1*RSUN;
				r2=ir2*RSUN;
				chi1=ichi1;
				chi2=ichi2;
				theta1=itheta1*DEG;
				theta2=itheta2*DEG;
				phi1=iphi1*DEG;
				phi2=iphi2*DEG;
				quadrupole=iquad;
				octupole=ioct;
				pericenter=iperi;
				spinorbit=iso;
				spinspin=iss;
				radiation=rad;
				y = new double[DIMENSION];
			}

		//destructor
		~kozai_struct() {delete y;}

		void initialize()
		{
			mu1 = m1*m2 / (m1+m2);
			mu2 = (m1+m2)*m3 / (m1+m2+m3);
			L1_no_a = mu1*sqrt(G*(m1+m2));
			L2 = mu2*sqrt(G*(m1+m2+m3)*a2);
			m1m2 = m1/m2;
			m2m1 = m2/m1;
			eta = mu1 / (m1+m2);

			//If no radius specified, set to the Schwarzschild ISCO 
			if(r1 == 0.) r1 = 2*1476.*m1/MSUN; 
			if(r2 == 0.) r2 = 2*1476.*m2/MSUN;

			//These can't be exactly 0
			if(chi1 == 0.0) chi1 = SMALLNUMBER;
			if(chi2 == 0.0) chi2 = SMALLNUMBER;
			if(theta1 == 0.0) theta1 = SMALLNUMBER;
			if(theta2 == 0.0) theta2 = SMALLNUMBER;

			//Using the angular momenta, figure out what the individual
			//inclinations are wrt the invariant plane
			//
			//NOTE: this means we're in a coordinate system where z-hat points
			//towards to total system angular momentum (before radiation
			//reaction)
			double ia = L1_no_a*sqrt(a1*(1-sqr(e1)));
			double oa = L2*sqrt(1-sqr(e2));


			double total_ang = sqrt(sqr(ia)+sqr(oa)+2*ia*oa*cos(inc));
			double inc1 = acos((ia + oa*cos(inc))/total_ang);
			double inc2 = acos((oa + ia*cos(inc))/total_ang);

			//Define the actual vectors
			j1_init = sqrt(1 - sqr(e1)) * vec(sin(inc1) * sin(Omega1), -sin(inc1) * cos(Omega1), cos(inc1));
			j2_init = sqrt(1 - sqr(e2)) * vec(sin(inc2) * sin(Omega2), -sin(inc2) * cos(Omega2), cos(inc2));
			e1_init = e1 * vec(cos(g1)*cos(Omega1) - sin(g1)*cos(inc1)*sin(Omega1),
					cos(g1)*sin(Omega1) + sin(g1)*cos(inc1)*cos(Omega1),
					sin(g1)*sin(inc1));
			e2_init = e2 * vec(cos(g2)*cos(Omega2) - sin(g2)*cos(inc2)*sin(Omega2),
					cos(g2)*sin(Omega2) + sin(g2)*cos(inc2)*cos(Omega2),
					sin(g2)*sin(inc2));

			//Also define the spin vectors (originally parallel to j1, then
			//rotate down to the specified intial orientation)
			s1_init =  m1*m1*(G/c)*chi1*j1_init / abs(j1_init);
			s2_init = m2*m2*(G/c)*chi2*j1_init / abs(j1_init);

			vec yhat = j1_init^e1_init;
			yhat /= abs(yhat);
			vec zhat = j1_init / abs(j1_init);

			//first, rotate by theta about j1 x e1
			s1_init = s1_init*cos(theta1) + (yhat^s1_init)*sin(theta1) + yhat*(yhat*s1_init)*(1-cos(theta1));
			s2_init = s2_init*cos(theta2) + (yhat^s2_init)*sin(theta2) + yhat*(yhat*s2_init)*(1-cos(theta2));

			//Then rotate about j1 by phi
			s1_init = s1_init*cos(phi1) + (zhat^s1_init)*sin(phi1) + zhat*(zhat*s1_init)*(1-cos(phi1));
			s2_init = s2_init*cos(phi2) + (zhat^s2_init)*sin(phi2) + zhat*(zhat*s2_init)*(1-cos(phi2));

			//And set the components into the y vector (which will be passed to
			//GSL)
			for(int i=0; i<3; i++){
				y[i] = j1_init[i];
				y[i+3] = e1_init[i];
				y[i+6] = j2_init[i];
				y[i+9] = e2_init[i];
				y[i+12] = s1_init[i];
				y[i+15] = s2_init[i];
			}
			y[18] = a1;
		}

		//The getters and setters for the class;
		//better to use this than the constructor
		//
		//The inner binary
		void set_a1(double a1_i) {a1=a1_i;}
		double get_a1() {return y[18];}

		void set_ecc1(double e1_i) {e1=e1_i;}
		double get_ecc1() {return abs(this->get_e1());}

		void set_g1(double g1_i) {g1=g1_i;}
		double get_g1() {
			double g1 = asin(this->get_e1()[2] / (sin(this->get_inc1())) / this->get_ecc1());
			if((this->get_e1()^this->get_j1())[2] > 0) g1 = PI - g1;
			if(g1 < 0) g1 += TWOPI; // need to put the angle in the right quadrant
			return g1;
		}

		void set_Omega1(double Omega1_i) {Omega1=Omega1_i;}
		double get_Omega1() {
			double norm_xy = sqrt(sqr(this->get_j1()[0]) + sqr(this->get_j1()[1]));
			if(norm_xy == 0.)
				return 0.;
			else{
				if(this->get_j1()[0] >= 0)
					return acos(-this->get_j1()[1] / norm_xy);
				else
		        	return TWOPI - acos(-this->get_j1()[1] / norm_xy);
			}
		}

		void set_m1(double m1_i) {m1=m1_i;}
		double get_m1() {return m1;}

		void set_m2(double m2_i) {m2=m2_i;}
		double get_m2() {return m2;}

		//Outer Binary
		void set_a2(double a2_i) {a2=a2_i;}
		double get_a2() {return a2;}

		void set_ecc2(double e2_i) {e2=e2_i;}
		double get_ecc2() {return abs(this->get_e2());}

		void set_g2(double g2_i) {g2=g2_i;}
		double get_g2() {
			double g2 = asin(this->get_e2()[2] / (sin(this->get_inc2())) / this->get_ecc2());
			if((this->get_e2()^this->get_j2())[2] > 0) g2 = PI- g2;
			if(g2 < 0) g2 += TWOPI; // need to put the angle in the right quadrant
			return g2;
		}

		void set_Omega2(double Omega2_i) {Omega2=Omega2_i;}
		double get_Omega2() {
			double norm_xy = sqrt(sqr(this->get_j2()[0]) + sqr(this->get_j2()[1]));
			if(norm_xy == 0.)
				return 0.;
			else{
				if(this->get_j2()[0] >= 0)
					return acos(-this->get_j2()[1] / norm_xy);
				else
		        	return TWOPI - acos(-this->get_j2()[1] / norm_xy);
			}
		}

		void set_m3(double m3_i) {m3=m3_i;}
		double get_m3() {return m3;}

		void set_inc(double inc_i) {inc=inc_i;}
		double get_inc() {
			return acos(this->get_j1()*this->get_j2()/(abs(this->get_j1())*abs(this->get_j2())));
		}

		double get_inc1(){
			return acos(this->get_j1()[2] / abs(this->get_j1()));
		}

		double get_inc2(){
			return acos(this->get_j2()[2] / abs(this->get_j2()));
		}

		//Physical Properties of Inner Binary (including spins)
		void set_r1(double r1_i) {r1=r1_i*RSUN;}
		double get_r1() {return r1;}

		void set_r2(double r2_i) {r2=r2_i*RSUN;}
		double get_r2() {return r2;}

		void set_chi1(double chi1_i) {chi1=chi1_i;}
		double get_chi1() {return chi1;}

		void set_chi2(double chi2_i) {chi2=chi2_i;}
		double get_chi2() {return chi2;}

		void set_theta1(double theta1_i) {theta1=theta1_i;}
		double get_theta1() {
			double dot = (this->get_j1()*this->get_s1()/(abs(this->get_j1())*abs(this->get_s1())));
			return acos(min(max(dot,-1.0),1.0)); // careful with rounding errors if theta = 0
		}

		void set_theta2(double theta2_i) {theta2=theta2_i;}
		double get_theta2() {
			double dot = (this->get_j1()*this->get_s2()/(abs(this->get_j1())*abs(this->get_s2())));
			return acos(min(max(dot,-1.0),1.0)); // careful with rounding errors if theta = 0
		}

		double get_Theta1() {
			double dot = this->get_s1()[2]/abs(this->get_s1());
			return acos(min(max(dot,-1.0),1.0)); // careful with rounding errors if theta = 0
		}

		double get_Theta2() {
			double dot = this->get_s2()[2]/abs(this->get_s2());
			return acos(min(max(dot,-1.0),1.0)); // careful with rounding errors if theta = 0
		}

		void set_phi1(double phi1_i) {phi1=phi1_i;}
		double get_phi1() {
			double dot = (this->get_e1()*this->get_s1()/(abs(this->get_e1())*abs(this->get_s1())));
			return acos(min(max(dot,-1.0),1.0)); // careful with rounding errors if phi = 0
		}

		void set_phi2(double phi2_i) {phi2=phi2_i;}
		double get_phi2() {
			double dot = (this->get_e1()*this->get_s2()/(abs(this->get_e1())*abs(this->get_s2())));
			return acos(min(max(dot,-1.0),1.0)); // careful with rounding errors if phi = 0
		}

		double get_deltaPhi(){
			vec s1h = this->get_s1() / abs(this->get_s1());
			vec s2h = this->get_s2() / abs(this->get_s2());
			vec j1h = this->get_j1() / abs(this->get_j1());
			vec first_proj = s1h-(s1h*j1h)*j1h;
			vec secon_proj = s2h-(s2h*j1h)*j1h;
			double dot = first_proj*secon_proj / (abs(first_proj)*abs(secon_proj));
			if(abs(dot) < 1.0)
				return acos(dot);
			else if(dot >= 1.0)
				return 0.;
			else
				return PI;
		}

		double get_s1s2(){
			vec s1h = this->get_s1() / abs(this->get_s1());
			vec s2h = this->get_s2() / abs(this->get_s2());
			return acos(min(max(s1h*s2h,-1.0),1.0));
		}

		double collision() {
			return (this->get_a1()*(1-this->get_ecc1()) < (this->get_r1() + this->get_r2()));
		}

		double get_forb() {
			return sqrt(G*(this->get_m1()+this->get_m2())/pow(this->get_a1(),3))/(2*M_PI);
		}

		// get the gravitational-wave frequency of the inner binary
		double gwave_freq() {
			double ecc = this->get_ecc1();
			double forb = this->get_forb();
			return 2*forb*pow(1+ecc,1.1954)/pow(1-ecc*ecc,1.5);
		}

		//Which terms to include in integration
		void set_quadrupole(bool quadrupole_i) {quadrupole=quadrupole_i;}
		bool get_quadrupole() {return quadrupole;}

		void set_octupole(bool octupole_i) {octupole=octupole_i;}
		bool get_octupole() {return octupole;}

		void set_pericenter(bool pericenter_i) {pericenter=pericenter_i;}
		bool get_pericenter() {return pericenter;}

		void set_spinorbit(bool spinorbit_i) {spinorbit=spinorbit_i;}
		bool get_spinorbit() {return spinorbit;}

		void set_spinspin(bool spinspin_i) {spinspin=spinspin_i;}
		bool get_spinspin() {return spinspin;}

		void set_radiation(bool radiation_i) {radiation=radiation_i;}
		bool get_radiation() {return radiation;}

		//precomputed constants
		//
		//except for L1, which changes as a1 shrinks
		//due to radiation reaction
		double get_L1() {return L1_no_a*sqrt(y[18]);}
		double get_L1_no_a() {return L1_no_a;}
		double get_L2() {return L2;}
		double get_mu1() {return mu1;}
		double get_mu2() {return mu2;}
		double get_m1m2() {return m1m2;}
		double get_m2m1() {return m2m1;}
		double get_eta() {return eta;}

		//Pull the vectors back out of y
		vec get_j1() {return vec(y[0],y[1],y[2]);}
		vec get_e1() {return vec(y[3],y[4],y[5]);}
		vec get_j2() {return vec(y[6],y[7],y[8]);}
		vec get_e2() {return vec(y[9],y[10],y[11]);}
		vec get_s1() {return vec(y[12],y[13],y[14]);}
		vec get_s2() {return vec(y[15],y[16],y[17]);}

		//Return the pointer to y
		double *get_y() {return y;}

};

