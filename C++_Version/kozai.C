#include "kozai.h"

const string help="\nCarl's Kozai Code\n\n\
For integrating the secular equations of motions\n\
for the three-body problem, including the quadrupole\n\
and octupole terms (see Tremaine+2009 and Liu+2015)\n\n\
Also includes relativistic terms for the inner binary, including:\n\
  Pericenter precession (ala Ford+2000)\n\
  Gravitational Wave emission (ala Blaes+2002)\n\
  Spin-Orbit coupling (ala me, following Barker and O'Connell 1975)\n\
  Spin-Spin coupling (same)\n\n\
Options:\n\
  --no-int  don't integrate the system, just print info\n\
  --help    print this\n\n\
Parameters:\n\
  --m1      mass 1 of inner binary (MSUN)\n\
  --m2      mass 2 of inner binary (MSUN)\n\
  --m3      mass of tertiary (MSUN)\n\
  --a1      Semi-major axis of inner binary (AU)\n\
  --a2      Semi-major axis of outer binary (AU)\n\
  --e1      Eccentricity of inner binary\n\
  --e2      Eccentricity of outer binary\n\
  --g1      Argument of pericenter for inner binary (deg)\n\
  --g2      Argument of pericenter for outer binary (deg)\n\
  --omega1  Longitude of ascending node for inner binary (deg)\n\
  --omega2  Longitude of ascending node for outer binary (deg)\n\
  --inc     Mutual Inclination (deg)\n\
  --rad1    Radius of mass 1 (RSUN)\n\
  --rad2    Radius of mass 2 (RSUN)\n\
  --chi1    Dimmensionless Spin of mass 1 (0-1)\n\
  --chi2    Dimmensionless Spin of mass 2 (0-1)\n\
  --theta1  Spin-orbit misalignment of spin 1 (deg)\n\
  --theta2  Spin-orbit misalignment of spin 2 (deg)\n\
  --phi1    Angle in orbital plane between pericenter and spin 1 (deg)\n\
  --phi2    Angle in orbital plane between pericenter and spin 2 (deg)\n\
  --time    Time to integrate (years)\n\
  --dt      How often to print output\n\
                Positive dt outputs every dt years\n\
                Negative dt outputs every 10^{-dt} adaptive timesteps\n\
                    (e.g. dt=-2 outputs every 100 timesteps)\n\n\
Flags:\n\
  --quad       Include quadrupole terms \n\
  --oct        Include octupole Terms\n\
  --peri       Include pericenter precession (1pN) terms\n\
  --spinorbit  Include spin-orbit (1.5pN) terms\n\
  --spinspin   Include spin-spin (2pN) terms\n\
  --rad        Include gravitational-wave emission (2.5pN)\n\
  --ignore_gsl Ignore errors from GSL (i.e. when a binary merges, the inner\n\
               binary can decouple from the outer s.t. the optimal timestep\n\
               for the inner is beyond machine-tolerance away from the outer)\n";

// Function for printing the state of a triple; placed here for easy
// modification
void print_state(double t_next, kozai_struct *kozai, ostream& stream = cerr){

	streamsize ss = stream.precision();

	stream << setprecision(14) << (t_next/YEAR) << "  " 
		 << kozai->get_a1()/AU << "  " 
		<< kozai->get_ecc1() << "  "
		<< setprecision(ss)
		<< kozai->get_ecc2() << "  " 
		<< kozai->get_inc()/DEG << "  " 
		<< kozai->get_inc1()/DEG << "  " 
		<< kozai->get_inc2()/DEG << "  " 
		<< kozai->get_g1()/DEG << "  " 
		<< kozai->get_g2()/DEG << "  "
		<< kozai->get_theta1()/DEG << "  " 
		<< kozai->get_theta2()/DEG << "  " 
		<< kozai->get_deltaPhi()/DEG << "  " 
		<< kozai->get_Theta1()/DEG << "  " 
		<< kozai->get_Theta2()/DEG << "  " 
		<< kozai->get_s1s2()/DEG << endl; 
}

int main(int argc, char **argv){

	double eps_abs = 1.e-12;
	double eps_rel = 1.e-12;
	double h = 1.e-10;
	double tmin=0, tmax=3.2e8*YEAR, delta_t=1e5*YEAR;
	bool collision=false;
	bool failure=false;
	bool IGNORE_GSL_ERRORS=false;
	bool TENhz_ecc=false;
	bool TENhz_circ=false;
	double t=tmin;
	double t_insp;
	double max_e=0.;
	int status;


	//First, initialize the Kozai structure
	kozai_struct *kozai = new kozai_struct;
	set_parameters(argc, argv, kozai, tmax, delta_t, IGNORE_GSL_ERRORS);
	kozai->initialize();

	//Then, set up the GSL integrator (we're using an 8th-order Runge Kutta)
	const gsl_odeiv2_step_type *type_ptr = gsl_odeiv2_step_rk8pd;

	gsl_odeiv2_step *step_ptr = gsl_odeiv2_step_alloc(type_ptr, DIMENSION);
	gsl_odeiv2_control *control_ptr = gsl_odeiv2_control_y_new(eps_abs, eps_rel);
	gsl_odeiv2_evolve *evolve_ptr = gsl_odeiv2_evolve_alloc(DIMENSION);
	gsl_odeiv2_system my_system;

	my_system.function = rhs;
	my_system.dimension = DIMENSION;
	my_system.params = kozai;

	print_header_and_initial_state(kozai);
	print_state(0.,kozai);


	//Then actually integrate the system forward over each timestep
	long step = 0;
	long delta_output = 1;
	bool output_every_nsteps = false;

	for (double t_next = delta_t; t_next < tmax; t_next += delta_t){

		// if printing every x timesteps, the while loop does the whole
		// integration
		if (delta_t <= 0){
			delta_output = pow(10,-delta_t);
			delta_t = 1e100;
			t_next = tmax;
			output_every_nsteps = true;
		}

		//Then integrate the system forward for this timestep
		while (t < t_next){
			status = gsl_odeiv2_evolve_apply (evolve_ptr, control_ptr, step_ptr,
                                    &my_system, &t, t_next, &h, kozai->get_y());


			//Check for collisions
			if(kozai->collision()){
				collision = true;
				break;
			}

			//Record the maximum ecc of the inner binary
			if (kozai->get_ecc1() > max_e)
				max_e = kozai->get_ecc1();
			
			//check for hitting the LIGO band
			if(kozai->gwave_freq() > 10 && TENhz_ecc == false){
				TENhz_ecc = true;
				cout << "LIGO band (eccentric): ";
				print_state(t,kozai,cout);
			}
			
			if(2*kozai->get_forb() > 10){
				TENhz_circ = true;
				cout << "LIGO band (circular): ";
				print_state(t,kozai,cout);
				break;
			}


			//check for some type of failure/NaN
			if ((status != GSL_SUCCESS && !IGNORE_GSL_ERRORS) || isnan(kozai->get_ecc1()) || kozai->get_ecc1() >= 1.){
				failure = true;
				break;
			}

			//otherwise, print the step (if we're printing on the integrator
			//timescale)
			if (output_every_nsteps && step % delta_output == 0)
				print_state(t,kozai);
			step += 1;
		}

		//Finally, print the state of the system
		print_state(t,kozai);

		if (collision == true){
			cout << "Collision!" << endl;
			break;
		}

		if (TENhz_circ == true){
			cout << "LIGO band!" << endl;
			break;
		}

		if (failure == true){
			cout << "NaNs (if near merger, try --ignore_gsl)!" << endl;
			break;
		}
	}

	streamsize ss = cout.precision();
	cout << setprecision(10);
	cout << "Maximum Eccentricity: " << max_e << endl;
	cout << setprecision(ss);

	/* all done; free up the gsl_odeiv stuff */
	gsl_odeiv2_evolve_free (evolve_ptr);
	gsl_odeiv2_control_free (control_ptr);
	gsl_odeiv2_step_free (step_ptr);
	delete kozai;
}

void print_header_and_initial_state(kozai_struct* kozai){

	cout << "Triple setup:\n"  
		<< "  m1 = " << kozai->get_m1()/MSUN << endl
		<< "  m2 = " << kozai->get_m2()/MSUN << endl
		<< "  m3 = " << kozai->get_m3()/MSUN << endl
		<< "  a1 = " << kozai->get_a1()/AU << endl
		<< "  a2 = " << kozai->get_a2()/AU << endl
		<< "  e1 = " << kozai->get_ecc1() << endl
		<< "  e2 = " << kozai->get_ecc2() << endl
		<< "  g1 = " << kozai->get_g1()/DEG << endl
		<< "  g2 = " << kozai->get_g2()/DEG << endl
		<< "  Omega1 = " << kozai->get_Omega1()/DEG << endl
		<< "  Omega2 = " << kozai->get_Omega2()/DEG << endl
		<< "  inc = " << kozai->get_inc()/DEG << endl
		<< "  inc 1 = " << kozai->get_inc1()/DEG << endl
		<< "  inc 2 = " << kozai->get_inc2()/DEG << endl
		<< "  radius 1 = " << kozai->get_r1()/RSUN << endl
		<< "  radius 2 = " << kozai->get_r2()/RSUN << endl
		<< "  theta 1 = " << kozai->get_theta1()/DEG << endl
		<< "  theta 2 = " << kozai->get_theta2()/DEG << endl
		<< "  delta Phi = " << kozai->get_deltaPhi()/DEG << endl;
	
	//Print the header with units and column names, and the initial state of the
	//system
	cerr << "#1:t[yr]  #2:a1[AU]  #3:e1  #4:e2  #5:inc[deg] #6:i1  #7:i2  #8:g1  #9:g2  #10:theta1  #11:theta2  #12:dPhi #13:Theta1 #14:Theta2 #15:SS" << endl;
}


int rhs(double t, const double y[], double f[], void *kozai_ptr){

	kozai_struct *kozai = (kozai_struct *) kozai_ptr;

	//First, convert the y array into vectors
	//
	//Note: you can't use the vectors in the kozai structure; for adaptive
	//timstepping, GSL sometimes uses a temporary y array, which won't be
	//updated in the kozai class
	vec j1 = vec(y[0],y[1],y[2]);
	vec e1 = vec(y[3],y[4],y[5]);
	vec j2 = vec(y[6],y[7],y[8]);
	vec e2 = vec(y[9],y[10],y[11]);
	vec s1 = vec(y[12],y[13],y[14]);
	vec s2 = vec(y[15],y[16],y[17]);

	//extract semi-major axes, masses and angular momenta
	double a1 = y[18];
	double a2 = kozai->get_a2();
	double m1 = kozai->get_m1();
	double m2 = kozai->get_m2();
	double m3 = kozai->get_m3();
	double L1 = kozai->get_L1_no_a()*sqrt(a1);
	double L2 = kozai->get_L2();

	//pre-compute a bunch of the vectorial quantities that are used repeatedly
	double j1n = abs(j1);
	double j2n = abs(j2);
	vec n1 = j1/j1n;
	vec n2 = j2/j2n;

	double e1n = abs(e1);
	double e2n = abs(e2);
	vec u1 = e1/e1n;
	vec u2 = e2/e2n;

	double j1n2 = j1*n2;
	double j1u2 = j1*u2;
	double e1n2 = e1*n2;
	double e1u2 = e1*u2;

	vec j1xn2 = j1^n2;
	vec j1xu2 = j1^u2;
	vec e1xn2 = e1^n2;
	vec e1xu2 = e1^u2;

	double L1L2 = L1 / L2;

	//Compute the (local) quadrupole timescale
	double tsec = (sqrt((m1+m2)/(G*pow(a1,3.)))*pow(a2,3.)*pow(1.-sqr(e2n),1.5) /m3);

	vec dj1dt, de1dt=0., dj2dt, de2dt=0., ds1dt, ds2dt;
	double dadt  = 0;

	//Add the quadrupole-order secular evolution equations
	if(kozai->get_quadrupole()){
		double quad_coef = 0.75 / tsec;
		dj1dt += quad_coef*(j1n2*j1xn2 - 5.*e1n2*e1xn2);
		de1dt += quad_coef*((j1n2*e1xn2) + (2.*(j1^e1)) - (5.*(e1n2*j1xn2)));
		dj2dt += quad_coef*L1L2*((-j1n2*j1xn2) + (5*e1n2*e1xn2));
		de2dt += quad_coef*(L1L2/j2n)*((j1n2*(e2^j1)) - (5.*e1n2*(e2^e1)) 
			   - ((0.5 - 3.*sqr(e1n) + 12.5*sqr(e1n2) - 2.5*sqr(j1n2))*(n2^e2)));
	}

	//Add the octupole-order secular evolution equations
	if(kozai->get_octupole() == true){
		double octo_coef = 1.171875*((m1-m2)/(m1+m2)*(a1/a2)*(e2n/sqr(j2n)))/tsec;
		dj1dt += octo_coef*((((2.*(e1u2*j1n2+e1n2*j1u2)*j1)
				+ 2.*(j1u2*j1n2-7.*e1u2*e1n2)*e1)^n2) + ((((2.*e1n2*j1n2)*j1)
				+ (1.6*sqr(e1n) - 0.2 - 7.*sqr(e1n2) + sqr(j1n2))*e1)^u2));
		de1dt += octo_coef*(((((2.*e1n2*j1n2)*e1) + (1.6*sqr(e1n) - 0.2 - 7.*sqr(e1n2)   
				+ sqr(j1n2))*j1)^u2) + (((2.*(e1u2*j1n2+e1n2*j1u2)*e1) 
				+ (2.*(j1n2*j1u2-7*e1n2*e1u2)*j1))^n2) + (3.2*e1u2*(j1^e1)));
		dj2dt += octo_coef*L1L2*((2*(((e1n2*j1u2)*n2) + ((e1u2*j1n2)*n2)
				+ (e1n2*j1n2)*u2)^j1) + ((((2*j1u2*j1n2)*n2) - ((14.*e1u2*e1n2)*n2) 
				+ ((1.6*sqr(e1n) - 0.2 - 7.*sqr(e1n2) + sqr(j1n2))*u2))^e1));
		de2dt += octo_coef*(L1L2/j2n)*(2*(((e1n2*((j1*e2)*u2))
				+ (j1n2*((e1*e2)*u2)) + (((sqr(j2n)/e2n)*e1n2*j1n2)*n2))^j1) + ((((2.*e2n*j1u2*j1n2)*u2)
				- ((14.*e1n2*(e1*e2))*u2) + (sqr(j2n)/e2n)*((1.6*sqr(e1n) - 0.2 - 7.*sqr(e1n2)
				+ sqr(j1n2))*n2))^e1) - (((2.*(0.2-1.6*sqr(e1n))*e1u2*e2) + ((14.*e1n2*j1u2*j1n2)*e2)
				+ (7.*e1u2*(1.6*sqr(e1n) - 0.2 - 7.*sqr(e1n2) + sqr(j1n2))*e2))^n2));
	}

	//Add the pericenter precession of the inner binary
	if (kozai->get_pericenter() == true)
		de1dt += (3./(c*c*a1*sqr(j1n)))*(pow(G*(m1+m2)/a1, 1.5) * (n1^e1));

	//Add spin-orbit coupling for the inner binary 
	if (kozai->get_spinorbit() == true){
		vec seff = (1+0.75*kozai->get_m2m1()) * s1 + (1+0.75*kozai->get_m1m2()) * s2;

		double sep_coef = Gc2/pow(a1*j1n,3.);

		dj1dt += 2*sep_coef*(seff^j1);
		de1dt += 2*sep_coef*((-3*(seff*n1)*n1 + seff)^e1);

		ds1dt += 2*sep_coef*(1+0.75*kozai->get_m2m1())*(j1^s1)*L1;
		ds2dt += 2*sep_coef*(1+0.75*kozai->get_m1m2())*(j1^s2)*L1;
	}

	//Add spin-spin interactions for the inner binary 
	if (kozai->get_spinspin() == true){
		vec s0 = (1+kozai->get_m2m1()) * s1 + (1+kozai->get_m1m2()) * s2;

		double sep_coef = Gc2/pow(a1*j1n,3.);

		double s0n1 = s0*n1;
		double s0s0 = s0*s0;

		dj1dt += -0.75*kozai->get_eta()*sep_coef*(2*s0n1*(s0^j1)) / (L1*j1n);
		de1dt += 0.75*kozai->get_eta()*sep_coef*((5*sqr(s0n1)-s0s0)*n1-2*s0n1*s0)^e1 / (L1*j1n);

		ds1dt += 0.5*(1+kozai->get_m2m1())*kozai->get_eta()*sep_coef*(s0 - 3*s0n1*n1)^s1;
		ds2dt += 0.5*(1+kozai->get_m1m2())*kozai->get_eta()*sep_coef*(s0 - 3*s0n1*n1)^s2;
	}

	//Add gravitational-wave emission for the inner binary
	if (kozai->get_radiation() == true){
		double dedt  = -(c304o15*e1n*(G3 * m1*m2*(m1+m2)) /
				 (c5 * pow(a1,4.) * pow(j1n,5.))) * (1. + c121o304*sqr(e1n));
		dadt += -c64o5 * ((G3 * m1*m2*(m1+m2)) / (c5
				* pow(a1,3.) * pow(j1n,7.))) * (1. + c73o24*sqr(e1n) + c37o96*pow(e1n,4.)); 

		de1dt += dedt * u1;
		dj1dt += (-e1n / j1n) * dedt * n1;
	}

	//Finally, copy the derivatives into the output array
	for(int i=0; i<3; i++){
		f[i] = dj1dt[i];
		f[i+3] = de1dt[i];
		f[i+6] = dj2dt[i];
		f[i+9] = de2dt[i];
		f[i+12] = ds1dt[i];
		f[i+15] = ds2dt[i];
	}
	f[18] = dadt;

	return GSL_SUCCESS;
}

// Integrate the peters equation to find the inspiral time for a binary
// Not actually used at the moment, but probably useful to have around...
double peters_t(kozai_struct *kozai){
	double beta = (64./5.) * G3c5 * kozai->get_m1()*kozai->get_m2()*(kozai->get_m1()+kozai->get_m2());
	double e = kozai->get_ecc1();
	double c0 = kozai->get_a1() * (1-sqr(e)) * pow(e,-c12o19) * pow(1+(c121o304)*sqr(e),-c870o2299);

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	gsl_function F;
	F.function = &peters_integral;
	double integral, error;

	gsl_integration_qags(&F,0,e,1e-7,1e-8,1000, w,&integral, &error);
	gsl_integration_workspace_free(w);

	return integral * c12o19 * sqr(c0)*sqr(c0) / beta;
}

double peters_integral(double e, void *params){
	return pow(e,c29o19) * pow(1+(c121o304)*sqr(e),c1181o2299) / pow(1-sqr(e),1.5);
}

void set_parameters(int argc, char **argv, kozai_struct *kozai, double &t_end, double &t_step, bool &IGNORE_GSL_ERRORS){
 	int c;
	int opterr = 0;
	int option_index = 0;
	double radius;

	const struct option longopts[] =
	{
		{"m1",    1,  0, 'm'},
		{"m2",    1,  0, 'M'},
		{"m3",    1,  0, 'n'},
		{"a1",    1,  0, 'a'},
		{"a2",    1,  0, 'A'},
		{"e1",    1,  0, 'e'},
		{"e2",    1,  0, 'E'},
		{"g1",    1,  0, 'g'},
		{"g2",    1,  0, 'G'},
		{"omega1",1,  0, 'l'},
		{"omega2",1,  0, 'L'},
		{"inc"   ,1,  0, 'i'},
		{"rad1"   ,1,  0, 'c'},
		{"rad2"   ,1,  0, 'C'},
		{"chi1"   ,1,  0, 'c'},
		{"chi2"   ,1,  0, 'C'},
		{"theta1"   ,1,  0, 't'},
		{"theta2"   ,1,  0, 'T'},
		{"phi1"   ,1,  0, 'u'},
		{"phi2"   ,1,  0, 'U'},
		{"quad"  ,0,  0, 'q'},
		{"oct"   ,0,  0, 'o'},
		{"peri"  ,0,  0, 'p'},
		{"spinorbit"  ,0,  0, 's'},
		{"spinspin"   ,0,  0, 'S'},
		{"rad"   ,0,  0, 'r'},
		{"time"   ,1,  0, 'd'},
		{"dt"   ,1,  0, 'D'},
		{"ignore_gsl"   ,0,  0, 'I'},
		{"help"   ,0,  0, 'h'},
		{0,0,0,0},
	};

	while ((c = getopt_long (argc, argv, 
					"m:M:n:a:A:e:E:g:G:l:L:i:c:C:t:T:u:U:qopsSrId:D:b:B:h",
					longopts,&option_index)) != -1)
	switch (c)
    {
		case 'm':
		  kozai->set_m1(atof(optarg)*MSUN); break;
		case 'M':
		  kozai->set_m2(atof(optarg)*MSUN); break;
		case 'n':
		  kozai->set_m3(atof(optarg)*MSUN); break;
		case 'a':
		  kozai->set_a1(atof(optarg)*AU); break;
		case 'A':
		  kozai->set_a2(atof(optarg)*AU); break;
		case 'e':
		  kozai->set_ecc1(atof(optarg)); break;
		case 'E':
		  kozai->set_ecc2(atof(optarg)); break;
		case 'g':
		  kozai->set_g1(atof(optarg)*DEG); break;
		case 'G':
		  kozai->set_g2(atof(optarg)*DEG); break;
		case 'l':
		  kozai->set_Omega1(atof(optarg)*DEG); break;
		case 'L':
		  kozai->set_Omega2(atof(optarg)*DEG); break;
		case 'i':
		  kozai->set_inc(atof(optarg)*DEG); break;
		case 'b':
		  kozai->set_r1(atof(optarg)*RSUN); break;
		case 'B':
		  kozai->set_r2(atof(optarg)*RSUN); break;
		case 'c':
		  kozai->set_chi1(atof(optarg)); break;
		case 'C':
		  kozai->set_chi2(atof(optarg)); break;
		case 't':
		  kozai->set_theta1(atof(optarg)*DEG); break;
		case 'T':
		  kozai->set_theta2(atof(optarg)*DEG); break;
		case 'u':
		  kozai->set_phi1(atof(optarg)*DEG); break;
		case 'U':
		  kozai->set_phi2(atof(optarg)*DEG); break;
		case 'q':
		  kozai->set_quadrupole(true); break;
		case 'o':
		  kozai->set_octupole(true); break;
		case 'p':
		  kozai->set_pericenter(true); break;
		case 's':
		  kozai->set_spinorbit(true); break;
		case 'S':
		  kozai->set_spinspin(true); break;
		case 'r':
		  kozai->set_radiation(true); break;
		case 'I':
		  IGNORE_GSL_ERRORS = true; break;
		case 'd':
		  t_end = atof(optarg)*YEAR; break;
		case 'D':
		  t_step = atof(optarg);
		  if (t_step <= 0) break;
		  t_step*=YEAR;
		  break;
		case 'h':
			cout << help;
			exit(1);
		case '?':
			cout << help;
			exit(1);
		default:
			exit(1);
      }

}
