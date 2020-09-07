#ifndef HEG_H
#define HEG_H

/* External libraries */
#include <math.h>
#include <complex>
#include <iostream>
#include <cmath>
#include "../../eigen3/Eigen/Dense"

/* Self-written code */
#include "system.h"
#include "../type_definitions.h"

class HEG : public System
{
	private:
		bool vecDelta(Eigen::VectorXi v1, Eigen::VectorXi v2);
	public:
		const double pi     = M_PI;
		double       m_m    = 0;
		double       m_L3   = 0;
		double       m_L2   = 0;
		double       m_L1   = 0;
		
		int          m_Nh   = 0;
		int          m_Nb   = 0;
		int          m_Ns   = 0;
		int          m_dk   = 0;
		
		Eigen::VectorXi below_fermi;
		Eigen::VectorXi above_fermi;
		Eigen::MatrixXi m_states;
		
		HEG(double m, double L3, double L2, double L1);
		void makeStateSpace();
		void set_Nh(int Nh);
		void set_Nb(int Nb);
		void retrieve_Ns(int *Ns);
		variable_type  assym           (int p, int q, int r, int s);
		variable_type  assym_single    (int p, int q);
		variable_type  h0              (int p);
		variable_type  f               (int p);
		int     kUnique1        (int p, int s1);
		int     kUnique2        (int p, int q, int s1, int s2);
		int     kUnique3        (int p, int q, int r, int s1, int s2, int s3);
		int     kUnique4        (int p, int q, int r, int s, int s1, int s2, int s3, int s4);
		int     kUnique5        (int p, int q, int r, int s, int t, int s1, int s2, int s3, int s4, int s5);
};

#endif // HEG_H