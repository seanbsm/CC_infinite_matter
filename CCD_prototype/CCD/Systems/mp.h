#ifndef MP_H
#define MP_H

#include "master.h"
#include "system.h"
#include <math.h>
#include <cmath>
#include <eigen3/Eigen/Dense>

class MP : public System
{
private:
    bool vecDelta(Eigen::VectorXi v1, Eigen::VectorXi v2);

    const double V_0R_fac = 200*pow(pi/kappa_R, 1.5);    //MeV
    const double V_0T_fac = 178*pow(pi/kappa_T, 1.5);    //MeV
    const double V_0S_fac = 91.85*pow(pi/kappa_S, 1.5);  //MeV
    const double kappa_R  = 1.487;  //fm^-2
    const double kappa_T  = 0.639;  //fm^-2
    const double kappa_S  = 0.465;  //fm^-2

public:
    double          pi     = M_PI;
    double          m_m    = 0;
    double          m_L3   = 0;
    double          m_L2   = 0;
    double          m_L1   = 0;

    int             m_Nh   = 0;
    int             m_Nb   = 0;
    int             m_Ns   = 0;

    Eigen::VectorXi below_fermi;
    Eigen::VectorXi above_fermi;
    Eigen::MatrixXi m_states;

    MP(class Master* master, double m, double L3, double L2, double L1);
    class   Master* m_master = nullptr;
    void    makeStateSpace  ();
    double  assym           (int p, int q, int r, int s);
    double  h0              (int p);
    double  f               (int p);
    int     kUnique2        (int p, int q);
};

#endif // MP_H
