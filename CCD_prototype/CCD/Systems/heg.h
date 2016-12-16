#ifndef HEG_H
#define HEG_H

#include "master.h"
#include "system.h"
#include <math.h>
#include <cmath>
#include <eigen3/Eigen/Dense>

class HEG : public System
{
private:
    bool vecDelta(Eigen::VectorXi v1, Eigen::VectorXi v2);
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

    HEG(class Master* master, double m, double L3, double L2, double L1);
    class   Master* m_master = nullptr;
    void    makeStateSpace  ();
    double  assym           (int p, int q, int r, int s);
  //  double  assym_single    (int p, int q);
    double  h0              (int p);
    double  f               (int p);
    int     kUnique2        (int p, int q);
};

#endif // HEG_H
