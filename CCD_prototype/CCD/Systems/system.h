#ifndef SYSTEM_H
#define SYSTEM_H

#include <master.h>
#include <eigen3/Eigen/Dense>

class System
{
public:
    System(class Master* master);
    int             m_Nh;        //number of particles
    int             m_Nb;        //number of closed-shells
    int             m_Ns;
    Eigen::VectorXi below_fermi;
    Eigen::VectorXi above_fermi;
    Eigen::MatrixXi m_states;                         //each column is a state where the rows are quantum numbers of that states
    class Master*   m_master        = nullptr;
    virtual void    makeStateSpace  ()                              = 0;
    virtual double  assym           (int p, int q, int r, int s)    = 0;
   // virtual double  assym_single    (int p, int q)                  = 0;
    virtual double  h0              (int p)                         = 0;
    virtual double  f               (int p)                         = 0;
    virtual int     kUnique2        (int p, int q)                  = 0;
};

#endif // SYSTEM_H
