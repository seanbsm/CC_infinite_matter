#ifndef DIAGRAMS_H
#define DIAGRAMS_H

#include <eigen3/Eigen/Dense>
#include <makeintmat.h>
#include <makeampmat.h>
#include <Systems/system.h>
#include <Systems/heg.h>
#include <Systems/mp.h>
#include <Systems/pm.h>

class Diagrams
{
public:
    MakeIntMat*  m_intClass = nullptr;
    MakeAmpMat*  m_ampClass = nullptr;
    System*      m_system   = nullptr;

    Diagrams();
    void setIntClass(class MakeIntMat* intClass);
    void setAmpClass(class MakeAmpMat* ampClass);
    void setSystem(class System* system);
    void La(); //I made this only with the ladder diagram in mind
    void Lb();
    void Lc(); //index is a temp variable i need to check something
    void Qa();
    void Qb();
    void Qc();
    void Qd();

    Eigen::MatrixXd Pab(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4);
    Eigen::MatrixXd Pij();
};

#endif // DIAGRAMS_H
