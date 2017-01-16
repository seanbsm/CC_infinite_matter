#ifndef MAKEAMPMAT_H
#define MAKEAMPMAT_H

#include <eigen3/Eigen/Dense>
#include <makeintmat.h>
#include <Systems/system.h>
#include <Systems/heg.h>
#include <Systems/mp.h>
#include <Systems/pm.h>

class MakeAmpMat
{
public:
    MakeAmpMat();
    std::vector<Eigen::MatrixXd> denomMat;

    MakeIntMat*  m_intClass = nullptr;
    System*      m_system   = nullptr;
    std::vector<Eigen::MatrixXd>  Amplitudes;

    std::map<int, double>           T_elements;

    Eigen::MatrixXd                 make3x1Block(int ku, int i1, int i2, int i3, int i4);
    Eigen::MatrixXd                 make2x2Block(int ku, int i1, int i2, int i3, int i4);

    void setIntClass(class MakeIntMat* intClass);
    void setSystem(class System* system);
    void setElements();
    Eigen::MatrixXd makeBlockMat(int index);
    void makeDenomMat();
};

#endif // MAKEAMPMAT_H
