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

    void setIntClass(class MakeIntMat* intClass);
    void setSystem(class System* system);
    void makeBlockMat();
    void makeDenomMat();
};

#endif // MAKEAMPMAT_H
