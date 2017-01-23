#ifndef DIAGRAMS_H
#define DIAGRAMS_H

#include <eigen3/Eigen/Dense>
#include <makeintmat.h>
#include <makeampmat.h>

class Diagrams
{
public:
    MakeIntMat*  m_intClass = nullptr;
    MakeAmpMat*  m_ampClass = nullptr;

    Diagrams();
    void setIntClass(class MakeIntMat* intClass);
    void setAmpClass(class MakeAmpMat* ampClass);
    void La(int ku); //I made this only with the ladder diagram in mind
    void Lb(int ku);
    void Lc(int ku);
    void Qa(int ku);
    void Qb(int ku);
    void Qc(int ku);
    void Qd(int ku);

    Eigen::MatrixXd Pab(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4);
    Eigen::MatrixXd Pij();
};

#endif // DIAGRAMS_H
