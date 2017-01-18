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
};

#endif // DIAGRAMS_H
