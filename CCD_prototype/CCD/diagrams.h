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
    Eigen::MatrixXf La(int i1, int i2); //I made this only with the ladder diagram in mind
    Eigen::MatrixXf Lb(int i1, int i2);
    Eigen::MatrixXf Lc(int i1, int i2);
    Eigen::MatrixXf Qa(int i1, int i2);
    Eigen::MatrixXf Qb(int i1, int i2);
    Eigen::MatrixXf Qc(int i1, int i2);
    Eigen::MatrixXf Qd(int i1, int i2);
};

#endif // DIAGRAMS_H
