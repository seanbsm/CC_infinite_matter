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
    Eigen::MatrixXf La(int index); //I made this only with the ladder diagram in mind
    Eigen::MatrixXf Lb(int index);
    Eigen::MatrixXf Lc(int index);
    Eigen::MatrixXf Qa(int index);
    Eigen::MatrixXf Qb(int index);
    Eigen::MatrixXf Qc(int index);
    Eigen::MatrixXf Qd(int index);
};

#endif // DIAGRAMS_H
