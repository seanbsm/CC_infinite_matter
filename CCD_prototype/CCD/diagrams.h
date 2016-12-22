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
    Eigen::MatrixXd La(int index); //I made this only with the ladder diagram in mind
    Eigen::MatrixXd Lb(int index);
    Eigen::MatrixXd Lc(int index);
    Eigen::MatrixXd Qa(int index);
    Eigen::MatrixXd Qb(int index);
    Eigen::MatrixXd Qc(int index);
    Eigen::MatrixXd Qd(int index);
};

#endif // DIAGRAMS_H
