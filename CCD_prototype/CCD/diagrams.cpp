#include "diagrams.h"

#include <eigen3/Eigen/Dense>

Diagrams::Diagrams()
{
}

void Diagrams::setIntClass(class MakeIntMat* intClass){
    m_intClass = intClass;
}

void Diagrams::setAmpClass(class MakeAmpMat* ampClass){
    m_ampClass = ampClass;
}

Eigen::MatrixXf Diagrams::La(int i1, int i2){
    Eigen::MatrixXf product = m_ampClass->Amplitudes[i1]*m_intClass->Vpppp[i2];
    return 0.5*product;
}

Eigen::MatrixXf Diagrams::Lb(int i1, int i2){
    Eigen::MatrixXf product = m_ampClass->Amplitudes[i1]*m_intClass->Vhhhh[i2];
    return 0.5*product;
}

Eigen::MatrixXf Diagrams::Lc(int i1, int i2){
    Eigen::MatrixXf product ;
    return -product;
}

Eigen::MatrixXf Diagrams::Qa(int i1, int i2){
    Eigen::MatrixXf product;
    return 0.25*product;
}

Eigen::MatrixXf Diagrams::Qb(int i1, int i2){
    Eigen::MatrixXf product;
    return 0.5*product;
}

Eigen::MatrixXf Diagrams::Qc(int i1, int i2){
    Eigen::MatrixXf product;
    return -0.5*product;
}

Eigen::MatrixXf Diagrams::Qd(int i1, int i2){
    Eigen::MatrixXf product;
    return -0.5*product;
}
