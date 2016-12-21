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

Eigen::MatrixXf Diagrams::La(int index){
    Eigen::MatrixXf product = m_ampClass->Amplitudes[index]*m_intClass->Vpppp[index];
    return 0.5*product;
}

Eigen::MatrixXf Diagrams::Lb(int index){
    Eigen::MatrixXf product = m_intClass->Vhhhh[index].transpose()*m_ampClass->Amplitudes[index];
    return 0.5*product;
}

Eigen::MatrixXf Diagrams::Lc(int index){
    Eigen::MatrixXf product ;
    return -product;
}

Eigen::MatrixXf Diagrams::Qa(int index){
    Eigen::MatrixXf product = ( m_ampClass->Amplitudes[index]*m_intClass->Vhhpp[index].transpose() )*m_ampClass->Amplitudes[index];
    return 0.25*product;
}

Eigen::MatrixXf Diagrams::Qb(int index){
    Eigen::MatrixXf product;
    return 0.5*product;
}

Eigen::MatrixXf Diagrams::Qc(int index){
    Eigen::MatrixXf product;
    return -0.5*product;
}

Eigen::MatrixXf Diagrams::Qd(int index){
    Eigen::MatrixXf product;
    return -0.5*product;
}
