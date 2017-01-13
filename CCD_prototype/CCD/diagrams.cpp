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

Eigen::MatrixXd Diagrams::La(int index){
    Eigen::MatrixXd product = m_ampClass->Amplitudes[index]*m_intClass->Vpppp[index];
    return 0.5*product;
}

Eigen::MatrixXd Diagrams::Lb(int index){
    Eigen::MatrixXd product = m_intClass->Vhhhh[index]*m_ampClass->Amplitudes[index];
    return 0.5*product;
}

Eigen::MatrixXd Diagrams::Lc(int index){
    Eigen::MatrixXd product;// = m_intClass->Vhphp[index]*m_ampClass->Amplitudes;
    return -product;
}

Eigen::MatrixXd Diagrams::Qa(int index){
    Eigen::MatrixXd product = ( m_ampClass->Amplitudes[index]*m_intClass->Vhhpp[index].transpose() )*m_ampClass->Amplitudes[index];
    return 0.25*product;
}

Eigen::MatrixXd Diagrams::Qb(int index){
    Eigen::MatrixXd product; //= m_intClass->make3x1Block(ku);
    return 0.5*product;
}

Eigen::MatrixXd Diagrams::Qc(int index){
    Eigen::MatrixXd product;
    return -0.5*product;
}

Eigen::MatrixXd Diagrams::Qd(int index){
    Eigen::MatrixXd product;
    return -0.5*product;
}
