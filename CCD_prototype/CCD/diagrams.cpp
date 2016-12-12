#include "diagrams.h"

#include <eigen3/Eigen/Dense>

Diagrams::Diagrams()
{
}

Eigen::MatrixXf Diagrams::D3(Eigen::MatrixXf ampMat, Eigen::MatrixXf intMat){
    Eigen::MatrixXf product = ampMat*intMat.transpose();
    return 0.5*product;
}
