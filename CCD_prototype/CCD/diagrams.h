#ifndef DIAGRAMS_H
#define DIAGRAMS_H

#include <eigen3/Eigen/Dense>

class Diagrams
{
public:
    Diagrams();
    Eigen::MatrixXf D3(Eigen::MatrixXf ampMat, Eigen::MatrixXf intMat); //I made this only with the ladder diagram in mind
};

#endif // DIAGRAMS_H
