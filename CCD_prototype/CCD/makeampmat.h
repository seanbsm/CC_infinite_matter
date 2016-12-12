#ifndef MAKEAMPMAT_H
#define MAKEAMPMAT_H

#include <eigen3/Eigen/Dense>

class MakeAmpMat
{
public:
    MakeAmpMat();
    std::vector<Eigen::MatrixXf>  Amplitudes;
    void makeBlockMat(std::vector<Eigen::MatrixXf> Vhhpp);
};

#endif // MAKEAMPMAT_H
