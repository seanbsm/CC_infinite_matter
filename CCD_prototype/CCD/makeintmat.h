#ifndef MAKEINTMAT_H
#define MAKEINTMAT_H

#include <eigen3/Eigen/Dense>
#include <Systems/system.h>
#include <Systems/heg.h>
#include <iostream>

class MakeIntMat
{
private:
    int m_Nh = 0;
    int m_Ns = 0;
    //index-keepers
    Eigen::VectorXi identify_pppp;
    Eigen::VectorXi identify_hhpp;
    Eigen::VectorXi identify_fock;
    //memory storage
    Eigen::MatrixXi intMat_pppp;
    Eigen::MatrixXi intMat_hhpp;
    Eigen::MatrixXi fockMat;


    //bool contractor(int i, int j){ return i==j; } //contracts repeated elements to a single edit
public:
    //blockArrays hold quantum numbers, sortVec holds distinct kUnique
    Eigen::MatrixXi blockArrays_hh;
    Eigen::MatrixXi blockArrays_hp;
    Eigen::MatrixXi blockArrays_ph;
    Eigen::MatrixXi blockArrays_pp;
    std::vector<int> sortVec_hh;
    std::vector<int> sortVec_hp;
    std::vector<int> sortVec_ph;
    std::vector<int> sortVec_pp;
    //indexHolder hold upper and lower bound of indices for a certain kUnique, same indexing as The corresponding matrices
    Eigen::MatrixXi boundsHolder_hhpp_hh;
    Eigen::MatrixXi boundsHolder_hhpp_pp;
    //Eigen::MatrixXi indexHolder_hp;
    //Eigen::MatrixXi indexHolder_ph;
    //Eigen::MatrixXi boundsHolder_pppp_pp;

    MakeIntMat();
    System* m_system = nullptr;
    void                            mapper(int i1, int i2);
    void                            makeBlockMat(System* system, int Nh, int Ns);
    Eigen::MatrixXd                 makeSquareBlock(Eigen::MatrixXi& array, int range_lower, int range_upper);
    Eigen::MatrixXd                 makeRektBlock(Eigen::MatrixXi& array1, Eigen::MatrixXi& array2, int range_lower1, int range_upper1, int range_lower2, int range_upper2);


    //interaction matrices for CCD
    std::vector<Eigen::MatrixXd>   Vhhhh;
    std::vector<Eigen::MatrixXd>   Vhphp;
    std::vector<Eigen::MatrixXd>   Vhhpp;
    std::vector<Eigen::MatrixXd>   Vpppp;

    //vectors with kUnique for each matrix, indices match ( that is, Vhhhh_i[h] <-> Vhhhh[h] )
    std::vector<int>               Vhhhh_i;
    std::vector<int>               Vhphp_i;
    std::vector<int>               Vhhpp_i;
    std::vector<int>               Vpppp_i;
};

#endif // MAKEINTMAT_H
