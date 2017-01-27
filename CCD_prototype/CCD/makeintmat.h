#ifndef MAKEINTMAT_H
#define MAKEINTMAT_H

#include <eigen3/Eigen/Dense>
#include <Systems/system.h>
#include <Systems/heg.h>
#include <iostream>
#include <map>

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
    //blockArrays hold quantum numbers
    Eigen::MatrixXi blockArrays_h;
    Eigen::MatrixXi blockArrays_p;
    Eigen::MatrixXi blockArrays_hh;
    Eigen::MatrixXi blockArrays_hp;
    Eigen::MatrixXi blockArrays_ph;
    Eigen::MatrixXi blockArrays_pp;
    Eigen::MatrixXi blockArrays_hhp;
    Eigen::MatrixXi blockArrays_pph;

    //sortVec holds all distinct kUnique for each index series
    std::vector<int> sortVec_h;
    std::vector<int> sortVec_p;
    std::vector<int> sortVec_hh;
    std::vector<int> sortVec_hp;
    std::vector<int> sortVec_ph;
    std::vector<int> sortVec_pp;
    std::vector<int> sortVec_hhp;
    std::vector<int> sortVec_pph;

    //indexHolder holds upper and lower bound of indices for a certain kUnique, same indexing as the corresponding matrices
    Eigen::MatrixXi indexHolder_h;
    Eigen::MatrixXi indexHolder_p;
    Eigen::MatrixXi indexHolder_hh;
    Eigen::MatrixXi indexHolder_hp;
    Eigen::MatrixXi indexHolder_ph;
    Eigen::MatrixXi indexHolder_pp;
    Eigen::MatrixXi indexHolder_hhp;
    Eigen::MatrixXi indexHolder_pph;

    // these additional boundsHolders are possibly not necessary after the implementation of make3x1- and make2x2Block
    //these are needed for Qa
    Eigen::MatrixXi boundsHolder_hhpp_hh;
    Eigen::MatrixXi boundsHolder_hhpp_pp;

    //these are needed for Qc and Qd
    Eigen::MatrixXi indexHolder_h_pph_hpp;
    Eigen::MatrixXi indexHolder_h_pph_h;
    Eigen::MatrixXi indexHolder_hhp_p_hhp;
    Eigen::MatrixXi indexHolder_hhp_p_p;

    //Eigen::MatrixXi indexHolder_hp;
    //Eigen::MatrixXi indexHolder_ph;
    //Eigen::MatrixXi boundsHolder_pppp_pp;

    MakeIntMat();
    System* m_system = nullptr;
    void                            mapper_1(int i1);
    void                            mapper_2(int i1, int i2);
    void                            mapper_3(int i1, int i2, int i3);
    void                            makeBlockMat(System* system, int Nh, int Ns);
    int                             numOfKu;        //number of blocks in V_hh_pp
    Eigen::MatrixXd                 makeSquareBlock(Eigen::MatrixXi& array, int range_lower, int range_upper);
    Eigen::MatrixXd                 makeRektBlock(Eigen::MatrixXi& array1, Eigen::MatrixXi& array2, int range_lower1, int range_upper1, int range_lower2, int range_upper2);

    Eigen::MatrixXd                 make3x1Block(int ku, int i1, int i2, int i3, int i4);
    Eigen::MatrixXd                 make2x2Block(int ku, int i1, int i2, int i3, int i4);
    Eigen::MatrixXd                 make2x2Block_alt(int channel);

    void                            makeMatMap(Eigen::MatrixXi& array1, Eigen::MatrixXi& array2, int range_lower1, int range_upper1, int range_lower2, int range_upper2);
    void                            makeMatVec(Eigen::MatrixXi& array1, Eigen::MatrixXi& array2, int range_lower1, int range_upper1, int range_lower2, int range_upper2);
    int                             Identity(int h1, int h2, int p1, int p2);
    std::map<int, double>           Vhhpp_elements;

    //Eigen::VectorXd                 Vhhpp_vector;
    std::vector<double>             Vhhpp_vector;
    int                             Vhhpp_counter = 0;

    //interaction matrices for CCD
    std::vector<Eigen::MatrixXd>   Vhphp;
    std::vector<Eigen::MatrixXd>   Vhhpp;

    //these are special
    std::vector<Eigen::MatrixXd>   Vhhhh;
    std::vector<Eigen::MatrixXd>   Vpppp;

    std::vector<Eigen::MatrixXi>   V_hp_hp;
    std::vector<Eigen::MatrixXi>   V_hh_pp;

    //vectors with kUnique for each matrix, indices match ( that is, Vhhhh_i[h] <-> Vhhhh[h] )
    std::vector<int>               Vhhhh_i;
    std::vector<int>               Vhphp_i;
    std::vector<int>               Vhhpp_i;
    std::vector<int>               Vpppp_i;
};

#endif // MAKEINTMAT_H
