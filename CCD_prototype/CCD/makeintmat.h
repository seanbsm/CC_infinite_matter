#ifndef MAKEINTMAT_H
#define MAKEINTMAT_H

#include "eigen3/Eigen/Dense"
#include <Systems/system.h>
#include <Systems/heg.h>
#include <Systems/mp.h>
#include <iostream>
#include <map>
#include <unordered_map>
#include <sparsepp/spp.h>

class MakeIntMat
{
private:

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
    int m_Nh        = 0;
    int m_Ns        = 0;
    int m_Np        = 0;
    int m_Nh2       = 0;
    int m_Nh3       = 0;
    int m_Nh4       = 0;
    int m_Nh5       = 0;
    int m_NhNs2     = 0;
    unsigned long int m_Nh2Ns     = 0;
    unsigned long int m_Nh3Ns     = 0;
    unsigned long int m_Nh3Ns2    = 0;

    /* Due to the momentum-conservation relation (kp+kq = kr+ks for <pq||rs>),
     * and due to alignment of matrices, we need several blockArrays, and corresponding
     * sortVecs and indexHolders. We'd like the blockArray generation to be pretty general,
     * otherwise the generation of blockArrays would require a nightmare of if-tests and whatnot.
     * So, since the CC diagrams require many different alignments, giving various rewritings of
     * the momentum relation, I've used the convention _p/m_p/h, where in the first, p/m stands for
     * plus or minus between p/h in the next part. I.e. _ppm_hpp means total momentum kh+kp-kp.
     * Probably there's a better way, but I can't muck about with a little problem for long.
     */

    //blockArrays hold quantum numbers
    Eigen::MatrixXi blockArrays_p_h;
    Eigen::MatrixXi blockArrays_p_p;
    Eigen::MatrixXi blockArrays_pp_hh;
    Eigen::MatrixXi blockArrays_pp_hp;
    Eigen::MatrixXi blockArrays_pp_ph;
    Eigen::MatrixXi blockArrays_pp_pp;
    Eigen::MatrixXi blockArrays_hp_s;   //for Vhphp
    Eigen::MatrixXi blockArrays_pm_hh;
    Eigen::MatrixXi blockArrays_pm_hp;
    Eigen::MatrixXi blockArrays_pm_ph;
    Eigen::MatrixXi blockArrays_pm_pp;
    Eigen::MatrixXi blockArrays_ppp_hhh;
    Eigen::MatrixXi blockArrays_ppm_hhh;
    Eigen::MatrixXi blockArrays_ppm_hhp;
    Eigen::MatrixXi blockArrays_ppm_pph;
    Eigen::MatrixXi blockArrays_ppp_ppp;
    Eigen::MatrixXi blockArrays_pppm_hhhp;
    Eigen::MatrixXi blockArrays_ppmm_hhpp;
    Eigen::MatrixXi blockArrays_ppmm_pphh;
    Eigen::MatrixXi blockArrays_pppm_ppph;
    Eigen::MatrixXi blockArrays_pppmm_hhhpp;
    Eigen::MatrixXi blockArrays_pppmm_ppphh;

    //T3 permutations, needed to restack t_ijk^abc matrices
    //same shape as corresponding blockArrays
    Eigen::MatrixXi blockArrays_ppp_hhh_Pij;    // ijk <-> jik
    Eigen::MatrixXi blockArrays_ppp_hhh_Pik;    // ijk <-> kji
    Eigen::MatrixXi blockArrays_ppp_hhh_Pjk;    // ijk <-> ikj
    Eigen::MatrixXi blockArrays_ppp_hhh_Pijik;  // ijk <-> kij
    Eigen::MatrixXi blockArrays_ppp_hhh_Pijjk;  // ijk <-> jki

    Eigen::MatrixXi blockArrays_ppp_ppp_Pab;    // abc <-> bac
    Eigen::MatrixXi blockArrays_ppp_ppp_Pac;    // abc <-> cba
    Eigen::MatrixXi blockArrays_ppp_ppp_Pbc;    // abc <-> acb
    Eigen::MatrixXi blockArrays_ppp_ppp_Pabac;  // abc <-> cab
    Eigen::MatrixXi blockArrays_ppp_ppp_Pabbc;  // abc <-> bca

    void makePermutations(); //fills all permutaion matrices

    //sortVec holds all distinct kUnique for each index series
    std::vector<int> sortVec_p_h;
    std::vector<int> sortVec_p_p;
    std::vector<int> sortVec_pp_hh;
    std::vector<int> sortVec_pp_hp;
    std::vector<int> sortVec_pp_ph;
    std::vector<int> sortVec_pp_pp;
    std::vector<int> sortVec_hp_s;      //for Vhphp
    std::vector<int> sortVec_pm_hh;
    std::vector<int> sortVec_pm_hp;
    std::vector<int> sortVec_pm_ph;
    std::vector<int> sortVec_pm_pp;
    std::vector<int> sortVec_ppp_hhh;
    std::vector<int> sortVec_ppm_hhh;
    std::vector<int> sortVec_ppm_hhp;
    std::vector<int> sortVec_ppm_pph;
    std::vector<int> sortVec_ppp_ppp;
    std::vector<int> sortVec_pppm_hhhp;
    std::vector<int> sortVec_ppmm_hhpp;
    std::vector<int> sortVec_ppmm_pphh;
    std::vector<int> sortVec_pppm_ppph;
    std::vector<int> sortVec_pppmm_hhhpp;
    std::vector<int> sortVec_pppmm_ppphh;

    //indexHolder holds upper and lower bound of indices for a certain kUnique, same indexing as the corresponding matrices
    Eigen::MatrixXi indexHolder_p_h;
    Eigen::MatrixXi indexHolder_p_p;
    Eigen::MatrixXi indexHolder_pp_hh;
    Eigen::MatrixXi indexHolder_pp_hp;
    Eigen::MatrixXi indexHolder_pp_ph;
    Eigen::MatrixXi indexHolder_pp_pp;
    Eigen::MatrixXi indexHolder_hp_s;   //for Vhphp
    Eigen::MatrixXi indexHolder_pm_hh;
    Eigen::MatrixXi indexHolder_pm_hp;
    Eigen::MatrixXi indexHolder_pm_ph;
    Eigen::MatrixXi indexHolder_pm_pp;
    Eigen::MatrixXi indexHolder_ppp_hhh;
    Eigen::MatrixXi indexHolder_ppm_hhh;
    Eigen::MatrixXi indexHolder_ppm_hhp;
    Eigen::MatrixXi indexHolder_ppm_pph;
    Eigen::MatrixXi indexHolder_ppp_ppp;
    Eigen::MatrixXi indexHolder_pppm_hhhp;
    Eigen::MatrixXi indexHolder_ppmm_hhpp;
    Eigen::MatrixXi indexHolder_ppmm_pphh;
    Eigen::MatrixXi indexHolder_pppm_ppph;
    Eigen::MatrixXi indexHolder_pppmm_hhhpp;
    Eigen::MatrixXi indexHolder_pppmm_ppphh;

    // these additional boundsHolders are possibly not necessary after the implementation of make3x1- and make2x2Block
    //these are needed for Qa
    Eigen::MatrixXi boundsHolder_hhpp_hh;
    Eigen::MatrixXi boundsHolder_hhpp_pp;
    Eigen::MatrixXi boundsHolder_hhhppp_hhh;
    Eigen::MatrixXi boundsHolder_hhhppp_ppp;

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
    void                            mapper_1(std::vector<int>& sortVecIn, Eigen::MatrixXi& blockArraysIn, int i1, int s1);
    void                            mapper_2(std::vector<int> &sortVecIn, Eigen::MatrixXi &blockArraysIn, int i1, int i2, int s1, int s2);
    void                            mapper_2_alt();
    void                            mapper_3(std::vector<int> &sortVecIn, Eigen::MatrixXi &blockArraysIn, int i1, int i2, int i3, int s1, int s2, int s3);
    void                            mapper_4(std::vector<int> &sortVecIn, Eigen::MatrixXi &blockArraysIn, int i1, int i2, int i3, int i4, int s1, int s2, int s3, int s4);
    void                            mapper_5(std::vector<int> &sortVecIn, Eigen::MatrixXi &blockArraysIn, int i1, int i2, int i3, int i4, int i5, int s1, int s2, int s3, int s4, int s5);
    void                            mapper_hp();        //special func made for Lc diagram, not necessary IF I fix sign convention for ph and hp

    void                            makeBlockMat(System* system, int Nh, int Ns);
    void                            setTriples(bool argument);
    bool                            m_triplesOn;

    int                             numOfKu;        //number of blocks in V_hh_pp
    int                             numOfKu3;       //number of blocks in V_hhh_ppp
    Eigen::MatrixXd                 makeSquareBlock(Eigen::MatrixXi& array, int range_lower, int range_upper);
    Eigen::MatrixXd                 makeSquareBlock_s(Eigen::MatrixXi& array, int range_lower, int range_upper);
    Eigen::MatrixXd                 makeRektBlock(Eigen::MatrixXi& array1, Eigen::MatrixXi& array2, int range_lower1, int range_upper1, int range_lower2, int range_upper2);


    //CCD terms
    Eigen::MatrixXd                 I1_makemat(int channel1, int channel2);
    Eigen::MatrixXd                 I2_makemat(int channel1, int channel2);
    Eigen::MatrixXd                 I3_makemat(int channel1, int channel2);
    Eigen::MatrixXd                 I4_makemat(int channel1, int channel2);

    //T3 contributions to T2
    Eigen::MatrixXd                 D10b_makemat(int channel1, int channel2);
    Eigen::MatrixXd                 D10c_makemat(int channel1, int channel2);

    //linear T2 terms in T3
    Eigen::MatrixXd                 T1a_makemat(int channel1, int channel2);
    Eigen::MatrixXd                 T1b_makemat(int channel1, int channel2);

    //linear T3 terms in T3 (these already stored directly, and need no "make" function)

    //quadratic T2 terms in T3
    Eigen::MatrixXd                 T3b_makemat(int channel1, int channel2);
    Eigen::MatrixXd                 T3c_makemat(int channel1, int channel2);
    Eigen::MatrixXd                 T3d_makemat(int channel1, int channel2);
    Eigen::MatrixXd                 T3e_makemat(int channel1, int channel2);

    //T2*T3 terms in T3
    Eigen::MatrixXd                 T5a_makemat(int channel1, int channel2);
    Eigen::MatrixXd                 T5b_makemat(int channel1, int channel2);
    Eigen::MatrixXd                 T5c_makemat(int channel1, int channel2);
    Eigen::MatrixXd                 T5d_makemat(int channel1, int channel2);
    Eigen::MatrixXd                 T5e_makemat(int channel1, int channel2);
    Eigen::MatrixXd                 T5f_makemat(int channel1, int channel2); //T5f and T5g are identical, but it's nice to think of them seperate
    Eigen::MatrixXd                 T5g_makemat(int channel1, int channel2);

    Eigen::MatrixXd                 make3x1Block(int ku, int i1, int i2, int i3, int i4);
    Eigen::MatrixXd                 make2x2Block(int ku, int i1, int i2, int i3, int i4);
    Eigen::MatrixXd                 make2x2Block_alt(int channel);

    void                            makeMatMap_hhhp(Eigen::MatrixXi& array1, Eigen::MatrixXi& array2, int range_lower1, int range_upper1, int range_lower2, int range_upper2);
    void                            makeMatMap_hhpp(Eigen::MatrixXi& array1, Eigen::MatrixXi& array2, int range_lower1, int range_upper1, int range_lower2, int range_upper2);
    void                            makeMatMap_ppph(Eigen::MatrixXi& array1, Eigen::MatrixXi& array2, int range_lower1, int range_upper1, int range_lower2, int range_upper2);

    void                            makeMatVec(Eigen::MatrixXi& array1, Eigen::MatrixXi& array2, int range_lower1, int range_upper1, int range_lower2, int range_upper2);
    int                             Identity_hhhp(int h1, int h2, int h3, int p1);
    int                             Identity_hhpp(int h1, int h2, int p1, int p2);
    int                             Identity_ppph(int p1, int p2, int p3, int h1);
    unsigned long int               Identity_hhhppp(int h1, int h2, int h3, int p1, int p2, int p3);
    unsigned long int               Identity_hhhhhp(int h1, int h2, int h3, int h4, int h5, int p1); //this is a special function for T3e alone
    spp::sparse_hash_map<int, double>           Vhhhp_elements; //needed for T3
    spp::sparse_hash_map<int, double>           Vhhpp_elements;
    spp::sparse_hash_map<int, double>           Vppph_elements; //needed for T3


    //Eigen::VectorXd                 Vhhpp_vector;
    std::vector<double>             Vhhpp_vector;
    int                             Vhhpp_counter = 0;

    //interaction matrices for CCD
    std::vector<Eigen::MatrixXd>   Vhhpp;

    //these are special
    std::vector<Eigen::MatrixXd>   Vpppp; //for La
    std::vector<Eigen::MatrixXd>   Vhhhh; //for Lb
    std::vector<Eigen::MatrixXd>   Vhphp; //for Lc

    std::vector<Eigen::MatrixXi>   V_hp_hp;
    std::vector<Eigen::MatrixXi>   V_hh_pp;

    //vectors with kUnique for each matrix, indices match ( that is, Vhhhh_i[h] <-> Vhhhh[h] )
    std::vector<int>               Vhhhh_i;
    std::vector<int>               Vhphp_i;
    std::vector<int>               Vhhpp_i;
    std::vector<int>               Vpppp_i;
    std::vector<int>               Vhhhppp_i;
};

#endif // MAKEINTMAT_H
