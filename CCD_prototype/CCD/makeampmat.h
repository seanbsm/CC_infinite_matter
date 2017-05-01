#ifndef MAKEAMPMAT_H
#define MAKEAMPMAT_H

#include "eigen3/Eigen/Dense"
#include <makeintmat.h>
#include <Systems/system.h>
#include <Systems/heg.h>
#include <Systems/mp.h>
#include <Systems/pm.h>

class MakeAmpMat
{
public:
    MakeAmpMat();

    void makeFockMaps();
    void emptyFockMaps();
    spp::sparse_hash_map<int, double> FockMap_h;
    spp::sparse_hash_map<int, double> FockMap_p;
    std::vector<Eigen::MatrixXd> denomMat;
    std::vector<Eigen::MatrixXd> denomMat3;

    int m_Counter = 0; //a test counter, used for debugging sessions

    void setThreads(int numthreads);
    int m_numThreads;

    MakeIntMat*  m_intClass = nullptr;
    System*      m_system   = nullptr;
    std::vector<Eigen::MatrixXd>  Amplitudes;

    spp::sparse_hash_map<int, double>           T2_elements;
    spp::sparse_hash_map<int, double>           T2_temp;
    spp::sparse_hash_map<int, double>           T2_elements_new;

    spp::sparse_hash_map<int, double>           T3_elements;
    spp::sparse_hash_map<int, double>           T3_temp;
    spp::sparse_hash_map<int, double>           T3_elements_new;

    //The T5b and T5c use very demanding remappings, so it's far more efficient, both in CPU and memory,
    //to have matrices storing the indices, rather than finding them on the go.
    //These matrices are made in mapper_5 in the intClass, because this is a late fix
    std::vector<Eigen::MatrixXi>              T3_T5b_indices;
    std::vector<Eigen::MatrixXi>              T3_T5c_indices;

    //std::unordered_map<int, int>              T2_elements_I;
    //std::vector<double>                       T2_elements_A;
    //std::map<unsigned long int, int> T3_elements_I_um;
    //std::map<unsigned long int, int>        T3_elements_I;        //holds indices to T3_elements_A, same for _new and _temp
    spp::sparse_hash_map<unsigned long int, unsigned long int>        T3_elements_I;
    //std::vector<spp::sparse_hash_map<unsigned long int, int>> T3_elements_IV;
    std::vector<double>                       T3_elements_A;        //holds T3 amplitudes
    std::vector<double>                       T3_elements_A_new;    //holds new T3 amplitudes
    std::vector<double>                       T3_elements_A_temp;   //holds temporary diagram contributions
    void                                      T3_makeMap(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, int i5, int i6);
    void                                      T3_makeDirectMat();
    Eigen::MatrixXd                           T3_buildDirectMat(int channel, std::vector<double>& T_vec);
    std::vector<Eigen::MatrixXi>              T3_directMat;     //holds indices for T3_elements_A to make t_ijk^abc

    Eigen::MatrixXd                 make3x1Block(int ku, int i1, int i2, int i3, int i4, spp::sparse_hash_map<int, double> &T_list);
    Eigen::MatrixXd                 make2x2Block(int ku, int i1, int i2, int i3, int i4, spp::sparse_hash_map<int, double>& T_list);

    Eigen::MatrixXd                 make3x3Block(int ku, int i1, int i2, int i3, int i4, int i5, int i6, spp::sparse_hash_map<int, double>& T_list);
    Eigen::MatrixXd                 make3x3Block_I(int ku, int i1, int i2, int i3, int i4, int i5, int i6, std::vector<double>& T_vec);
    Eigen::MatrixXd                 make3x3Block_I_D10c(int ku, int i1, int i2, int i3, int i4, int i5, int i6, std::vector<double>& T_vec);



    //CCD terms
    Eigen::MatrixXd                 I1_makemat_1(int channel1, int channel2);
    Eigen::MatrixXd                 I1_makemat_2(int channel1, int channel2);
    Eigen::MatrixXd                 I2_makemat_1(int channel1, int channel2);
    Eigen::MatrixXd                 I2_makemat_2(int channel1, int channel2);
    Eigen::MatrixXd                 I3_makemat_1(int channel1, int channel2);
    Eigen::MatrixXd                 I3_makemat_2(int channel1, int channel2);
    Eigen::MatrixXd                 I4_makemat_1(int channel1, int channel2);
    Eigen::MatrixXd                 I4_makemat_2(int channel1, int channel2);

    //T3 contributions to T2
    Eigen::MatrixXd                 D10b_makemat(int channel1, int channel2);
    Eigen::MatrixXd                 D10c_makemat(int channel1, int channel2);

    //linear T2 terms in T3
    Eigen::MatrixXd                 T1a_makemat(int channel1, int channel2);
    Eigen::MatrixXd                 T1b_makemat(int channel1, int channel2);

    //linear T3 terms in T3
    Eigen::MatrixXd                 T2c_makemat(int channel1, int channel2);
    Eigen::MatrixXd                 T2d_makemat(int channel1, int channel2);
    Eigen::MatrixXd                 T2e_makemat(int channel1, int channel2);


    //quadratic T2 terms in T3
    Eigen::MatrixXd                 T3b_makemat_1(int channel1, int channel2);
    Eigen::MatrixXd                 T3b_makemat_2(int channel1, int channel2);
    Eigen::MatrixXd                 T3b_makemat_3(int channel1, int channel2);  //The convention makemat_1, makemat_2, and makemat_3 comes from the diagrams
    Eigen::MatrixXd                 T3c_makemat_1(int channel1, int channel2);  //there are three matrices in these, so you need three constructors
    Eigen::MatrixXd                 T3c_makemat_2(int channel1, int channel2);  //a pain, trust me, I know
    Eigen::MatrixXd                 T3c_makemat_3(int channel1, int channel2);
    Eigen::MatrixXd                 T3d_makemat_1(int channel1, int channel2);
    Eigen::MatrixXd                 T3d_makemat_2(int channel1, int channel2);
    Eigen::MatrixXd                 T3d_makemat_3(int channel1, int channel2);
    Eigen::MatrixXd                 T3e_makemat_1(int channel1, int channel2);
    Eigen::MatrixXd                 T3e_makemat_2(int channel1, int channel2);
    Eigen::MatrixXd                 T3e_makemat_3(int channel1, int channel2);

    //T2*T3 terms in T3
    Eigen::MatrixXd                 T5a_makemat_1(int channel1, int channel2);
    Eigen::MatrixXd                 T5a_makemat_2(int channel1, int channel2);
    Eigen::MatrixXd                 T5b_makemat_1(int channel1, int channel2);
    Eigen::MatrixXd                 T5b_makemat_2(int channel1, int channel2);
    Eigen::MatrixXi                 T5b_makemat_2_I(int channel1, int channel2); //index, not double
    Eigen::MatrixXd                 T5c_makemat_1(int channel1, int channel2);
    Eigen::MatrixXd                 T5c_makemat_2(int channel1, int channel2);
    Eigen::MatrixXi                 T5c_makemat_2_I(int channel1, int channel2); //index, not double
    Eigen::MatrixXd                 T5d_makemat_1(int channel1, int channel2);
    Eigen::MatrixXd                 T5d_makemat_2(int channel1, int channel2);
    Eigen::MatrixXd                 T5e_makemat_1(int channel1, int channel2);
    Eigen::MatrixXd                 T5e_makemat_2(int channel1, int channel2);
    Eigen::MatrixXd                 T5f_makemat_1(int channel1, int channel2);
    Eigen::MatrixXd                 T5f_makemat_2(int channel1, int channel2);
    Eigen::MatrixXd                 T5g_makemat_1(int channel1, int channel2);
    Eigen::MatrixXd                 T5g_makemat_2(int channel1, int channel2);

    //CCD terms
    void                            I1_inverse(Eigen::MatrixXd& inMat, int channel1, int channel2);
    void                            I2_inverse(Eigen::MatrixXd& inMat, int channel1, int channel2);
    void                            I3_inverse(Eigen::MatrixXd inMat, int channel1, int channel2);
    void                            I4_inverse(Eigen::MatrixXd inMat, int channel1, int channel2);

    //T3 contributions to T2
    void                            D10b_inverse(Eigen::MatrixXd& inMat, int channel1, int channel2);
    void                            D10c_inverse(Eigen::MatrixXd& inMat, int channel1, int channel2);

    //linear T2 terms in T3
    void                            T1a_inverse(Eigen::MatrixXd& inMat, int channel1, int channel2);
    void                            T1b_inverse(Eigen::MatrixXd& inMat, int channel1, int channel2);

    //linear T3 terms in T3
    void                            T2c_inverse(Eigen::MatrixXd& inMat, int channel1, int channel2);
    void                            T2d_inverse(Eigen::MatrixXd& inMat, int channel1, int channel2);
    void                            T2e_inverse(Eigen::MatrixXd& inMat, int channel1, int channel2);

    //since the T3b-e diagrams are more finicky than the rest, we need a remapper for the first product
    spp::sparse_hash_map<unsigned long int, double> T3D_remap;
    void                            T3b_Inverse_temp(Eigen::MatrixXd& inMat, int channel1, int channel2);
    void                            T3c_Inverse_temp(Eigen::MatrixXd inMat, int channel1, int channel2);
    void                            T3d_Inverse_temp(Eigen::MatrixXd& inMat, int channel1, int channel2);
    void                            T3e_Inverse_temp(Eigen::MatrixXd& inMat, int channel1, int channel2);

    //quadratic T2 terms in T3
    void                            T3b_inverse(Eigen::MatrixXd& inMat, int channel1, int channel2);
    void                            T3c_inverse(Eigen::MatrixXd& inMat, int channel1, int channel2);
    void                            T3d_inverse(Eigen::MatrixXd& inMat, int channel1, int channel2);
    void                            T3e_inverse(Eigen::MatrixXd& inMat, int channel1, int channel2);

    //T2*T3 terms in T3
    void                            T5a_inverse(Eigen::MatrixXd& inMat, int channel1, int channel2);
    void                            T5b_inverse(Eigen::MatrixXd& inMat, int channel1, int channel2);
    void                            T5c_inverse(Eigen::MatrixXd& inMat, int channel1, int channel2);
    void                            T5b_inverse_I(Eigen::MatrixXd& inMat, int channel);
    void                            T5c_inverse_I(Eigen::MatrixXd& inMat, int channel);
    void                            T5d_inverse(Eigen::MatrixXd& inMat, int channel1, int channel2);
    void                            T5e_inverse(Eigen::MatrixXd& inMat, int channel1, int channel2);
    void                            T5f_inverse(Eigen::MatrixXd inMat, int channel1, int channel2);
    void                            T5g_inverse(Eigen::MatrixXd& inMat, int channel1, int channel2);



    void                            make3x1Block_inverse(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, spp::sparse_hash_map<int, double>& T_list, bool add);
    void                            make3x1Block_inverse_D10b(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, spp::sparse_hash_map<int, double>& T_list, bool add);
    void                            make2x2Block_inverse(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, spp::sparse_hash_map<int, double>& T_list, bool add);

    void                            make3x3Block_inverse(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, int i5, int i6, spp::sparse_hash_map<int, double>& T_list, bool add);
    void                            make3x3Block_inverse_I(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, int i5, int i6, std::vector<double>& T_vec, bool add);
    void                            make3x3Block_inverse_I_T1a(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, int i5, int i6, std::vector<double>& T_vec, bool add);
    void                            make3x3Block_inverse_I_T1b(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, int i5, int i6, std::vector<double>& T_vec, bool add);

    void                            addElementsT2(bool Pij, bool Pab);
    void                            addElementsT3_T1a();
    void                            addElementsT3_T1b();
    void                            addElementsT3_T2c();
    void                            addElementsT3_T2d();
    void                            addElementsT3_T2e();
    void                            addElementsT3_T3b();
    void                            addElementsT3_T3c();
    void                            addElementsT3_T3d();
    void                            addElementsT3_T3e();
    void                            addElementsT3_T5a();
    void                            addElementsT3_T5b();
    void                            addElementsT3_T5c();
    void                            addElementsT3_T5d();
    void                            addElementsT3_T5e();
    void                            addElementsT3_T5f();
    void                            addElementsT3_T5g();
    void                            addElementsT3(bool Pij, bool Pik, bool Pjk, bool Pab, bool Pac, bool Pbc);

    spp::sparse_hash_map<int, int> permuteT3(int index, spp::sparse_hash_map<int, int> indices);

    void setIntClass(class MakeIntMat* intClass);
    void setSystem(class System* system);
    void setElements_T2();
    void setElements_T3();
    Eigen::MatrixXd makeBlockMat(int index);
    void makeDenomMat();
    void makeDenomMat3();
};

#endif // MAKEAMPMAT_H
