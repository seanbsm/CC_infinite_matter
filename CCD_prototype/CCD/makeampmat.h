#ifndef MAKEAMPMAT_H
#define MAKEAMPMAT_H

#include <eigen3/Eigen/Dense>
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
    std::unordered_map<int, double> FockMap_h;
    std::unordered_map<int, double> FockMap_p;
    std::vector<Eigen::MatrixXd> denomMat;
    std::vector<Eigen::MatrixXd> denomMat3;

    MakeIntMat*  m_intClass = nullptr;
    System*      m_system   = nullptr;
    std::vector<Eigen::MatrixXd>  Amplitudes;

    std::unordered_map<int, double>           T2_elements;
    std::unordered_map<int, double>           T2_temp;
    std::unordered_map<int, double>           T2_elements_new;

    std::unordered_map<int, double>           T3_elements;
    std::unordered_map<int, double>           T3_temp;
    std::unordered_map<int, double>           T3_elements_new;

    //std::unordered_map<int, int>              T2_elements_I;
    //std::vector<double>                       T2_elements_A;
    std::unordered_map<int, int>              T3_elements_I;        //holds indices to T3_elements_A, same for _new and _temp
    std::vector<double>                       T3_elements_A;        //holds T3 amplitudes
    std::vector<double>                       T3_elements_A_new;    //holds new T3 amplitudes
    std::vector<double>                       T3_elements_A_temp;   //holds temporary diagram contributions
    void                                      T3_makeMap(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, int i5, int i6);
    void                                      T3_makeDirectMat();
    std::vector<Eigen::MatrixXi>              T3_directMat;     //holds indices for T3_elements_A to make t_ijk^abc

    std::vector<double>             T2_elements_alt;
    std::vector<double>             T2_elements_alt_new;

    Eigen::MatrixXd                 make3x1Block(int ku, int i1, int i2, int i3, int i4, std::unordered_map<int, double>& T_list);
    Eigen::MatrixXd                 make2x2Block(int ku, int i1, int i2, int i3, int i4, std::unordered_map<int, double>& T_list);

    Eigen::MatrixXd                 make3x3Block(int ku, int i1, int i2, int i3, int i4, int i5, int i6, std::unordered_map<int, double>& T_list);
    Eigen::MatrixXd                 make3x3Block_I(int ku, int i1, int i2, int i3, int i4, int i5, int i6, std::vector<double>& T_vec);

    void                            make3x1Block_inverse(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, std::unordered_map<int, double>& T_list, bool add);
    void                            make2x2Block_inverse(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, std::unordered_map<int, double>& T_list, bool add);

    void                            make3x3Block_inverse(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, int i5, int i6, std::unordered_map<int, double>& T_list, bool add);
    void                            make3x3Block_inverse_I(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, int i5, int i6, std::vector<double>& T_vec, bool add);

    void                            addElementsT2(bool Pij, bool Pab);
    void                            addElementsT3_T1a();
    void                            addElementsT3_T1b();
    void                            addElementsT3(bool Pij, bool Pik, bool Pjk, bool Pab, bool Pac, bool Pbc);

    std::unordered_map<int, int> permuteT3(int index, std::unordered_map<int, int> indices);

    void setIntClass(class MakeIntMat* intClass);
    void setSystem(class System* system);
    void setElements_T2();
    void setElements_T3();
    Eigen::MatrixXd makeBlockMat(int index);
    void makeDenomMat();
    void makeDenomMat3();
};

#endif // MAKEAMPMAT_H
