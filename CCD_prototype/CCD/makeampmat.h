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
    std::vector<Eigen::MatrixXd> denomMat;

    MakeIntMat*  m_intClass = nullptr;
    System*      m_system   = nullptr;
    std::vector<Eigen::MatrixXd>  Amplitudes;

    std::map<int, double>           T2_elements;
    std::map<int, double>           T2_temp;
    std::map<int, double>           T2_elements_new;

    std::map<int, double>           T3_elements;
    std::map<int, double>           T3_temp;
    std::map<int, double>           T3_elements_new;

    std::vector<double>             T2_elements_alt;
    std::vector<double>             T2_elements_alt_new;

    Eigen::MatrixXd                 make3x1Block(int ku, int i1, int i2, int i3, int i4, std::map<int, double>& T_list);
    Eigen::MatrixXd                 make2x2Block(int ku, int i1, int i2, int i3, int i4, std::map<int, double>& T_list);

    void                            make3x1Block_inverse(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, std::map<int, double>& T_list, bool add);
    void                            make2x2Block_inverse(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, std::map<int, double>& T_list, bool add);

    void                            make3x3Block_inverse(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, int i5, int i6, std::map<int, double>& T_list, bool add);

    void                            addElementsT2(bool Pij, bool Pab);
    void                            addElementsT3(bool Pij, bool Pik, bool Pjk, bool Pab, bool Pac, bool Pbc);

    void setIntClass(class MakeIntMat* intClass);
    void setSystem(class System* system);
    void setElements_T2();
    void setElements_T3();
    Eigen::MatrixXd makeBlockMat(int index);
    void makeDenomMat();
};

#endif // MAKEAMPMAT_H
