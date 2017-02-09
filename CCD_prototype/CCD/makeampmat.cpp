#include "makeampmat.h"
#include <iostream>

MakeAmpMat::MakeAmpMat()
{
}

void MakeAmpMat::setIntClass(class MakeIntMat* intClass){
    m_intClass = intClass;
}

void MakeAmpMat::setSystem(class System* system){
    m_system = system;
}

void MakeAmpMat::setElements_T2(){
    Amplitudes = m_intClass->Vhhpp;
    T2_elements = m_intClass->Vhhpp_elements;
}

void MakeAmpMat::setElements_T3(){

}

// "index" is the index for kUnique of Vhhpp[index],
// thus we use makeBlockMat to reconstruct block "index" for some kUnique from T2_elements
Eigen::MatrixXd MakeAmpMat::makeBlockMat(int index){
    int range_lower_hh = m_intClass->boundsHolder_hhpp_hh(0,index);
    int range_upper_hh = m_intClass->boundsHolder_hhpp_hh(1,index);
    int range_lower_pp = m_intClass->boundsHolder_hhpp_pp(0,index);
    int range_upper_pp = m_intClass->boundsHolder_hhpp_pp(1,index);

    int dim_hh = range_upper_hh - range_lower_hh;
    int dim_pp = range_upper_pp - range_lower_pp;
    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(dim_hh, dim_pp);
    for (int i = range_lower_hh; i<range_upper_hh; i++){
        for (int j = range_lower_pp; j<range_upper_pp; j++){
            returnMat(i-range_lower_hh, j-range_lower_pp) = m_intClass->Vhhpp_elements[m_intClass->Identity_hhpp((m_intClass->blockArrays_pp_hh)(1,i),
                                                                                                            (m_intClass->blockArrays_pp_hh)(2,i),
                                                                                                            (m_intClass->blockArrays_pp_pp)(1,j),
                                                                                                            (m_intClass->blockArrays_pp_pp)(2,j))];
        }
    }
    return returnMat;
}

//for now we store the Fock matrix, but perhaps later it would be better to calculate it "on the fly"
//I deliberately declared a lot of ints here, but merely for easier debugging and readability
void MakeAmpMat::makeDenomMat(){
    //for (int i=0; i<m_intClass->sortVec_hh.size(); i++){   //remember, Vhhpp_i holds all kUnique for Vhhpp
    for (int i=0; i<m_intClass->numOfKu; i++){
        int lowBound_hh  = m_intClass->boundsHolder_hhpp_hh(0,i);
        int highBound_hh = m_intClass->boundsHolder_hhpp_hh(1,i);
        int lowBound_pp  = m_intClass->boundsHolder_hhpp_pp(0,i);
        int highBound_pp = m_intClass->boundsHolder_hhpp_pp(1,i);

        int dim_hh = highBound_hh - lowBound_hh;
        int dim_pp = highBound_pp - lowBound_pp;
        Eigen::MatrixXd newMat;
        newMat.conservativeResize(dim_hh, dim_pp);

        for (int hh=lowBound_hh; hh<highBound_hh; hh++){
            for (int pp=lowBound_pp; pp<highBound_pp; pp++){
                int ii = m_intClass->blockArrays_pp_hh(1,hh);
                int jj = m_intClass->blockArrays_pp_hh(2,hh);
                int aa = m_intClass->blockArrays_pp_pp(1,pp);
                int bb = m_intClass->blockArrays_pp_pp(2,pp);
                newMat(hh-lowBound_hh, pp-lowBound_pp) = 1/( (double)(m_system->f(ii) + m_system->f(jj) - m_system->f(aa) - m_system->f(bb)) );
            }
        }
        denomMat.push_back( newMat );
    }
}

void MakeAmpMat::make3x1Block_inverse(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, std::map<int, double>& T_list, bool add){

    bool cond_hhp = (i1 == 0 && i2 == 0 && i3==1);
    bool cond_pph = (i1 == 1 && i2 == 1 && i3==0);

    bool cond_h   = (i4 == 0);
    bool cond_p   = (i4 == 1);


    Eigen::MatrixXi blockArrays1_pointer;
    Eigen::MatrixXi blockArrays2_pointer;

    std::vector<int> sortVec1;
    std::vector<int> sortVec2;

    Eigen::MatrixXi indexHolder1_pointer;
    Eigen::MatrixXi indexHolder2_pointer;


    // 0 0 1
    if (cond_hhp){
        blockArrays1_pointer = m_intClass->blockArrays_ppm_hhp;
        sortVec1             = m_intClass->sortVec_ppm_hhp;
        indexHolder1_pointer = m_intClass->indexHolder_ppm_hhp;
    }
    // 0 1 1
    else if (cond_pph){
        blockArrays1_pointer = m_intClass->blockArrays_ppm_pph;
        sortVec1             = m_intClass->sortVec_ppm_pph;
        indexHolder1_pointer = m_intClass->indexHolder_ppm_pph;
    }

    // 0
    if (cond_h){
        blockArrays2_pointer = m_intClass->blockArrays_p_h;
        sortVec2             = m_intClass->sortVec_p_h;
        indexHolder2_pointer = m_intClass->indexHolder_p_h;
    }
    // 1
    else if (cond_p){
        blockArrays2_pointer = m_intClass->blockArrays_p_p;
        sortVec2             = m_intClass->sortVec_p_p;
        indexHolder2_pointer = m_intClass->indexHolder_p_p;
    }

    int length1 = indexHolder1_pointer.cols();
    int length2 = indexHolder2_pointer.cols();

    int index1; int index2;

    auto it1 = std::find(sortVec1.begin(), sortVec1.end(), ku);
    if (it1 == sortVec1.end()){
        std::cout << "make3x1Block_inverse in MakeAmpMat, kUnique not found for rows" << std::endl;
        std::exit;
    }
    else{
      index1 = distance(sortVec1.begin(), it1);
    }

    auto it2 = std::find(sortVec2.begin(), sortVec2.end(), ku);
    if (it2 == sortVec2.end()){
        std::cout << "make3x1Block_inverse in MakeAmpMat, kUnique not found for columns" << std::endl;
        std::exit;
    }
    else{
      index2 = distance(sortVec2.begin(), it2);
    }

    int range_lower1 = indexHolder1_pointer(0,index1);
    int range_upper1 = indexHolder1_pointer(1,index1);
    int range_lower2 = indexHolder2_pointer(0,index2);
    int range_upper2 = indexHolder2_pointer(1,index2);

    int dim1 = range_upper1 - range_lower1;
    int dim2 = range_upper2 - range_lower2;

    int id;
    int ii; int jj;
    int aa; int bb;

    if (add == true){
        if (cond_hhp && cond_p){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    ii = (blockArrays1_pointer)(1,i);
                    jj = (blockArrays1_pointer)(2,i);
                    aa = (blockArrays1_pointer)(3,i);
                    bb = (blockArrays2_pointer)(1,j);
                    //id = m_intClass->Identity_hhpp((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays1_pointer)(3,i), (blockArrays2_pointer)(1,j));
                    //T_list[id] += inMat(i-range_lower1,j-range_lower2);
                    id = m_intClass->Identity_hhpp(ii,jj,aa,bb);
                    T2_temp[id] =  inMat(i-range_lower1,j-range_lower2);
                }
            }
        }
        else if (cond_pph && cond_h){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    ii = (blockArrays1_pointer)(3,i);
                    jj = (blockArrays2_pointer)(1,j);
                    aa = (blockArrays1_pointer)(1,i);
                    bb = (blockArrays1_pointer)(2,i);
                    //id = m_intClass->Identity_hhpp((blockArrays1_pointer)(3,i), (blockArrays2_pointer)(1,j), (blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i));
                    //T_list[id] += inMat(i-range_lower1,j-range_lower2);
                    id = m_intClass->Identity_hhpp(ii,jj,aa,bb);
                    T2_temp[id] =  inMat(i-range_lower1,j-range_lower2);
                }
            }
        }
    }
    else{
        if (cond_hhp && cond_p){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    id = m_intClass->Identity_hhpp((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays1_pointer)(3,i), (blockArrays2_pointer)(1,j));
                    T_list[id] = inMat(i-range_lower1,j-range_lower2);
                }
            }
        }
        else if (cond_pph && cond_h){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    id = m_intClass->Identity_hhpp((blockArrays1_pointer)(3,i), (blockArrays2_pointer)(1,j), (blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i));
                    T_list[id] = inMat(i-range_lower1,j-range_lower2);
                }
            }
        }
    }
}

void MakeAmpMat::make2x2Block_inverse(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, std::map<int, double>& T_list, bool add){

    //std::cout << "hey" << std::endl;
    bool cond_hh1 = (i1 == 0 && i2 == 0);
    bool cond_hp1 = (i1 == 0 && i2 == 1);
    bool cond_ph1 = (i1 == 1 && i2 == 0);
    bool cond_pp1 = (i1 == 1 && i2 == 1);

    bool cond_hh2 = (i3 == 0 && i4 == 0);
    bool cond_hp2 = (i3 == 0 && i4 == 1);
    bool cond_ph2 = (i3 == 1 && i4 == 0);
    bool cond_pp2 = (i3 == 1 && i4 == 1);


    Eigen::MatrixXi blockArrays1_pointer;
    Eigen::MatrixXi blockArrays2_pointer;

    std::vector<int> sortVec1;
    std::vector<int> sortVec2;

    Eigen::MatrixXi indexHolder1_pointer;
    Eigen::MatrixXi indexHolder2_pointer;


    // 0 0
    if (cond_hh1){
        blockArrays1_pointer = m_intClass->blockArrays_pp_hh;
        sortVec1             = m_intClass->sortVec_pp_hh;
        indexHolder1_pointer = m_intClass->indexHolder_pp_hh;
    }
    // 0 1
    else if (cond_hp1){
        blockArrays1_pointer = m_intClass->blockArrays_pm_hp;
        sortVec1             = m_intClass->sortVec_pm_hp;
        indexHolder1_pointer = m_intClass->indexHolder_pm_hp;
    }
    // 1 0
    else if (cond_ph1){
        blockArrays1_pointer = m_intClass->blockArrays_pm_ph;
        sortVec1             = m_intClass->sortVec_pm_ph;
        indexHolder1_pointer = m_intClass->indexHolder_pm_ph;
    }
    // 1 1
    else {
        blockArrays1_pointer = m_intClass->blockArrays_pp_pp;
        sortVec1             = m_intClass->sortVec_pp_pp;
        indexHolder1_pointer = m_intClass->indexHolder_pp_pp;
    }

    // 0 0
    if (cond_hh2){
        blockArrays2_pointer = m_intClass->blockArrays_pp_hh;
        sortVec2             = m_intClass->sortVec_pp_hh;
        indexHolder2_pointer = m_intClass->indexHolder_pp_hh;
    }
    // 0 1
    else if (cond_hp2){
        blockArrays2_pointer = m_intClass->blockArrays_pm_hp;
        sortVec2             = m_intClass->sortVec_pm_hp;
        indexHolder2_pointer = m_intClass->indexHolder_pm_hp;
    }
    // 1 0
    else if (cond_ph2){
        blockArrays2_pointer = m_intClass->blockArrays_pm_ph;
        sortVec2             = m_intClass->sortVec_pm_ph;
        indexHolder2_pointer = m_intClass->indexHolder_pm_ph;
    }
    // 1 1
    else {
        blockArrays2_pointer = m_intClass->blockArrays_pp_pp;
        sortVec2             = m_intClass->sortVec_pp_pp;
        indexHolder2_pointer = m_intClass->indexHolder_pp_pp;
    }

    int length1 = indexHolder1_pointer.cols();
    int length2 = indexHolder2_pointer.cols();

    int index1; int index2;

    auto it1 = std::find(sortVec1.begin(), sortVec1.end(), ku);
    if (it1 == sortVec1.end()){
        std::cout << "make2x2Block_inverse in MakeAmpMat, kUnique not found for rows" << std::endl;
        std::exit;
    }
    else{
      index1 = distance(sortVec1.begin(), it1);
    }

    auto it2 = std::find(sortVec2.begin(), sortVec2.end(), ku);
    if (it2 == sortVec2.end()){
        std::cout << "make2x2Block_inverse in MakeAmpMat, kUnique not found for columns" << std::endl;
        std::exit;
    }
    else{
      index2 = distance(sortVec2.begin(), it2);
    }

    int range_lower1 = indexHolder1_pointer(0,index1);
    int range_upper1 = indexHolder1_pointer(1,index1);
    int range_lower2 = indexHolder2_pointer(0,index2);
    int range_upper2 = indexHolder2_pointer(1,index2);


    int id;
    int ii; int jj;
    int aa; int bb;
    //std::map<int, double> T2_temp;

    if (add == true){
        if (cond_hp1 && cond_ph2){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    ii = (blockArrays1_pointer)(1,i);
                    jj = (blockArrays2_pointer)(2,j);
                    aa = (blockArrays1_pointer)(2,i);
                    bb = (blockArrays2_pointer)(1,j);
                    //id = m_intClass->Identity_hhpp((blockArrays1_pointer)(1,i), (blockArrays2_pointer)(2,j), (blockArrays1_pointer)(2,i), (blockArrays2_pointer)(1,j));
                    id = m_intClass->Identity_hhpp(ii,jj,aa,bb);
                    T2_temp[id] =  inMat(i-range_lower1,j-range_lower2);
                    //T_list[id] += inMat(i-range_lower1,j-range_lower2);//bad, T_list is MUCH bigger than T_temp
                }
            }
        }
        else if (cond_hh1 && cond_pp2){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    ii = (blockArrays1_pointer)(1,i);
                    jj = (blockArrays1_pointer)(2,i);
                    aa = (blockArrays2_pointer)(1,j);
                    bb = (blockArrays2_pointer)(2,j);
                    //id = m_intClass->Identity_hhpp((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays2_pointer)(1,j), (blockArrays2_pointer)(2,j));
                    id = m_intClass->Identity_hhpp(ii,jj,aa,bb);
                    T2_temp[id] =  inMat(i-range_lower1,j-range_lower2);
                    //T_list[id] += inMat(i-range_lower1,j-range_lower2);
                }
            }
        }
    }
    else{
        if (cond_hp1 && cond_ph2){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    ii = (blockArrays1_pointer)(1,i);
                    jj = (blockArrays2_pointer)(2,j);
                    aa = (blockArrays1_pointer)(2,i);
                    bb = (blockArrays2_pointer)(1,j);
                    //id = m_intClass->Identity_hhpp((blockArrays1_pointer)(1,i), (blockArrays2_pointer)(2,j), (blockArrays1_pointer)(2,i), (blockArrays2_pointer)(1,j));
                    id = m_intClass->Identity_hhpp(ii,jj,aa,bb);
                    T_list[id] = inMat(i-range_lower1,j-range_lower2);
                }
            }
        }
        else if (cond_hh1 && cond_pp2){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    ii = (blockArrays1_pointer)(1,i);
                    jj = (blockArrays1_pointer)(2,i);
                    aa = (blockArrays2_pointer)(1,j);
                    bb = (blockArrays2_pointer)(2,j);
                    //id = m_intClass->Identity_hhpp((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays2_pointer)(1,j), (blockArrays2_pointer)(2,j));
                    id = m_intClass->Identity_hhpp(ii,jj,aa,bb);
                    T_list[id] = inMat(i-range_lower1,j-range_lower2);
                }
            }
        }
    }
}

void MakeAmpMat::make3x3Block_inverse(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, int i5, int i6, std::map<int, double>& T_list, bool add){

    bool cond_hhh1 = (i1 == 0 && i2 == 0 && i3==0);
    bool cond_pph1 = (i1 == 1 && i2 == 1 && i3==0);
    bool cond_hhp1 = (i1 == 0 && i2 == 0 && i3==1);
    bool cond_ppp1 = (i1 == 1 && i2 == 1 && i3==1);

    bool cond_hhh2 = (i4 == 0 && i5 == 0 && i6==0);
    bool cond_pph2 = (i4 == 1 && i5 == 1 && i6==0);
    bool cond_hhp2 = (i4 == 0 && i5 == 0 && i6==1);
    bool cond_ppp2 = (i4 == 1 && i5 == 1 && i6==1);


    Eigen::MatrixXi blockArrays1_pointer;
    Eigen::MatrixXi blockArrays2_pointer;

    std::vector<int> sortVec1;
    std::vector<int> sortVec2;

    Eigen::MatrixXi indexHolder1_pointer;
    Eigen::MatrixXi indexHolder2_pointer;


    // 0 0 0
    if (cond_hhh1){
        blockArrays1_pointer = m_intClass->blockArrays_ppp_hhh;
        sortVec1             = m_intClass->sortVec_ppp_hhh;
        indexHolder1_pointer = m_intClass->indexHolder_ppp_hhh;
    }
    // 0 0 1
    else if (cond_hhp1){
        blockArrays1_pointer = m_intClass->blockArrays_ppm_hhp;
        sortVec1             = m_intClass->sortVec_ppm_hhp;
        indexHolder1_pointer = m_intClass->indexHolder_ppm_hhp;
    }
    // 1 1 0
    else if (cond_pph1){
        blockArrays1_pointer = m_intClass->blockArrays_ppm_pph;
        sortVec1             = m_intClass->sortVec_ppm_pph;
        indexHolder1_pointer = m_intClass->indexHolder_ppm_pph;
    }
    // 1 1 1
    else {
        blockArrays1_pointer = m_intClass->blockArrays_ppp_ppp;
        sortVec1             = m_intClass->sortVec_ppp_ppp;
        indexHolder1_pointer = m_intClass->indexHolder_ppp_ppp;
    }

    // 0 0 0
    if (cond_hhh2){
        blockArrays2_pointer = m_intClass->blockArrays_ppp_hhh;
        sortVec2             = m_intClass->sortVec_ppp_hhh;
        indexHolder2_pointer = m_intClass->indexHolder_ppp_hhh;
    }
    // 0 0 1
    else if (cond_hhp2){
        blockArrays2_pointer = m_intClass->blockArrays_ppm_hhp;
        sortVec2             = m_intClass->sortVec_ppm_hhp;
        indexHolder2_pointer = m_intClass->indexHolder_ppm_hhp;
    }
    // 1 1 0
    else if (cond_pph2){
        blockArrays2_pointer = m_intClass->blockArrays_ppm_pph;
        sortVec2             = m_intClass->sortVec_ppm_pph;
        indexHolder2_pointer = m_intClass->indexHolder_ppm_pph;
    }
    // 1 1 1
    else {
        blockArrays2_pointer = m_intClass->blockArrays_ppp_ppp;
        sortVec2             = m_intClass->sortVec_ppp_ppp;
        indexHolder2_pointer = m_intClass->indexHolder_ppp_ppp;
    }

    int length1 = indexHolder1_pointer.cols();
    int length2 = indexHolder2_pointer.cols();

    int index1; int index2;

    auto it1 = std::find(sortVec1.begin(), sortVec1.end(), ku);
    if (it1 == sortVec1.end()){
        std::cout << "make3x1Block_inverse in MakeAmpMat, kUnique not found for rows" << std::endl;
        std::exit;
    }
    else{
      index1 = distance(sortVec1.begin(), it1);
    }

    auto it2 = std::find(sortVec2.begin(), sortVec2.end(), ku);
    if (it2 == sortVec2.end()){
        std::cout << "make3x1Block_inverse in MakeAmpMat, kUnique not found for columns" << std::endl;
        std::exit;
    }
    else{
      index2 = distance(sortVec2.begin(), it2);
    }

    int range_lower1 = indexHolder1_pointer(0,index1);
    int range_upper1 = indexHolder1_pointer(1,index1);
    int range_lower2 = indexHolder2_pointer(0,index2);
    int range_upper2 = indexHolder2_pointer(1,index2);

    int dim1 = range_upper1 - range_lower1;
    int dim2 = range_upper2 - range_lower2;

    int id;
    int ii; int jj; int kk;
    int aa; int bb; int cc;

    if (add == true){
        if (cond_hhp1 && cond_pph2){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    ii = (blockArrays1_pointer)(1,i);
                    jj = (blockArrays1_pointer)(2,i);
                    kk = (blockArrays2_pointer)(3,j);
                    aa = (blockArrays2_pointer)(1,j);
                    bb = (blockArrays2_pointer)(2,j);
                    cc = (blockArrays1_pointer)(3,i);
                    id = m_intClass->Identity_hhhppp(ii,jj,kk,aa,bb,cc);
                    T3_temp[id] =  inMat(i-range_lower1,j-range_lower2);
                }
            }
        }
        else if (cond_hhh1 && cond_ppp2){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    ii = (blockArrays1_pointer)(1,i);
                    jj = (blockArrays1_pointer)(2,j);
                    kk = (blockArrays1_pointer)(3,j);
                    aa = (blockArrays2_pointer)(1,i);
                    bb = (blockArrays2_pointer)(2,i);
                    cc = (blockArrays2_pointer)(3,i);
                    id = m_intClass->Identity_hhhppp(ii,jj,kk,aa,bb,cc);
                    T3_temp[id] =  inMat(i-range_lower1,j-range_lower2);
                }
            }
        }
    }
    else{
        if (cond_hhp1 && cond_pph2){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    id = m_intClass->Identity_hhpp((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays1_pointer)(3,i), (blockArrays2_pointer)(1,j));
                    T_list[id] = inMat(i-range_lower1,j-range_lower2);
                }
            }
        }
        else if (cond_hhh1 && cond_ppp2){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    id = m_intClass->Identity_hhpp((blockArrays1_pointer)(3,i), (blockArrays2_pointer)(1,j), (blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i));
                    T_list[id] = inMat(i-range_lower1,j-range_lower2);
                }
            }
        }
    }
}

void MakeAmpMat::addElementsT2(bool Pij, bool Pab){
    for (int channel = 0; channel<m_intClass->numOfKu; channel++){
        int range_lower1 = m_intClass->boundsHolder_hhpp_hh(0,channel);
        int range_upper1 = m_intClass->boundsHolder_hhpp_hh(1,channel);
        int range_lower2 = m_intClass->boundsHolder_hhpp_pp(0,channel);
        int range_upper2 = m_intClass->boundsHolder_hhpp_pp(1,channel);

        int id; int id_Pij; int id_Pab; int id_Pijab;

        int ii; int jj;
        int aa; int bb;

        double val;

        for (int hh = range_lower1; hh<range_upper1; hh++){
            for (int pp = range_lower2; pp<range_upper2; pp++){
                ii = (m_intClass->blockArrays_pp_hh)(1,hh);
                jj = (m_intClass->blockArrays_pp_hh)(2,hh);
                aa = (m_intClass->blockArrays_pp_pp)(1,pp);
                bb = (m_intClass->blockArrays_pp_pp)(2,pp);

                id = m_intClass->Identity_hhpp(ii,jj,aa,bb);

                val = 0;
                val += T2_temp[id];

                if (Pab){
                    id_Pab = m_intClass->Identity_hhpp(ii,jj,bb,aa);
                    val -= T2_temp[id_Pab];
                }
                if (Pij){
                    id_Pij = m_intClass->Identity_hhpp(jj,ii,aa,bb);
                    val -= T2_temp[id_Pij];
                }
                if (Pij && Pab){
                    id_Pijab = m_intClass->Identity_hhpp(jj,ii,bb,aa);
                    val += T2_temp[id_Pijab];
                }

                T2_elements_new[id] += val;
            }
        }
    }
}

void MakeAmpMat::addElementsT3(bool Pij, bool Pik, bool Pjk, bool Pab, bool Pac, bool Pbc){
    for (int channel = 0; channel<m_intClass->numOfKu; channel++){
        int range_lower1 = m_intClass->boundsHolder_hhpp_hh(0,channel);
        int range_upper1 = m_intClass->boundsHolder_hhpp_hh(1,channel);
        int range_lower2 = m_intClass->boundsHolder_hhpp_pp(0,channel);
        int range_upper2 = m_intClass->boundsHolder_hhpp_pp(1,channel);

        int id; int id_Pij; int id_Pab; int id_Pijab;

        int ii; int jj;
        int aa; int bb;

        double val;

        for (int hh = range_lower1; hh<range_upper1; hh++){
            for (int pp = range_lower2; pp<range_upper2; pp++){
                ii = (m_intClass->blockArrays_pp_hh)(1,hh);
                jj = (m_intClass->blockArrays_pp_hh)(2,hh);
                aa = (m_intClass->blockArrays_pp_pp)(1,pp);
                bb = (m_intClass->blockArrays_pp_pp)(2,pp);

                id = m_intClass->Identity_hhpp(ii,jj,aa,bb);

                val = 0;
                val += T2_temp[id];

                if (Pab){
                    id_Pab = m_intClass->Identity_hhpp(ii,jj,bb,aa);
                    val -= T2_temp[id_Pab];
                }
                if (Pij){
                    id_Pij = m_intClass->Identity_hhpp(jj,ii,aa,bb);
                    val -= T2_temp[id_Pij];
                }
                if (Pij && Pab){
                    id_Pijab = m_intClass->Identity_hhpp(jj,ii,bb,aa);
                    val += T2_temp[id_Pijab];
                }

                T2_elements_new[id] += val;
            }
        }
    }
}

//returns a block matrix of dimensions 3x1, currently only made for Vhhpp
// i1,i2,i3,i4 specify whether there is a hole or particle (by a 0 or 1)  index at index ij, for j=1-4
Eigen::MatrixXd MakeAmpMat::make3x1Block(int ku, int i1, int i2, int i3, int i4, std::map<int, double>& T_list){

    bool cond_hhp = (i1 == 0 && i2 == 0 && i3==1);
    bool cond_pph = (i1 == 1 && i2 == 1 && i3==0);

    bool cond_h   = (i4 == 0);
    bool cond_p   = (i4 == 1);


    Eigen::MatrixXi blockArrays1_pointer;
    Eigen::MatrixXi blockArrays2_pointer;

    std::vector<int> sortVec1;
    std::vector<int> sortVec2;

    Eigen::MatrixXi indexHolder1_pointer;
    Eigen::MatrixXi indexHolder2_pointer;


    // 0 0 1
    if (cond_hhp){
        blockArrays1_pointer = m_intClass->blockArrays_ppm_hhp;
        sortVec1             = m_intClass->sortVec_ppm_hhp;
        indexHolder1_pointer = m_intClass->indexHolder_ppm_hhp;
    }
    // 0 1 1
    else if (cond_pph){
        blockArrays1_pointer = m_intClass->blockArrays_ppm_pph;
        sortVec1             = m_intClass->sortVec_ppm_pph;
        indexHolder1_pointer = m_intClass->indexHolder_ppm_pph;
    }

    // 0
    if (cond_h){
        blockArrays2_pointer = m_intClass->blockArrays_p_h;
        sortVec2             = m_intClass->sortVec_p_h;
        indexHolder2_pointer = m_intClass->indexHolder_p_h;
    }
    // 1
    else if (cond_p){
        blockArrays2_pointer = m_intClass->blockArrays_p_p;
        sortVec2             = m_intClass->sortVec_p_p;
        indexHolder2_pointer = m_intClass->indexHolder_p_p;
    }

    int length1 = indexHolder1_pointer.cols();
    int length2 = indexHolder2_pointer.cols();

    int index1; int index2;

    Eigen::MatrixXd returnMat;

    auto it1 = std::find(sortVec1.begin(), sortVec1.end(), ku);
    if (it1 == sortVec1.end()){
        returnMat.conservativeResize(1,1);
        returnMat(0,0) = 0;
        std::cout << "make3x1Block in MakeAmpMat, kUnique not found for rows" << std::endl;
      return returnMat;
    }
    else{
      index1 = distance(sortVec1.begin(), it1);
    }

    auto it2 = std::find(sortVec2.begin(), sortVec2.end(), ku);
    if (it2 == sortVec2.end()){
        returnMat.conservativeResize(1,1);
        returnMat(0,0) = 0;
        std::cout << "make3x1Block in MakeAmpMat, kUnique not found for columns" << std::endl;
      return returnMat;
    }
    else{
      index2 = distance(sortVec2.begin(), it2);
    }

    int range_lower1 = indexHolder1_pointer(0,index1);
    int range_upper1 = indexHolder1_pointer(1,index1);
    int range_lower2 = indexHolder2_pointer(0,index2);
    int range_upper2 = indexHolder2_pointer(1,index2);

    int dim1 = range_upper1 - range_lower1;
    int dim2 = range_upper2 - range_lower2;

    returnMat.conservativeResize(dim1, dim2);
    int id;
    if (cond_hhp && cond_p){
        for (int i = range_lower1; i<range_upper1; i++){
            for (int j = range_lower2; j<range_upper2; j++){
                id = m_intClass->Identity_hhpp((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays1_pointer)(3,i), (blockArrays2_pointer)(1,j));
                returnMat(i-range_lower1, j-range_lower2) = T_list[id];
            }
        }
    }
    else if (cond_pph && cond_h){
        for (int i = range_lower1; i<range_upper1; i++){
            for (int j = range_lower2; j<range_upper2; j++){
                id = m_intClass->Identity_hhpp((blockArrays1_pointer)(3,i), (blockArrays2_pointer)(1,j), (blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i));
                returnMat(i-range_lower1, j-range_lower2) = T_list[id];
            }
        }
    }
    return returnMat;
}

//returns a block matrix of dimensions 2x2, currently only made for Vhhpp
// i1,i2,i3,i4 specify whether there is a hole or particle (by a 0 or 1) index at index ij, for j=1-4
Eigen::MatrixXd MakeAmpMat::make2x2Block(int ku, int i1, int i2, int i3, int i4, std::map<int, double>& T_list){

    //std::cout << "hey" << std::endl;
    bool cond_hh1 = (i1 == 0 && i2 == 0);
    bool cond_hp1 = (i1 == 0 && i2 == 1);
    bool cond_ph1 = (i1 == 1 && i2 == 0);
    bool cond_pp1 = (i1 == 1 && i2 == 1);

    bool cond_hh2 = (i3 == 0 && i4 == 0);
    bool cond_hp2 = (i3 == 0 && i4 == 1);
    bool cond_ph2 = (i3 == 1 && i4 == 0);
    bool cond_pp2 = (i3 == 1 && i4 == 1);


    Eigen::MatrixXi blockArrays1_pointer;
    Eigen::MatrixXi blockArrays2_pointer;

    std::vector<int> sortVec1;
    std::vector<int> sortVec2;

    Eigen::MatrixXi indexHolder1_pointer;
    Eigen::MatrixXi indexHolder2_pointer;


    // 0 0
    if (cond_hh1){
        blockArrays1_pointer = m_intClass->blockArrays_pp_hh;
        sortVec1             = m_intClass->sortVec_pp_hh;
        indexHolder1_pointer = m_intClass->indexHolder_pp_hh;
    }
    // 0 1
    else if (cond_hp1){
        blockArrays1_pointer = m_intClass->blockArrays_pm_hp;
        sortVec1             = m_intClass->sortVec_pm_hp;
        indexHolder1_pointer = m_intClass->indexHolder_pm_hp;
    }
    // 1 0
    else if (cond_ph1){
        blockArrays1_pointer = m_intClass->blockArrays_pm_ph;
        sortVec1             = m_intClass->sortVec_pm_ph;
        indexHolder1_pointer = m_intClass->indexHolder_pm_ph;
    }
    // 1 1
    else {
        blockArrays1_pointer = m_intClass->blockArrays_pp_pp;
        sortVec1             = m_intClass->sortVec_pp_pp;
        indexHolder1_pointer = m_intClass->indexHolder_pp_pp;
    }

    // 0 0
    if (cond_hh2){
        blockArrays2_pointer = m_intClass->blockArrays_pp_hh;
        sortVec2             = m_intClass->sortVec_pp_hh;
        indexHolder2_pointer = m_intClass->indexHolder_pp_hh;
    }
    // 0 1
    else if (cond_hp2){
        blockArrays2_pointer = m_intClass->blockArrays_pm_hp;
        sortVec2             = m_intClass->sortVec_pm_hp;
        indexHolder2_pointer = m_intClass->indexHolder_pm_hp;
    }
    // 1 0
    else if (cond_ph2){
        blockArrays2_pointer = m_intClass->blockArrays_pm_ph;
        sortVec2             = m_intClass->sortVec_pm_ph;
        indexHolder2_pointer = m_intClass->indexHolder_pm_ph;
    }
    // 1 1
    else {
        blockArrays2_pointer = m_intClass->blockArrays_pp_pp;
        sortVec2             = m_intClass->sortVec_pp_pp;
        indexHolder2_pointer = m_intClass->indexHolder_pp_pp;
    }

    int length1 = indexHolder1_pointer.cols();
    int length2 = indexHolder2_pointer.cols();

    int index1; int index2;

    Eigen::MatrixXd returnMat;

    auto it1 = std::find(sortVec1.begin(), sortVec1.end(), ku);
    if (it1 == sortVec1.end()){
        returnMat.conservativeResize(1,1);
        returnMat(0,0) = 0;
        std::cout << "make2x2Block in MakeAmpMat, kUnique not found for rows" << std::endl;
      return returnMat;
    }
    else{
      index1 = distance(sortVec1.begin(), it1);
    }

    auto it2 = std::find(sortVec2.begin(), sortVec2.end(), ku);
    if (it2 == sortVec2.end()){
        returnMat.conservativeResize(1,1);
        returnMat(0,0) = 0;
        std::cout << "make2x2Block in MakeAmpMat, kUnique not found for columns" << std::endl;
      return returnMat;
    }
    else{
      index2 = distance(sortVec2.begin(), it2);
    }

    int range_lower1 = indexHolder1_pointer(0,index1);
    int range_upper1 = indexHolder1_pointer(1,index1);
    int range_lower2 = indexHolder2_pointer(0,index2);
    int range_upper2 = indexHolder2_pointer(1,index2);

    int dim1 = range_upper1 - range_lower1;
    int dim2 = range_upper2 - range_lower2;

    returnMat.conservativeResize(dim1, dim2);
    if (cond_hp1 && cond_ph2){
        for (int i = range_lower1; i<range_upper1; i++){
            for (int j = range_lower2; j<range_upper2; j++){
                returnMat(i-range_lower1, j-range_lower2) = T_list[m_intClass->Identity_hhpp((blockArrays1_pointer)(1,i), (blockArrays2_pointer)(2,j), (blockArrays1_pointer)(2,i), (blockArrays2_pointer)(1,j))];
                /*int ii = (blockArrays1_pointer)(1,i);
                int jj = (blockArrays2_pointer)(2,j);
                int aa = (blockArrays1_pointer)(2,i);
                int bb = (blockArrays2_pointer)(1,j);
                if (m_system->kUnique2(ii,jj,1,1) == m_system->kUnique2(aa,bb,1,1) ){
                    if ( T_list.find(m_system->kUnique2(ii,jj,1,1)) != T_list.end() ) {
                      // not found
                        std::cout << returnMat(i-range_lower1, j-range_lower2) << std::endl;
                    }
                }*/
                //std::cout << (blockArrays1_pointer)(1,i)<<" "<< (blockArrays2_pointer)(2,j)<<" "<< (blockArrays1_pointer)(2,i) <<" "<< (blockArrays2_pointer)(1,j) << std::endl;
            }
        }
    }
    else if (cond_hh1 && cond_pp2){
        for (int i = range_lower1; i<range_upper1; i++){
            for (int j = range_lower2; j<range_upper2; j++){
                returnMat(i-range_lower1, j-range_lower2) = T_list[m_intClass->Identity_hhpp((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays2_pointer)(1,j), (blockArrays2_pointer)(2,j))];
                //std::cout << (blockArrays1_pointer)(1,i)<<" "<< (blockArrays1_pointer)(2,i)<<" "<< (blockArrays2_pointer)(1,j)<<" "<< (blockArrays2_pointer)(2,j) << std::endl;
            }
        }
    }
    return returnMat;
}
