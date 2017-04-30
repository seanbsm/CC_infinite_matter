#include "makeampmat.h"
#include <iostream>
#include <stdlib.h>

#include <omp.h>

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

    for (int hh = 0; hh<m_intClass->numOfKu; hh++){
        int ku = m_intClass->Vhhpp_i[hh];

        Eigen::MatrixXd Vhhpp = m_intClass->make2x2Block(ku,0,0,1,1);
        Eigen::MatrixXd temp = (Vhhpp).array()*denomMat[hh].array();

        make2x2Block_inverse(temp, ku, 0,0,1,1, T2_elements_new, false);
    }
    T2_elements = T2_elements_new;
}

void MakeAmpMat::setElements_T3(){

}

Eigen::MatrixXd MakeAmpMat::T3_buildDirectMat(int channel, std::vector<double>& T_vec){
    //int ku = m_intClass->Vhhhppp_i[channel];

    //Eigen::MatrixXi tempIMat = (m_ampClass->T3_directMat[channel]);
    Eigen::MatrixXd outMat;
    outMat.conservativeResize( (T3_directMat[channel]).rows(), (T3_directMat[channel]).cols() );

    //should be easily parallizable?
    for (int col=0; col<(T3_directMat[channel]).cols(); col++){
        for (int row=0; row<(T3_directMat[channel]).rows(); row++){
            outMat(row, col) = T_vec[ (T3_directMat[channel])(row, col) ];
        }
    }

    return outMat;
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

void MakeAmpMat::makeFockMaps(){
    for (int h=0; h<m_intClass->m_Nh; h++){
        FockMap_h[h] = m_system->f(h);
    }

    for (int p=m_intClass->m_Nh; p<m_intClass->m_Ns; p++){
        FockMap_p[p] = m_system->f(p);
    }
}

void MakeAmpMat::emptyFockMaps(){
    FockMap_h.clear();
    FockMap_p.clear();
}

//for now we store the Fock matrix, but perhaps later it would be better to calculate it "on the fly"
//I deliberately declared a lot of ints here, but merely for easier debugging and readability
void MakeAmpMat::makeDenomMat(){


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
                newMat(hh-lowBound_hh, pp-lowBound_pp) = 1./( (double)(FockMap_h[ii] + FockMap_h[jj] - FockMap_p[aa] - FockMap_p[bb] ) );
            }
        }
        denomMat.push_back( newMat );
    }
}

void MakeAmpMat::makeDenomMat3(){

    spp::sparse_hash_map<int, double> FockMap_h_copy = FockMap_h;
    spp::sparse_hash_map<int, double> FockMap_p_copy = FockMap_p;

    for (int i=0; i<m_intClass->numOfKu3; i++){
        int lowBound_hhh  = m_intClass->boundsHolder_hhhppp_hhh(0,i);
        int highBound_hhh = m_intClass->boundsHolder_hhhppp_hhh(1,i);
        int lowBound_ppp  = m_intClass->boundsHolder_hhhppp_ppp(0,i);
        int highBound_ppp = m_intClass->boundsHolder_hhhppp_ppp(1,i);

        int dim_hhh = highBound_hhh - lowBound_hhh;
        int dim_ppp = highBound_ppp - lowBound_ppp;
        Eigen::MatrixXd newMat;
        newMat.conservativeResize(dim_hhh, dim_ppp);

        int n = omp_get_max_threads();
        #pragma omp parallel for /*num_threads(n)*/ firstprivate(FockMap_h_copy, FockMap_p_copy)
        for (int ppp=lowBound_ppp; ppp<highBound_ppp; ppp++){
            for (int hhh=lowBound_hhh; hhh<highBound_hhh; hhh++){
                int ii = m_intClass->blockArrays_ppp_hhh(1,hhh);
                int jj = m_intClass->blockArrays_ppp_hhh(2,hhh);
                int kk = m_intClass->blockArrays_ppp_hhh(3,hhh);
                int aa = m_intClass->blockArrays_ppp_ppp(1,ppp);
                int bb = m_intClass->blockArrays_ppp_ppp(2,ppp);
                int cc = m_intClass->blockArrays_ppp_ppp(3,ppp);
                newMat(hhh-lowBound_hhh, ppp-lowBound_ppp) = 1./( (double)(FockMap_h_copy[ii] + FockMap_h_copy[jj] + FockMap_h_copy[kk] - FockMap_p_copy[aa] - FockMap_p_copy[bb] - FockMap_p_copy[cc] ) );
            }
        }
        denomMat3.push_back( newMat );
    }
}

void MakeAmpMat::make3x1Block_inverse(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, spp::sparse_hash_map<int, double> &T_list, bool add){

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

    int index1=-1; int index2=-1;

    auto it1 = std::find(sortVec1.begin(), sortVec1.end(), ku);
    if (it1 == sortVec1.end()){
        std::cout << "make3x1Block_inverse in MakeAmpMat, kUnique not found for rows" << std::endl;
    }
    else{
      index1 = distance(sortVec1.begin(), it1);
    }

    auto it2 = std::find(sortVec2.begin(), sortVec2.end(), ku);
    if (it2 == sortVec2.end()){
        std::cout << "make3x1Block_inverse in MakeAmpMat, kUnique not found for columns" << std::endl;
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

void MakeAmpMat::make3x1Block_inverse_D10b(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, spp::sparse_hash_map<int, double> &T_list, bool add){

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

    int index1=-1; int index2=-1;

    auto it1 = std::find(sortVec1.begin(), sortVec1.end(), ku);
    if (it1 == sortVec1.end()){
        std::cout << "make3x1Block_inverse in MakeAmpMat, kUnique not found for rows" << std::endl;
    }
    else{
      index1 = distance(sortVec1.begin(), it1);
    }

    auto it2 = std::find(sortVec2.begin(), sortVec2.end(), ku);
    if (it2 == sortVec2.end()){
        std::cout << "make3x1Block_inverse in MakeAmpMat, kUnique not found for columns" << std::endl;
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

    if (add == true){
        if (cond_hhp && cond_p){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    ii = (blockArrays1_pointer)(1,i);
                    jj = (blockArrays1_pointer)(2,i);
                    bb = (blockArrays1_pointer)(3,i);
                    aa = (blockArrays2_pointer)(1,j);
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

void MakeAmpMat::make2x2Block_inverse(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, spp::sparse_hash_map<int, double> &T_list, bool add){

    //std::cout << "hey" << std::endl;
    bool cond_hh1 = (i1 == 0 && i2 == 0);
    bool cond_hp1 = (i1 == 0 && i2 == 1);
    bool cond_ph1 = (i1 == 1 && i2 == 0);
    //bool cond_pp1 = (i1 == 1 && i2 == 1); //no test done for this

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

    int index1=-1; int index2=-1;

    auto it1 = std::find(sortVec1.begin(), sortVec1.end(), ku);
    if (it1 == sortVec1.end()){
        std::cout << "make2x2Block_inverse in MakeAmpMat, kUnique not found for rows" << std::endl;
    }
    else{
      index1 = distance(sortVec1.begin(), it1);
    }

    auto it2 = std::find(sortVec2.begin(), sortVec2.end(), ku);
    if (it2 == sortVec2.end()){
        std::cout << "make2x2Block_inverse in MakeAmpMat, kUnique not found for columns" << std::endl;
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
    //spp::sparse_hash_map<int, double> T2_temp;

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

                    //T2_elements_A.push_back( inMat(i-range_lower1,j-range_lower2) );
                    //T2_elements_I[id] = T2_elements_A.size()-1;
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

                    //T2_elements_A.push_back( inMat(i-range_lower1,j-range_lower2) );
                    //T2_elements_I[id] = T2_elements_A.size()-1;
                }
            }
        }
    }
}

void MakeAmpMat::make3x3Block_inverse(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, int i5, int i6, spp::sparse_hash_map<int, double> &T_list, bool add){

    bool cond_hhh1 = (i1 == 0 && i2 == 0 && i3==0);
    bool cond_pph1 = (i1 == 1 && i2 == 1 && i3==0);
    bool cond_hhp1 = (i1 == 0 && i2 == 0 && i3==1);
    //bool cond_ppp1 = (i1 == 1 && i2 == 1 && i3==1);

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

    int index1=-1; int index2=-1;

    auto it1 = std::find(sortVec1.begin(), sortVec1.end(), ku);
    if (it1 == sortVec1.end()){
        std::cout << "make3x1Block_inverse in MakeAmpMat, kUnique not found for rows" << std::endl;
    }
    else{
      index1 = distance(sortVec1.begin(), it1);
    }

    auto it2 = std::find(sortVec2.begin(), sortVec2.end(), ku);
    if (it2 == sortVec2.end()){
        std::cout << "make3x1Block_inverse in MakeAmpMat, kUnique not found for columns" << std::endl;
    }
    else{
      index2 = distance(sortVec2.begin(), it2);
    }

    int range_lower1 = indexHolder1_pointer(0,index1);
    int range_upper1 = indexHolder1_pointer(1,index1);
    int range_lower2 = indexHolder2_pointer(0,index2);
    int range_upper2 = indexHolder2_pointer(1,index2);

    unsigned long int id;
    int ii; int jj; int kk;
    int aa; int bb; int cc;

    double val;

    //std::cout << inMat << std::endl;
    if (add == true){
        if (cond_hhp1 && cond_pph2){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    val = inMat(i-range_lower1,j-range_lower2);
                    if (val != 0){
                        ii = (blockArrays1_pointer)(1,i);
                        jj = (blockArrays1_pointer)(2,i);
                        kk = (blockArrays2_pointer)(3,j);
                        aa = (blockArrays2_pointer)(1,j);
                        bb = (blockArrays2_pointer)(2,j);
                        cc = (blockArrays1_pointer)(3,i);
                        id = m_intClass->Identity_hhhppp(ii,jj,kk,aa,bb,cc);
                        T3_temp[id] = val;
                    }
                    //std::cout << T3_temp[id] << std::endl;
                }
            }
        }
        else if (cond_hhh1 && cond_ppp2){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    val = inMat(i-range_lower1,j-range_lower2);
                    if (val != 0){
                        ii = (blockArrays1_pointer)(1,i);
                        jj = (blockArrays1_pointer)(2,i);
                        kk = (blockArrays1_pointer)(3,i);
                        aa = (blockArrays2_pointer)(1,j);
                        bb = (blockArrays2_pointer)(2,j);
                        cc = (blockArrays2_pointer)(3,j);
                        id = m_intClass->Identity_hhhppp(ii,jj,kk,aa,bb,cc);
                        T3_temp[id] = val;
                    }
                }
            }
        }
    }
    else{
        if (cond_hhp1 && cond_pph2){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    val = inMat(i-range_lower1,j-range_lower2);
                    if (val != 0){
                        id = m_intClass->Identity_hhhppp((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays2_pointer)(3,j), (blockArrays2_pointer)(1,j), (blockArrays2_pointer)(2,j),  (blockArrays1_pointer)(3,i));
                        T_list[id] = val;
                    }
                }
            }
        }
        else if (cond_hhh1 && cond_ppp2){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    val = inMat(i-range_lower1,j-range_lower2);
                    if (val != 0){
                        id = m_intClass->Identity_hhhppp((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays1_pointer)(3,i), (blockArrays2_pointer)(1,j), (blockArrays2_pointer)(2,j), (blockArrays2_pointer)(3,j));
                        T_list[id] = val;
                    }
                }
            }
        }
    }
}

void MakeAmpMat::T3_makeDirectMat(){
    for (int channel = 0; channel<m_intClass->numOfKu3; channel++){
        int range_lower1 = m_intClass->boundsHolder_hhhppp_hhh(0,channel);
        int range_upper1 = m_intClass->boundsHolder_hhhppp_hhh(1,channel);
        int range_lower2 = m_intClass->boundsHolder_hhhppp_ppp(0,channel);
        int range_upper2 = m_intClass->boundsHolder_hhhppp_ppp(1,channel);

        Eigen::MatrixXi insertMat;
        insertMat.conservativeResize(range_upper1-range_lower1, range_upper2-range_lower2);

        int index; unsigned long int id;
        for (int ppp = range_lower2; ppp<range_upper2; ppp++){
            int aa = (m_intClass->blockArrays_ppp_ppp)(1,ppp);
            int bb = (m_intClass->blockArrays_ppp_ppp)(2,ppp);
            int cc = (m_intClass->blockArrays_ppp_ppp)(3,ppp);

            //int th = omp_get_thread_num();

            for (int hhh = range_lower1; hhh<range_upper1; hhh++){

                int ii = (m_intClass->blockArrays_ppp_hhh)(1,hhh);
                int jj = (m_intClass->blockArrays_ppp_hhh)(2,hhh);
                int kk = (m_intClass->blockArrays_ppp_hhh)(3,hhh);

                id = m_intClass->Identity_hhhppp(ii,jj,kk,aa,bb,cc);
                //id = ii + jj*N1 + kk*N2 + aa*N3 + bb*N4 + cc*N5;

                index = T3_elements_I.find(id)->second;
                //if (index == 0){std::cout << index << std::endl;}
                insertMat(hhh-range_lower1, ppp-range_lower2) = index;
            }
        }
        T3_directMat.push_back( insertMat );
    }
}

void MakeAmpMat::T3_makeMap(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, int i5, int i6){
    bool cond_hhh1 = (i1 == 0 && i2 == 0 && i3==0);
    bool cond_pph1 = (i1 == 1 && i2 == 1 && i3==0);
    bool cond_hhp1 = (i1 == 0 && i2 == 0 && i3==1);
    //bool cond_ppp1 = (i1 == 1 && i2 == 1 && i3==1);

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

    int index1=-1; int index2=-1;

    auto it1 = std::find(sortVec1.begin(), sortVec1.end(), ku);
    if (it1 == sortVec1.end()){
        std::cout << "make3x1Block_inverse in MakeAmpMat, kUnique not found for rows" << std::endl;
    }
    else{
      index1 = distance(sortVec1.begin(), it1);
    }

    auto it2 = std::find(sortVec2.begin(), sortVec2.end(), ku);
    if (it2 == sortVec2.end()){
        std::cout << "make3x1Block_inverse in MakeAmpMat, kUnique not found for columns" << std::endl;
    }
    else{
      index2 = distance(sortVec2.begin(), it2);
    }

    int range_lower1 = indexHolder1_pointer(0,index1);
    int range_upper1 = indexHolder1_pointer(1,index1);
    int range_lower2 = indexHolder2_pointer(0,index2);
    int range_upper2 = indexHolder2_pointer(1,index2);

    unsigned long int id; int index;
    int ii; int jj; int kk;
    int aa; int bb; int cc;

    double val;

    if (cond_hhp1 && cond_pph2){
        for (int i = range_lower1; i<range_upper1; i++){
            for (int j = range_lower2; j<range_upper2; j++){
                val = inMat(i-range_lower1,j-range_lower2);
                if (val != 0){
                    ii = (blockArrays1_pointer)(1,i);
                    jj = (blockArrays1_pointer)(2,i);
                    kk = (blockArrays2_pointer)(3,j);
                    aa = (blockArrays2_pointer)(1,j);
                    bb = (blockArrays2_pointer)(2,j);
                    cc = (blockArrays1_pointer)(3,i);
                    id = m_intClass->Identity_hhhppp(ii,jj,kk,aa,bb,cc);

                    //remember to remove first loops in make T3 in diagrams if you want to use this
                    //T3_elements_I.find(id)->second = index;
                    //T3_elements_A.push_back(val);

                    if (T3_elements_I.find(id) != T3_elements_I.end() ){
                        index = T3_elements_I.find(id)->second;
                        T3_elements_A[index] += val;
                    }
                }
            }
        }
    }
    else if (cond_hhh1 && cond_ppp2){
        for (int i = range_lower1; i<range_upper1; i++){
            for (int j = range_lower2; j<range_upper2; j++){
                val = inMat(i-range_lower1,j-range_lower2);
                if (val != 0){
                    ii = (blockArrays1_pointer)(1,i);
                    jj = (blockArrays1_pointer)(2,i);
                    kk = (blockArrays1_pointer)(3,i);
                    aa = (blockArrays2_pointer)(1,j);
                    bb = (blockArrays2_pointer)(2,j);
                    cc = (blockArrays2_pointer)(3,j);
                    id = m_intClass->Identity_hhhppp(ii,jj,kk,aa,bb,cc);

                    if (T3_elements_I.find(id) != T3_elements_I.end() ){
                        index = T3_elements_I.find(id)->second;
                        T3_elements_A[index] += val;
                    }
                }
            }
        }
    }
}

void MakeAmpMat::make3x3Block_inverse_I(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, int i5, int i6, std::vector<double> &T_vec, bool add){

    bool cond_hhh1 = (i1 == 0 && i2 == 0 && i3==0);
    bool cond_pph1 = (i1 == 1 && i2 == 1 && i3==0);
    bool cond_hhp1 = (i1 == 0 && i2 == 0 && i3==1);
    //bool cond_ppp1 = (i1 == 1 && i2 == 1 && i3==1);

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

    int index1=-1; int index2=-1;

    auto it1 = std::find(sortVec1.begin(), sortVec1.end(), ku);
    if (it1 == sortVec1.end()){
        std::cout << "make3x1Block_inverse in MakeAmpMat, kUnique not found for rows" << std::endl;
    }
    else{
      index1 = distance(sortVec1.begin(), it1);
    }

    auto it2 = std::find(sortVec2.begin(), sortVec2.end(), ku);
    if (it2 == sortVec2.end()){
        std::cout << "make3x1Block_inverse in MakeAmpMat, kUnique not found for columns" << std::endl;
    }
    else{
      index2 = distance(sortVec2.begin(), it2);
    }

    int range_lower1 = indexHolder1_pointer(0,index1);
    int range_upper1 = indexHolder1_pointer(1,index1);
    int range_lower2 = indexHolder2_pointer(0,index2);
    int range_upper2 = indexHolder2_pointer(1,index2);

    unsigned long int id; int index;
    int ii; int jj; int kk;
    int aa; int bb; int cc;

    double val;

    //std::cout << inMat << std::endl;
    if (add == true){
        if (cond_hhp1 && cond_pph2){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    val = inMat(i-range_lower1,j-range_lower2);
                    //if (val != 0){
                        ii = (blockArrays1_pointer)(1,i);
                        jj = (blockArrays1_pointer)(2,i);
                        kk = (blockArrays2_pointer)(3,j);
                        aa = (blockArrays2_pointer)(1,j);
                        bb = (blockArrays2_pointer)(2,j);
                        cc = (blockArrays1_pointer)(3,i);
                        id = m_intClass->Identity_hhhppp(ii,jj,kk,aa,bb,cc);

                        index = T3_elements_I.find(id)->second;

                        //std::cout << T3_temp[id] << std::endl;

                        T_vec[index] += val;
                        //std::cout << T_vec[index] <<" " << val << std::endl;
                    //}
                    //std::cout << T3_temp[id] << std::endl;
                }
            }
        }
        else if (cond_hhh1 && cond_ppp2){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    val = inMat(i-range_lower1,j-range_lower2);
                    //if (val != 0){
                        ii = (blockArrays1_pointer)(1,i);
                        jj = (blockArrays1_pointer)(2,i);
                        kk = (blockArrays1_pointer)(3,i);
                        aa = (blockArrays2_pointer)(1,j);
                        bb = (blockArrays2_pointer)(2,j);
                        cc = (blockArrays2_pointer)(3,j);
                        id = m_intClass->Identity_hhhppp(ii,jj,kk,aa,bb,cc);

                        index = T3_elements_I.find(id)->second;

                        //std::cout << id << std::endl;
                        //std::cout << index << std::endl;
                        T_vec[index] += val;
                    //}
                }
            }
        }
    }
    else{
        if (cond_hhp1 && cond_pph2){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    val = inMat(i-range_lower1,j-range_lower2);
                    id = m_intClass->Identity_hhhppp((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays2_pointer)(3,j), (blockArrays2_pointer)(1,j), (blockArrays2_pointer)(2,j),  (blockArrays1_pointer)(3,i));

                    index = T3_elements_I.find(id)->second;
                    T_vec[index] = val;
                }
            }
        }
        else if (cond_hhh1 && cond_ppp2){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    val = inMat(i-range_lower1,j-range_lower2);
                    id = m_intClass->Identity_hhhppp((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays1_pointer)(3,i), (blockArrays2_pointer)(1,j), (blockArrays2_pointer)(2,j), (blockArrays2_pointer)(3,j));

                    index = T3_elements_I.find(id)->second;
                    T_vec[index] = val;
                    //m_Counter ++;
                }
            }
        }
        //std::cout << counter << std::endl;
    }
}

void MakeAmpMat::make3x3Block_inverse_I_T1a(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, int i5, int i6, std::vector<double> &T_vec, bool add){

    bool cond_hhh1 = (i1 == 0 && i2 == 0 && i3==0);
    bool cond_pph1 = (i1 == 1 && i2 == 1 && i3==0);
    bool cond_hhp1 = (i1 == 0 && i2 == 0 && i3==1);
    //bool cond_ppp1 = (i1 == 1 && i2 == 1 && i3==1);

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

    int index1=-1; int index2=-1;

    auto it1 = std::find(sortVec1.begin(), sortVec1.end(), ku);
    if (it1 == sortVec1.end()){
        std::cout << "make3x1Block_inverse in MakeAmpMat, kUnique not found for rows" << std::endl;
    }
    else{
      index1 = distance(sortVec1.begin(), it1);
    }

    auto it2 = std::find(sortVec2.begin(), sortVec2.end(), ku);
    if (it2 == sortVec2.end()){
        std::cout << "make3x1Block_inverse in MakeAmpMat, kUnique not found for columns" << std::endl;
    }
    else{
      index2 = distance(sortVec2.begin(), it2);
    }

    int range_lower1 = indexHolder1_pointer(0,index1);
    int range_upper1 = indexHolder1_pointer(1,index1);
    int range_lower2 = indexHolder2_pointer(0,index2);
    int range_upper2 = indexHolder2_pointer(1,index2);

    unsigned long int id; int index;
    int ii; int jj; int kk;
    int aa; int bb; int cc;

    double val;

    //std::cout << inMat << std::endl;
    if (add == true){
        if (cond_hhp1 && cond_pph2){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    val = inMat(i-range_lower1,j-range_lower2);
                    //if (val != 0){
                        ii = (blockArrays1_pointer)(1,i);
                        jj = (blockArrays1_pointer)(2,i);
                        aa = (blockArrays1_pointer)(3,i);
                        bb = (blockArrays2_pointer)(1,j);
                        cc = (blockArrays2_pointer)(2,j);
                        kk = (blockArrays2_pointer)(3,j);
                        id = m_intClass->Identity_hhhppp(ii,jj,kk,aa,bb,cc);

                        index = T3_elements_I.find(id)->second;

                        //std::cout << T3_temp[id] << std::endl;

                        T_vec[index] += val;
                        //std::cout << T_vec[index] <<" " << val << std::endl;
                    //}
                    //std::cout << T3_temp[id] << std::endl;
                }
            }
        }
        else if (cond_hhh1 && cond_ppp2){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    val = inMat(i-range_lower1,j-range_lower2);
                    //if (val != 0){
                        ii = (blockArrays1_pointer)(1,i);
                        jj = (blockArrays1_pointer)(2,i);
                        kk = (blockArrays1_pointer)(3,i);
                        aa = (blockArrays2_pointer)(1,j);
                        bb = (blockArrays2_pointer)(2,j);
                        cc = (blockArrays2_pointer)(3,j);
                        id = m_intClass->Identity_hhhppp(ii,jj,kk,aa,bb,cc);

                        index = T3_elements_I.find(id)->second;
                        T_vec[index] += val;
                    //}
                }
            }
        }
    }
    else{
        if (cond_hhp1 && cond_pph2){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    val = inMat(i-range_lower1,j-range_lower2);
                    id = m_intClass->Identity_hhhppp((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays2_pointer)(3,j), (blockArrays2_pointer)(1,j), (blockArrays2_pointer)(2,j),  (blockArrays1_pointer)(3,i));

                    index = T3_elements_I.find(id)->second;
                    T_vec[index] = val;
                }
            }
        }
        else if (cond_hhh1 && cond_ppp2){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    val = inMat(i-range_lower1,j-range_lower2);
                    id = m_intClass->Identity_hhhppp((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays1_pointer)(3,i), (blockArrays2_pointer)(1,j), (blockArrays2_pointer)(2,j), (blockArrays2_pointer)(3,j));

                    index = T3_elements_I.find(id)->second;
                    T_vec[index] = val;
                    //m_Counter ++;
                }
            }
        }
        //std::cout << counter << std::endl;
    }
}

void MakeAmpMat::make3x3Block_inverse_I_T1b(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4, int i5, int i6, std::vector<double> &T_vec, bool add){

    bool cond_hhh1 = (i1 == 0 && i2 == 0 && i3==0);
    bool cond_pph1 = (i1 == 1 && i2 == 1 && i3==0);
    bool cond_hhp1 = (i1 == 0 && i2 == 0 && i3==1);
    //bool cond_ppp1 = (i1 == 1 && i2 == 1 && i3==1);

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

    int index1=-1; int index2=-1;

    auto it1 = std::find(sortVec1.begin(), sortVec1.end(), ku);
    if (it1 == sortVec1.end()){
        std::cout << "make3x1Block_inverse in MakeAmpMat, kUnique not found for rows" << std::endl;
    }
    else{
      index1 = distance(sortVec1.begin(), it1);
    }

    auto it2 = std::find(sortVec2.begin(), sortVec2.end(), ku);
    if (it2 == sortVec2.end()){
        std::cout << "make3x1Block_inverse in MakeAmpMat, kUnique not found for columns" << std::endl;
    }
    else{
      index2 = distance(sortVec2.begin(), it2);
    }

    int range_lower1 = indexHolder1_pointer(0,index1);
    int range_upper1 = indexHolder1_pointer(1,index1);
    int range_lower2 = indexHolder2_pointer(0,index2);
    int range_upper2 = indexHolder2_pointer(1,index2);

    unsigned long int id; int index;
    int ii; int jj; int kk;
    int aa; int bb; int cc;

    double val;

    //std::cout << inMat << std::endl;
    if (add == true){
        if (cond_hhp1 && cond_pph2){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    val = inMat(i-range_lower1,j-range_lower2);
                    //if (val != 0){
                        jj = (blockArrays1_pointer)(1,i);
                        kk = (blockArrays1_pointer)(2,i);
                        cc = (blockArrays1_pointer)(3,i);
                        aa = (blockArrays2_pointer)(1,j);
                        bb = (blockArrays2_pointer)(2,j);
                        ii = (blockArrays2_pointer)(3,j);
                        id = m_intClass->Identity_hhhppp(ii,jj,kk,aa,bb,cc);

                        index = T3_elements_I.find(id)->second;

                        //std::cout << T3_temp[id] << std::endl;

                        T_vec[index] += val;
                        //std::cout << T_vec[index] <<" " << val << std::endl;
                    //}
                    //std::cout << T3_temp[id] << std::endl;
                }
            }
        }
        else if (cond_hhh1 && cond_ppp2){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    val = inMat(i-range_lower1,j-range_lower2);
                    //if (val != 0){
                        ii = (blockArrays1_pointer)(1,i);
                        jj = (blockArrays1_pointer)(2,i);
                        kk = (blockArrays1_pointer)(3,i);
                        aa = (blockArrays2_pointer)(1,j);
                        bb = (blockArrays2_pointer)(2,j);
                        cc = (blockArrays2_pointer)(3,j);
                        id = m_intClass->Identity_hhhppp(ii,jj,kk,aa,bb,cc);

                        index = T3_elements_I.find(id)->second;
                        T_vec[index] += val;
                    //}
                }
            }
        }
    }
    else{
        if (cond_hhp1 && cond_pph2){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    val = inMat(i-range_lower1,j-range_lower2);
                    id = m_intClass->Identity_hhhppp((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays2_pointer)(3,j), (blockArrays2_pointer)(1,j), (blockArrays2_pointer)(2,j),  (blockArrays1_pointer)(3,i));

                    index = T3_elements_I.find(id)->second;
                    T_vec[index] = val;
                }
            }
        }
        else if (cond_hhh1 && cond_ppp2){
            for (int i = range_lower1; i<range_upper1; i++){
                for (int j = range_lower2; j<range_upper2; j++){
                    val = inMat(i-range_lower1,j-range_lower2);
                    id = m_intClass->Identity_hhhppp((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays1_pointer)(3,i), (blockArrays2_pointer)(1,j), (blockArrays2_pointer)(2,j), (blockArrays2_pointer)(3,j));

                    index = T3_elements_I.find(id)->second;
                    T_vec[index] = val;
                    //m_Counter ++;
                }
            }
        }
        //std::cout << counter << std::endl;
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


void MakeAmpMat::addElementsT3_T1a(){
    int ku;

    for (int channel = 0; channel<m_intClass->numOfKu3; channel++){
        ku = m_intClass->Vhhhppp_i[channel];

        int range_lower1 = m_intClass->boundsHolder_hhhppp_hhh(0,channel);
        int range_upper1 = m_intClass->boundsHolder_hhhppp_hhh(1,channel);
        int range_lower2 = m_intClass->boundsHolder_hhhppp_ppp(0,channel);
        int range_upper2 = m_intClass->boundsHolder_hhhppp_ppp(1,channel);

        Eigen::MatrixXi tempIMat = T3_directMat[channel]; //this holds indices for T3_elements_A (as well as _new and _temp)
        int rows = tempIMat.rows();
        int cols = tempIMat.cols();

        Eigen::MatrixXd tempAMat;
        tempAMat.conservativeResize( rows, cols );

        Eigen::MatrixXd tempAMat1;
        tempAMat1.conservativeResize( rows, cols );

        //should be easily parallelizable?
        for (int col=0; col<cols; col++){
            for (int row=0; row<rows; row++){
                tempAMat(row, col) = T3_elements_A_temp[ tempIMat(row, col) ];
            }
        }

        Eigen::MatrixXd Pab(rows,cols);
        Eigen::MatrixXd Pac(rows,cols);

        for (int i=range_lower2; i<range_upper2; i++){
            Pab.col( i-range_lower2 )  = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pab(1,i)-range_lower2 );
            Pac.col( i-range_lower2 )  = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pac(1,i)-range_lower2 );
        }

        tempAMat1 = tempAMat - Pab - Pac;

        Eigen::MatrixXd tempAMat2;
        tempAMat2.conservativeResize( rows, cols );

        Eigen::MatrixXd Pik(rows,cols);
        Eigen::MatrixXd Pjk(rows,cols);

        for (int i=range_lower1; i<range_upper1; i++){
            Pik.row( i - range_lower1)  = tempAMat1.row( m_intClass->blockArrays_ppp_hhh_Pik(1,i) - range_lower1 );
            Pjk.row( i - range_lower1 ) = tempAMat1.row( m_intClass->blockArrays_ppp_hhh_Pjk(1,i) - range_lower1 );
        }

        tempAMat2 = tempAMat1 -Pik - Pjk;

        for (int col=0; col<cols; col++){
            for (int row=0; row<rows; row++){
                T3_elements_A_new[ tempIMat(row, col) ] += tempAMat2(row, col);
            }
        }

        //make3x3Block_inverse_I(tempAMat2, ku, 0,0,0,1,1,1, T3_elements_A_new, true);
    }


    /*int Ns = m_intClass->m_Ns;
    int Nh = m_intClass->m_Nh;
    int N1 = Nh;
    int N2 = N1*Nh; //Nh*Nh;
    int N3 = N2*Nh; //Nh*Nh*Nh;
    int N4 = N3*Ns; //Nh*Nh*Nh*Np;
    int N5 = N4*Ns; //Nh*Nh*Nh*Np*Np;
    unsigned long int id; double val; int th; int index2;
    Eigen::MatrixXi mat_ID;
    Eigen::MatrixXd mat_VAL;
    Eigen::VectorXi vec_END;

    int size = m_intClass->blockArrays_ppp_hhh.cols()*m_intClass->blockArrays_ppp_ppp.cols();
    //std::cout << size << std::endl;
    int index = 0;

    int n = 1;//omp_get_max_threads();
    mat_ID.conservativeResize(n,size); mat_ID.setZero(n,size);
    mat_VAL.conservativeResize(n,size); mat_VAL.setZero(n,size);
    vec_END.conservativeResize(n); vec_END.setZero(n);

    //should use dynamic here
    //#pragma omp parallel for num_threads(n) private(N1,N2,N3,N4,N5, id, val, th, index2) firstprivate(index) //shared(index)
    for (int channel = 0; channel<m_intClass->numOfKu3; channel++){
        int range_lower1 = m_intClass->boundsHolder_hhhppp_hhh(0,channel);
        int range_upper1 = m_intClass->boundsHolder_hhhppp_hhh(1,channel);
        int range_lower2 = m_intClass->boundsHolder_hhhppp_ppp(0,channel);
        int range_upper2 = m_intClass->boundsHolder_hhhppp_ppp(1,channel);
        th =0;// omp_get_thread_num();

        //int size = (range_upper1-range_lower1)*(range_upper2-range_lower2);

        //int index = 0; //-1
        //#pragma omp parallel for num_threads(n) private(N1,N2,N3,N4,N5, id) firstprivate(index) //shared(index)
        for (int ppp = range_lower2; ppp<range_upper2; ppp++){
            int aa = (m_intClass->blockArrays_ppp_ppp)(1,ppp);
            int bb = (m_intClass->blockArrays_ppp_ppp)(2,ppp);
            int cc = (m_intClass->blockArrays_ppp_ppp)(3,ppp);


            for (int hhh = range_lower1; hhh<range_upper1; hhh++){

                //int id = 6015861;
                val = 0;

                int ii = (m_intClass->blockArrays_ppp_hhh)(1,hhh);
                int jj = (m_intClass->blockArrays_ppp_hhh)(2,hhh);
                int kk = (m_intClass->blockArrays_ppp_hhh)(3,hhh);

                //id = m_intClass->Identity_hhhppp(ii,jj,kk,bb,aa,cc);
                id = ii + jj*N1 + kk*N2 + bb*N3 + aa*N4 + cc*N5; // Pab
                //val -= T3_temp[id];
                if (T3_elements_I.find(id) != T3_elements_I.end() ){
                    index2 = T3_elements_I.find(id)->second;
                    val -= T3_elements_A_temp[index2];
                }

                //id = m_intClass->Identity_hhhppp(ii,jj,kk,cc,bb,aa);
                id = ii + jj*N1 + kk*N2 + cc*N3 + bb*N4 + aa*N5; // Pac
                //val -= T3_temp[id];
                if (T3_elements_I.find(id) != T3_elements_I.end() ){
                    index2 = T3_elements_I.find(id)->second;
                    val -= T3_elements_A_temp[index2];
                }

                //id = m_intClass->Identity_hhhppp(kk,jj,ii,aa,bb,cc);
                id = kk + jj*N1 + ii*N2 + aa*N3 + bb*N4 + cc*N5; // Pki
                //val -= T3_temp[id];
                if (T3_elements_I.find(id) != T3_elements_I.end() ){
                    index2 = T3_elements_I.find(id)->second;
                    val -= T3_elements_A_temp[index2];
                }

                //id = m_intClass->Identity_hhhppp(ii,kk,jj,aa,bb,cc);
                id = ii + kk*N1 + jj*N2 + aa*N3 + bb*N4 + cc*N5; // Pkj
                //val -= T3_temp[id];
                if (T3_elements_I.find(id) != T3_elements_I.end() ){
                    index2 = T3_elements_I.find(id)->second;
                    val -= T3_elements_A_temp[index2];
                }

                id = kk + jj*N1 + ii*N2 + bb*N3 + aa*N4 + cc*N5; // Pab Pki
                //val += T3_temp[id];
                if (T3_elements_I.find(id) != T3_elements_I.end() ){
                    index2 = T3_elements_I.find(id)->second;
                    val += T3_elements_A_temp[index2];
                }

                id = ii + kk*N1 + jj*N2 + bb*N3 + aa*N4 + cc*N5; // Pab Pkj
                //val += T3_temp[id];
                if (T3_elements_I.find(id) != T3_elements_I.end() ){
                    index2 = T3_elements_I.find(id)->second;
                    val += T3_elements_A_temp[index2];
                }

                id = kk + jj*N1 + ii*N2 + cc*N3 + bb*N4 + aa*N5; // Pac Pki
                //val += T3_temp[id];
                if (T3_elements_I.find(id) != T3_elements_I.end() ){
                    index2 = T3_elements_I.find(id)->second;
                    //std::cout << "sup" << std::endl;
                    val += T3_elements_A_temp[index2];
                }

                id = ii + kk*N1 + jj*N2 + cc*N3 + bb*N4 + aa*N5; // Pac Pkj
                //val += T3_temp[id];
                if (T3_elements_I.find(id) != T3_elements_I.end() ){
                    index2 = T3_elements_I.find(id)->second;
                    val += T3_elements_A_temp[index2];
                }

                id = ii + jj*N1 + kk*N2 + aa*N3 + bb*N4 + cc*N5; // direct
                //val += T3_temp[id];
                if (T3_elements_I.find(id) != T3_elements_I.end() ){
                    index2 = T3_elements_I.find(id)->second;
                    val += T3_elements_A_temp[index2];
                }
                //std::cout << T3_elements_I.size() << std::endl;

                if (val != 0){
                    //std::cout << val << std::endl;
                    mat_ID(th, index) = id;
                    mat_VAL(th, index) = val;
                    vec_END(th) = index+1;
                    index ++;
                }
            }
        }
    }

    for (int thread=0; thread<n; thread++){
        for (int index = 0; index<vec_END(thread); index++){
            //std::cout << th << " " << index << std::endl;
            val = mat_VAL(thread,index);
            id  = mat_ID(thread,index);
            index2 = T3_elements_I.find(id)->second;
            //std::cout << index2 << std::endl;
            //val -= T3_elements_A_temp[index2];
            T3_elements_A_new[index2] += val;
        }
    }*/
}

void MakeAmpMat::addElementsT3_T1b(){

    int ku;

    for (int channel = 0; channel<m_intClass->numOfKu3; channel++){
        ku = m_intClass->Vhhhppp_i[channel];

        int range_lower1 = m_intClass->boundsHolder_hhhppp_hhh(0,channel);
        int range_upper1 = m_intClass->boundsHolder_hhhppp_hhh(1,channel);
        int range_lower2 = m_intClass->boundsHolder_hhhppp_ppp(0,channel);
        int range_upper2 = m_intClass->boundsHolder_hhhppp_ppp(1,channel);

        Eigen::MatrixXi tempIMat = T3_directMat[channel]; //this holds indices for T3_elements_A (as well as _new and _temp)
        int rows = tempIMat.rows();
        int cols = tempIMat.cols();

        Eigen::MatrixXd tempAMat;
        tempAMat.conservativeResize( rows, cols );

        Eigen::MatrixXd tempAMat1;
        tempAMat1.conservativeResize( rows, cols );

        //should be easily parallizable?
        for (int col=0; col<cols; col++){
            for (int row=0; row<rows; row++){
                tempAMat(row, col) = T3_elements_A_temp[ tempIMat(row, col) ];
            }
        }

        Eigen::MatrixXd Pbc(rows,cols);
        Eigen::MatrixXd Pac(rows,cols);

        for (int i=range_lower2; i<range_upper2; i++){
            Pbc.col( i-range_lower2 )  = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pbc(1,i)-range_lower2 );
            Pac.col( i-range_lower2 )  = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pac(1,i)-range_lower2 );
        }

        tempAMat1 = tempAMat - Pbc - Pac;

        Eigen::MatrixXd tempAMat2;
        tempAMat2.conservativeResize( rows, cols );

        Eigen::MatrixXd Pik(rows,cols);
        Eigen::MatrixXd Pij(rows,cols);

        for (int i=range_lower1; i<range_upper1; i++){
            Pik.row( i - range_lower1)  = tempAMat1.row( m_intClass->blockArrays_ppp_hhh_Pik(1,i) - range_lower1 );
            Pij.row( i - range_lower1 )  = tempAMat1.row( m_intClass->blockArrays_ppp_hhh_Pij(1,i) - range_lower1 );
        }

        tempAMat2 = tempAMat1 -Pik - Pij;

        for (int col=0; col<cols; col++){
            for (int row=0; row<rows; row++){
                T3_elements_A_new[ tempIMat(row, col) ] += tempAMat2(row, col);
            }
        }

        //make3x3Block_inverse_I(tempAMat2, ku, 0,0,0,1,1,1, T3_elements_A_new, true);
    }


    /*int Np = m_intClass->m_Ns;
    int Nh = m_intClass->m_Nh;
    int N1 = Nh;
    int N2 = Nh*Nh; //Nh*Nh;
    int N3 = N2*Nh; //Nh*Nh*Nh;
    int N4 = N3*Np; //Nh*Nh*Nh*Np;
    int N5 = N4*Np; //Nh*Nh*Nh*Np*Np;
    unsigned long int id; double val; int th; int index2;
    Eigen::MatrixXi mat_ID;
    Eigen::MatrixXd mat_VAL;
    Eigen::VectorXi vec_END;

    int size = m_intClass->blockArrays_ppp_hhh.cols()*m_intClass->blockArrays_ppp_ppp.cols();
    int index = 0;

    int n = 1;//omp_get_max_threads();
    mat_ID.conservativeResize(n,size); mat_ID.setZero(n,size);
    mat_VAL.conservativeResize(n,size); mat_VAL.setZero(n,size);
    vec_END.conservativeResize(n); vec_END.setZero(n);


    int ii; int jj; int kk;
    int aa; int bb; int cc;

    //should use dynamic here
    //#pragma omp parallel for num_threads(n) private(N1,N2,N3,N4,N5, id, val, th, index2) firstprivate(index) //shared(index)
    for (int channel = 0; channel<m_intClass->numOfKu3; channel++){
        int range_lower1 = m_intClass->boundsHolder_hhhppp_hhh(0,channel);
        int range_upper1 = m_intClass->boundsHolder_hhhppp_hhh(1,channel);
        int range_lower2 = m_intClass->boundsHolder_hhhppp_ppp(0,channel);
        int range_upper2 = m_intClass->boundsHolder_hhhppp_ppp(1,channel);
        th =0;// omp_get_thread_num();

        //int size = (range_upper1-range_lower1)*(range_upper2-range_lower2);

        //int index = 0; //-1
        //#pragma omp parallel for num_threads(n) private(N1,N2,N3,N4,N5, id) firstprivate(index) //shared(index)
        for (int hhh = range_lower1; hhh<range_upper1; hhh++){
            for (int ppp = range_lower2; ppp<range_upper2; ppp++){

                ii = (m_intClass->blockArrays_ppp_hhh)(1,hhh);
                jj = (m_intClass->blockArrays_ppp_hhh)(2,hhh);
                kk = (m_intClass->blockArrays_ppp_hhh)(3,hhh);
                aa = (m_intClass->blockArrays_ppp_ppp)(1,ppp);
                bb = (m_intClass->blockArrays_ppp_ppp)(2,ppp);
                cc = (m_intClass->blockArrays_ppp_ppp)(3,ppp);

                val = 0;

                //id = m_intClass->Identity_hhhppp(ii,jj,kk,bb,aa,cc);
                id = ii + jj*N1 + kk*N2 + cc*N3 + bb*N4 + aa*N5; // Pca
                //val -= T3_temp[id];
                if (T3_elements_I.find(id) != T3_elements_I.end() ){
                    index2 = T3_elements_I.find(id)->second;
                    val -= T3_elements_A_temp[index2];
                }

                //id = m_intClass->Identity_hhhppp(ii,jj,kk,cc,bb,aa);
                id = ii + jj*N1 + kk*N2 + aa*N3 + cc*N4 + bb*N5; // Pcb
                //val -= T3_temp[id];
                if (T3_elements_I.find(id) != T3_elements_I.end() ){
                    index2 = T3_elements_I.find(id)->second;
                    val -= T3_elements_A_temp[index2];
                }

                //id = m_intClass->Identity_hhhppp(kk,jj,ii,aa,bb,cc);
                id = jj + ii*N1 + kk*N2 + aa*N3 + bb*N4 + cc*N5; // Pij
                //val -= T3_temp[id];
                if (T3_elements_I.find(id) != T3_elements_I.end() ){
                    index2 = T3_elements_I.find(id)->second;
                    val -= T3_elements_A_temp[index2];
                }

                //id = m_intClass->Identity_hhhppp(ii,kk,jj,aa,bb,cc);
                id = kk + jj*N1 + ii*N2 + aa*N3 + bb*N4 + cc*N5; // Pik
                //val -= T3_temp[id];
                if (T3_elements_I.find(id) != T3_elements_I.end() ){
                    index2 = T3_elements_I.find(id)->second;
                    val -= T3_elements_A_temp[index2];
                }

                id = jj + ii*N1 + kk*N2 + cc*N3 + bb*N4 + aa*N5; // Pca Pij
                //val += T3_temp[id];
                if (T3_elements_I.find(id) != T3_elements_I.end() ){
                    index2 = T3_elements_I.find(id)->second;
                    val += T3_elements_A_temp[index2];
                }

                id = kk + jj*N1 + ii*N2 + cc*N3 + bb*N4 + aa*N5; // Pca Pik
                //val += T3_temp[id];
                if (T3_elements_I.find(id) != T3_elements_I.end() ){
                    index2 = T3_elements_I.find(id)->second;
                    val += T3_elements_A_temp[index2];
                }

                id = jj + ii*N1 + kk*N2 + aa*N3 + cc*N4 + bb*N5; // Pcb Pij
                //val += T3_temp[id];
                if (T3_elements_I.find(id) != T3_elements_I.end() ){
                    index2 = T3_elements_I.find(id)->second;
                    val += T3_elements_A_temp[index2];
                }

                id = kk + jj*N1 + ii*N2 + aa*N3 + cc*N4 + bb*N5; // Pcb Pik
                //val += T3_temp[id];
                if (T3_elements_I.find(id) != T3_elements_I.end() ){
                    index2 = T3_elements_I.find(id)->second;
                    val += T3_elements_A_temp[index2];
                }

                id = ii + jj*N1 + kk*N2 + aa*N3 + bb*N4 + cc*N5; // direct
                //val += T3_temp[id];
                if (T3_elements_I.find(id) != T3_elements_I.end() ){
                    index2 = T3_elements_I.find(id)->second;
                    val += T3_elements_A_temp[index2];
                }

                if (val != 0){
                    mat_ID(th, index) = id;
                    mat_VAL(th, index) = val;
                    vec_END(th) = index+1;
                    index ++;
                }
            }
        }
    }

    for (int thread=0; thread<n; thread++){
        for (int index = 0; index<vec_END(thread); index++){
            //std::cout << th << " " << index << std::endl;
            val = mat_VAL(thread,index);
            id  = mat_ID(thread,index);
            index2 = T3_elements_I.find(id)->second;
            //std::cout << index2 << std::endl;
            //val -= T3_elements_A_temp[index2];
            T3_elements_A_new[index2] += val;
        }
    }*/
}

void MakeAmpMat::addElementsT3_T2c(){

    int ku;

    for (int channel = 0; channel<m_intClass->numOfKu3; channel++){
        ku = m_intClass->Vhhhppp_i[channel];

        int range_lower2 = m_intClass->boundsHolder_hhhppp_ppp(0,channel);
        int range_upper2 = m_intClass->boundsHolder_hhhppp_ppp(1,channel);

        Eigen::MatrixXi tempIMat = T3_directMat[channel]; //this holds indices for T3_elements_A (as well as _new and _temp)
        int rows = tempIMat.rows();
        int cols = tempIMat.cols();

        Eigen::MatrixXd tempAMat;
        tempAMat.conservativeResize( rows, cols );

        Eigen::MatrixXd tempAMat1;
        tempAMat1.conservativeResize( rows, cols );

        //should be easily parallizable?
        for (int col=0; col<cols; col++){
            for (int row=0; row<rows; row++){
                tempAMat(row, col) = T3_elements_A_temp[ tempIMat(row, col) ];
            }
        }

        //P(c/ab)
        Eigen::MatrixXd Pac(rows,cols);
        Eigen::MatrixXd Pbc(rows,cols);

        for (int i=range_lower2; i<range_upper2; i++){
            Pac.col( i-range_lower2 )   = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pac(1,i)-range_lower2 );
            Pbc.col( i-range_lower2 )   = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pbc(1,i)-range_lower2 );
        }

        tempAMat1 = tempAMat - Pac - Pbc;

        make3x3Block_inverse_I(tempAMat1, ku, 0,0,0,1,1,1, T3_elements_A_new, true);
    }
}

void MakeAmpMat::addElementsT3_T2d(){

    int ku;

    for (int channel = 0; channel<m_intClass->numOfKu3; channel++){
        ku = m_intClass->Vhhhppp_i[channel];

        int range_lower1 = m_intClass->boundsHolder_hhhppp_hhh(0,channel);
        int range_upper1 = m_intClass->boundsHolder_hhhppp_hhh(1,channel);

        Eigen::MatrixXi tempIMat = T3_directMat[channel]; //this holds indices for T3_elements_A (as well as _new and _temp)
        int rows = tempIMat.rows();
        int cols = tempIMat.cols();

        Eigen::MatrixXd tempAMat;
        tempAMat.conservativeResize( rows, cols );

        Eigen::MatrixXd tempAMat1;
        tempAMat1.conservativeResize( rows, cols );

        //should be easily parallizable?
        for (int col=0; col<cols; col++){
            for (int row=0; row<rows; row++){
                tempAMat(row, col) = T3_elements_A_temp[ tempIMat(row, col) ];
            }
        }

        //P(k/ij)
        Eigen::MatrixXd Pik(rows,cols);
        Eigen::MatrixXd Pjk(rows,cols);

        for (int i=range_lower1; i<range_upper1; i++){
            Pik.row( i-range_lower1 )   = tempAMat.row( m_intClass->blockArrays_ppp_hhh_Pik(1,i)-range_lower1 );
            Pjk.row( i-range_lower1 )   = tempAMat.row( m_intClass->blockArrays_ppp_hhh_Pjk(1,i)-range_lower1 );
        }

        tempAMat1 = tempAMat - Pik - Pjk;

        make3x3Block_inverse_I(tempAMat1, ku, 0,0,0,1,1,1, T3_elements_A_new, true);
    }
}

void MakeAmpMat::addElementsT3_T2e(){

    int ku;

    for (int channel = 0; channel<m_intClass->numOfKu3; channel++){
        ku = m_intClass->Vhhhppp_i[channel];

        int range_lower1 = m_intClass->boundsHolder_hhhppp_hhh(0,channel);
        int range_upper1 = m_intClass->boundsHolder_hhhppp_hhh(1,channel);
        int range_lower2 = m_intClass->boundsHolder_hhhppp_ppp(0,channel);
        int range_upper2 = m_intClass->boundsHolder_hhhppp_ppp(1,channel);

        Eigen::MatrixXi tempIMat = T3_directMat[channel]; //this holds indices for T3_elements_A (as well as _new and _temp)
        int rows = tempIMat.rows();
        int cols = tempIMat.cols();

        Eigen::MatrixXd tempAMat;
        tempAMat.conservativeResize( rows, cols );

        Eigen::MatrixXd tempAMat1;
        tempAMat1.conservativeResize( rows, cols );

        //should be easily parallizable?
        for (int col=0; col<cols; col++){
            for (int row=0; row<rows; row++){
                tempAMat(row, col) = T3_elements_A_temp[ tempIMat(row, col) ];
            }
        }

        //P(a/bc)
        Eigen::MatrixXd Pab(rows,cols);
        Eigen::MatrixXd Pac(rows,cols);

        for (int i=range_lower2; i<range_upper2; i++){
            Pab.col( i-range_lower2 )   = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pab(1,i)-range_lower2 );
            Pac.col( i-range_lower2 )   = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pac(1,i)-range_lower2 );
        }

        tempAMat1 = tempAMat - Pab - Pac;

        Eigen::MatrixXd tempAMat2;
        tempAMat2.conservativeResize( rows, cols );

        //P(i/jk)
        Eigen::MatrixXd Pij(rows,cols);
        Eigen::MatrixXd Pik(rows,cols);

        for (int i=range_lower1; i<range_upper1; i++){
            Pij.row( i - range_lower1 ) = tempAMat1.row( m_intClass->blockArrays_ppp_hhh_Pij(1,i) - range_lower1 );
            Pik.row( i - range_lower1)  = tempAMat1.row( m_intClass->blockArrays_ppp_hhh_Pik(1,i) - range_lower1 );
        }

        tempAMat2 = tempAMat1 - Pij - Pik;

        //std::cout << tempAMat2 << std::endl;

        make3x3Block_inverse_I(tempAMat2, ku, 0,0,0,1,1,1, T3_elements_A_new, true);
    }
}

void MakeAmpMat::addElementsT3_T3b(){

    int ku;

    for (int channel = 0; channel<m_intClass->numOfKu3; channel++){
        ku = m_intClass->Vhhhppp_i[channel];

        int range_lower1 = m_intClass->boundsHolder_hhhppp_hhh(0,channel);
        int range_upper1 = m_intClass->boundsHolder_hhhppp_hhh(1,channel);
        int range_lower2 = m_intClass->boundsHolder_hhhppp_ppp(0,channel);
        int range_upper2 = m_intClass->boundsHolder_hhhppp_ppp(1,channel);

        Eigen::MatrixXi tempIMat = T3_directMat[channel]; //this holds indices for T3_elements_A (as well as _new and _temp)
        int rows = tempIMat.rows();
        int cols = tempIMat.cols();

        Eigen::MatrixXd tempAMat;
        tempAMat.conservativeResize( rows, cols );

        Eigen::MatrixXd tempAMat1;
        tempAMat1.conservativeResize( rows, cols );

        //should be easily parallizable?
        for (int col=0; col<cols; col++){
            for (int row=0; row<rows; row++){
                tempAMat(row, col) = T3_elements_A_temp[ tempIMat(row, col) ];
            }
        }

        //P(abc)
        Eigen::MatrixXd Pab(rows,cols);
        Eigen::MatrixXd Pac(rows,cols);
        Eigen::MatrixXd Pbc(rows,cols);
        Eigen::MatrixXd Pabac(rows,cols);
        Eigen::MatrixXd Pabbc(rows,cols);

        for (int i=range_lower2; i<range_upper2; i++){
            Pab.col( i-range_lower2 )   = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pab(1,i)  - range_lower2 );
            Pac.col( i-range_lower2 )   = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pac(1,i)  - range_lower2 );
            Pbc.col( i-range_lower2 )   = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pbc(1,i)  - range_lower2 );
            Pabac.col( i-range_lower2 ) = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pabac(1,i)- range_lower2 );
            Pabbc.col( i-range_lower2 ) = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pabbc(1,i)- range_lower2 );
        }

        tempAMat1 = tempAMat - Pab - Pbc - Pac + Pabac + Pabbc;

        Eigen::MatrixXd tempAMat2;
        tempAMat2.conservativeResize( rows, cols );

        //P(i/jk)
        Eigen::MatrixXd Pik(rows,cols);
        Eigen::MatrixXd Pij(rows,cols);

        for (int i=range_lower1; i<range_upper1; i++){
            Pik.row( i - range_lower1)  = tempAMat1.row( m_intClass->blockArrays_ppp_hhh_Pik(1,i) - range_lower1 );
            Pij.row( i - range_lower1 ) = tempAMat1.row( m_intClass->blockArrays_ppp_hhh_Pij(1,i) - range_lower1 );
        }

        tempAMat2 = tempAMat1 - Pik - Pij;

        //std::cout << tempAMat2 << std::endl;

        make3x3Block_inverse_I(tempAMat2, ku, 0,0,0,1,1,1, T3_elements_A_new, true);
    }
}

void MakeAmpMat::addElementsT3_T3c(){

    int ku;

    for (int channel = 0; channel<m_intClass->numOfKu3; channel++){
        ku = m_intClass->Vhhhppp_i[channel];

        int range_lower1 = m_intClass->boundsHolder_hhhppp_hhh(0,channel);
        int range_upper1 = m_intClass->boundsHolder_hhhppp_hhh(1,channel);
        int range_lower2 = m_intClass->boundsHolder_hhhppp_ppp(0,channel);
        int range_upper2 = m_intClass->boundsHolder_hhhppp_ppp(1,channel);

        Eigen::MatrixXi tempIMat = T3_directMat[channel]; //this holds indices for T3_elements_A (as well as _new and _temp)
        int rows = tempIMat.rows();
        int cols = tempIMat.cols();

        Eigen::MatrixXd tempAMat;
        tempAMat.conservativeResize( rows, cols );

        Eigen::MatrixXd tempAMat1;
        tempAMat1.conservativeResize( rows, cols );

        //should be easily parallizable?
        for (int col=0; col<cols; col++){
            for (int row=0; row<rows; row++){
                tempAMat(row, col) = T3_elements_A_temp[ tempIMat(row, col) ];
            }
        }

        //P(a/bc)
        Eigen::MatrixXd Pab(rows,cols);
        Eigen::MatrixXd Pac(rows,cols);

        for (int i=range_lower2; i<range_upper2; i++){
            Pab.col( i-range_lower2 )   = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pab(1,i)-range_lower2 );
            Pac.col( i-range_lower2 )   = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pac(1,i)-range_lower2 );
        }

        tempAMat1 = tempAMat - Pab - Pac;

        Eigen::MatrixXd tempAMat2;
        tempAMat2.conservativeResize( rows, cols );

        //P(ijk)
        Eigen::MatrixXd Pik(rows,cols);
        Eigen::MatrixXd Pij(rows,cols);
        Eigen::MatrixXd Pjk(rows,cols);
        Eigen::MatrixXd Pijik(rows,cols);
        Eigen::MatrixXd Pijjk(rows,cols);

        for (int i=range_lower1; i<range_upper1; i++){
            Pik.row( i - range_lower1)  = tempAMat1.row( m_intClass->blockArrays_ppp_hhh_Pik(1,i) - range_lower1 );
            Pij.row( i - range_lower1 ) = tempAMat1.row( m_intClass->blockArrays_ppp_hhh_Pij(1,i) - range_lower1 );
            Pjk.row( i - range_lower1)  = tempAMat1.row( m_intClass->blockArrays_ppp_hhh_Pjk(1,i) - range_lower1 );
            Pijik.row( i - range_lower1 ) = tempAMat1.row( m_intClass->blockArrays_ppp_hhh_Pijik(1,i) - range_lower1 );
            Pijjk.row( i - range_lower1)  = tempAMat1.row( m_intClass->blockArrays_ppp_hhh_Pijjk(1,i) - range_lower1 );
        }

        tempAMat2 = tempAMat1 - Pik - Pij - Pjk + Pijik + Pijjk;

        //std::cout << tempAMat2 << std::endl;

        make3x3Block_inverse_I(tempAMat2, ku, 0,0,0,1,1,1, T3_elements_A_new, true);
    }
}

void MakeAmpMat::addElementsT3_T3d(){

    int ku;

    for (int channel = 0; channel<m_intClass->numOfKu3; channel++){
        ku = m_intClass->Vhhhppp_i[channel];

        int range_lower1 = m_intClass->boundsHolder_hhhppp_hhh(0,channel);
        int range_upper1 = m_intClass->boundsHolder_hhhppp_hhh(1,channel);
        int range_lower2 = m_intClass->boundsHolder_hhhppp_ppp(0,channel);
        int range_upper2 = m_intClass->boundsHolder_hhhppp_ppp(1,channel);

        Eigen::MatrixXi tempIMat = T3_directMat[channel]; //this holds indices for T3_elements_A (as well as _new and _temp)
        int rows = tempIMat.rows();
        int cols = tempIMat.cols();

        Eigen::MatrixXd tempAMat;
        tempAMat.conservativeResize( rows, cols );

        Eigen::MatrixXd tempAMat1;
        tempAMat1.conservativeResize( rows, cols );

        //should be easily parallizable?
        for (int col=0; col<cols; col++){
            for (int row=0; row<rows; row++){
                tempAMat(row, col) = T3_elements_A_temp[ tempIMat(row, col) ];
            }
        }

        //P(c/ab)
        Eigen::MatrixXd Pac(rows,cols);
        Eigen::MatrixXd Pbc(rows,cols);

        for (int i=range_lower2; i<range_upper2; i++){
            Pac.col( i-range_lower2 )   = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pac(1,i)-range_lower2 );
            Pbc.col( i-range_lower2 )   = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pbc(1,i)-range_lower2 );
        }

        tempAMat1 = tempAMat - Pac - Pbc;

        Eigen::MatrixXd tempAMat2;
        tempAMat2.conservativeResize( rows, cols );

        //P(i/jk)
        Eigen::MatrixXd Pij(rows,cols);
        Eigen::MatrixXd Pik(rows,cols);

        for (int i=range_lower1; i<range_upper1; i++){
            Pij.row( i - range_lower1 ) = tempAMat1.row( m_intClass->blockArrays_ppp_hhh_Pij(1,i) - range_lower1 );
            Pik.row( i - range_lower1)  = tempAMat1.row( m_intClass->blockArrays_ppp_hhh_Pik(1,i) - range_lower1 );
        }

        tempAMat2 = tempAMat1 - Pij - Pik;

        //std::cout << tempAMat2 << std::endl;

        make3x3Block_inverse_I(tempAMat2, ku, 0,0,0,1,1,1, T3_elements_A_new, true);
    }
}

void MakeAmpMat::addElementsT3_T3e(){

    int ku;

    for (int channel = 0; channel<m_intClass->numOfKu3; channel++){
        ku = m_intClass->Vhhhppp_i[channel];

        int range_lower1 = m_intClass->boundsHolder_hhhppp_hhh(0,channel);
        int range_upper1 = m_intClass->boundsHolder_hhhppp_hhh(1,channel);
        int range_lower2 = m_intClass->boundsHolder_hhhppp_ppp(0,channel);
        int range_upper2 = m_intClass->boundsHolder_hhhppp_ppp(1,channel);

        Eigen::MatrixXi tempIMat = T3_directMat[channel]; //this holds indices for T3_elements_A (as well as _new and _temp)
        int rows = tempIMat.rows();
        int cols = tempIMat.cols();

        Eigen::MatrixXd tempAMat;
        tempAMat.conservativeResize( rows, cols );

        Eigen::MatrixXd tempAMat1;
        tempAMat1.conservativeResize( rows, cols );

        //should be easily parallizable?
        for (int col=0; col<cols; col++){
            for (int row=0; row<rows; row++){
                tempAMat(row, col) = T3_elements_A_temp[ tempIMat(row, col) ];
            }
        }

        //P(a/bc)
        Eigen::MatrixXd Pab(rows,cols);
        Eigen::MatrixXd Pac(rows,cols);

        for (int i=range_lower2; i<range_upper2; i++){
            Pab.col( i-range_lower2 )   = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pab(1,i)-range_lower2 );
            Pac.col( i-range_lower2 )   = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pac(1,i)-range_lower2 );
        }

        tempAMat1 = tempAMat - Pab - Pac;

        Eigen::MatrixXd tempAMat2;
        tempAMat2.conservativeResize( rows, cols );

        //P(k/ij)
        Eigen::MatrixXd Pik(rows,cols);
        Eigen::MatrixXd Pjk(rows,cols);

        for (int i=range_lower1; i<range_upper1; i++){
            Pik.row( i - range_lower1 ) = tempAMat1.row( m_intClass->blockArrays_ppp_hhh_Pik(1,i) - range_lower1 );
            Pjk.row( i - range_lower1)  = tempAMat1.row( m_intClass->blockArrays_ppp_hhh_Pjk(1,i) - range_lower1 );
        }

        tempAMat2 = tempAMat1 - Pik - Pjk;

        //std::cout << tempAMat2 << std::endl;

        make3x3Block_inverse_I(tempAMat2, ku, 0,0,0,1,1,1, T3_elements_A_new, true);
    }
}

void MakeAmpMat::addElementsT3_T5a(){

    int ku;

    for (int channel = 0; channel<m_intClass->numOfKu3; channel++){
        ku = m_intClass->Vhhhppp_i[channel];

        int range_lower1 = m_intClass->boundsHolder_hhhppp_hhh(0,channel);
        int range_upper1 = m_intClass->boundsHolder_hhhppp_hhh(1,channel);
        int range_lower2 = m_intClass->boundsHolder_hhhppp_ppp(0,channel);
        int range_upper2 = m_intClass->boundsHolder_hhhppp_ppp(1,channel);

        Eigen::MatrixXi tempIMat = T3_directMat[channel]; //this holds indices for T3_elements_A (as well as _new and _temp)
        int rows = tempIMat.rows();
        int cols = tempIMat.cols();

        Eigen::MatrixXd tempAMat;
        tempAMat.conservativeResize( rows, cols );

        Eigen::MatrixXd tempAMat1;
        tempAMat1.conservativeResize( rows, cols );

        //should be easily parallizable?
        for (int col=0; col<cols; col++){
            for (int row=0; row<rows; row++){
                tempAMat(row, col) = T3_elements_A_temp[ tempIMat(row, col) ];
            }
        }

        //P(a/bc)
        Eigen::MatrixXd Pab(rows,cols);
        Eigen::MatrixXd Pac(rows,cols);

        for (int i=range_lower2; i<range_upper2; i++){
            Pab.col( i-range_lower2 )   = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pab(1,i)-range_lower2 );
            Pac.col( i-range_lower2 )   = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pac(1,i)-range_lower2 );
        }

        tempAMat1 = tempAMat - Pab - Pac;

        Eigen::MatrixXd tempAMat2;
        tempAMat2.conservativeResize( rows, cols );

        //P(i/jk)
        Eigen::MatrixXd Pij(rows,cols);
        Eigen::MatrixXd Pik(rows,cols);

        for (int i=range_lower1; i<range_upper1; i++){
            Pij.row( i - range_lower1)  = tempAMat1.row( m_intClass->blockArrays_ppp_hhh_Pij(1,i) - range_lower1 );
            Pik.row( i - range_lower1 ) = tempAMat1.row( m_intClass->blockArrays_ppp_hhh_Pik(1,i) - range_lower1 );
        }

        tempAMat2 = tempAMat1 - Pij - Pik;

        //std::cout << tempAMat2 << std::endl;

        make3x3Block_inverse_I(tempAMat2, ku, 0,0,0,1,1,1, T3_elements_A_new, true);
    }
}

void MakeAmpMat::addElementsT3_T5b(){

    int ku;

    for (int channel = 0; channel<m_intClass->numOfKu3; channel++){
        ku = m_intClass->Vhhhppp_i[channel];

        int range_lower1 = m_intClass->boundsHolder_hhhppp_hhh(0,channel);
        int range_upper1 = m_intClass->boundsHolder_hhhppp_hhh(1,channel);

        Eigen::MatrixXi tempIMat = T3_directMat[channel]; //this holds indices for T3_elements_A (as well as _new and _temp)
        int rows = tempIMat.rows();
        int cols = tempIMat.cols();

        Eigen::MatrixXd tempAMat;
        tempAMat.conservativeResize( rows, cols );

        Eigen::MatrixXd tempAMat1;
        tempAMat1.conservativeResize( rows, cols );

        //should be easily parallizable?
        for (int col=0; col<cols; col++){
            for (int row=0; row<rows; row++){
                tempAMat(row, col) = T3_elements_A_temp[ tempIMat(row, col) ];
            }
        }

        //P(i/jk)
        Eigen::MatrixXd Pij(rows,cols);
        Eigen::MatrixXd Pik(rows,cols);

        for (int i=range_lower1; i<range_upper1; i++){
            Pij.row( i - range_lower1)  = tempAMat.row( m_intClass->blockArrays_ppp_hhh_Pij(1,i) - range_lower1 );
            Pik.row( i - range_lower1 ) = tempAMat.row( m_intClass->blockArrays_ppp_hhh_Pik(1,i) - range_lower1 );
        }

        tempAMat1 = tempAMat - Pij - Pik;

        make3x3Block_inverse_I(tempAMat1, ku, 0,0,0,1,1,1, T3_elements_A_new, true);
    }
}

void MakeAmpMat::addElementsT3_T5c(){

    int ku;

    for (int channel = 0; channel<m_intClass->numOfKu3; channel++){
        ku = m_intClass->Vhhhppp_i[channel];

        int range_lower2 = m_intClass->boundsHolder_hhhppp_ppp(0,channel);
        int range_upper2 = m_intClass->boundsHolder_hhhppp_ppp(1,channel);

        Eigen::MatrixXi tempIMat = T3_directMat[channel]; //this holds indices for T3_elements_A (as well as _new and _temp)
        int rows = tempIMat.rows();
        int cols = tempIMat.cols();

        Eigen::MatrixXd tempAMat;
        tempAMat.conservativeResize( rows, cols );

        Eigen::MatrixXd tempAMat1;
        tempAMat1.conservativeResize( rows, cols );

        //should be easily parallizable?
        for (int col=0; col<cols; col++){
            for (int row=0; row<rows; row++){
                tempAMat(row, col) = T3_elements_A_temp[ tempIMat(row, col) ];
            }
        }

        //P(a/bc)
        Eigen::MatrixXd Pab(rows,cols);
        Eigen::MatrixXd Pac(rows,cols);

        for (int i=range_lower2; i<range_upper2; i++){
            Pab.col( i - range_lower2)  = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pab(1,i) - range_lower2 );
            Pac.col( i - range_lower2 ) = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pac(1,i) - range_lower2 );
        }

        tempAMat1 = tempAMat - Pab - Pac;

        make3x3Block_inverse_I(tempAMat1, ku, 0,0,0,1,1,1, T3_elements_A_new, true);
    }
}

void MakeAmpMat::addElementsT3_T5d(){

    int ku;

    for (int channel = 0; channel<m_intClass->numOfKu3; channel++){
        ku = m_intClass->Vhhhppp_i[channel];

        int range_lower1 = m_intClass->boundsHolder_hhhppp_hhh(0,channel);
        int range_upper1 = m_intClass->boundsHolder_hhhppp_hhh(1,channel);
        int range_lower2 = m_intClass->boundsHolder_hhhppp_ppp(0,channel);
        int range_upper2 = m_intClass->boundsHolder_hhhppp_ppp(1,channel);

        Eigen::MatrixXi tempIMat = T3_directMat[channel]; //this holds indices for T3_elements_A (as well as _new and _temp)
        int rows = tempIMat.rows();
        int cols = tempIMat.cols();

        Eigen::MatrixXd tempAMat;
        tempAMat.conservativeResize( rows, cols );

        Eigen::MatrixXd tempAMat1;
        tempAMat1.conservativeResize( rows, cols );

        //should be easily parallizable?
        for (int col=0; col<cols; col++){
            for (int row=0; row<rows; row++){
                tempAMat(row, col) = T3_elements_A_temp[ tempIMat(row, col) ];
            }
        }

        //P(a/bc)
        Eigen::MatrixXd Pab(rows,cols);
        Eigen::MatrixXd Pac(rows,cols);

        for (int i=range_lower2; i<range_upper2; i++){
            Pab.col( i-range_lower2 )   = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pab(1,i)-range_lower2 );
            Pac.col( i-range_lower2 )   = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pac(1,i)-range_lower2 );
        }

        tempAMat1 = tempAMat - Pab - Pac;

        Eigen::MatrixXd tempAMat2;
        tempAMat2.conservativeResize( rows, cols );

        //P(k/ij)
        Eigen::MatrixXd Pik(rows,cols);
        Eigen::MatrixXd Pjk(rows,cols);

        for (int i=range_lower1; i<range_upper1; i++){
            Pik.row( i - range_lower1 ) = tempAMat1.row( m_intClass->blockArrays_ppp_hhh_Pik(1,i) - range_lower1 );
            Pjk.row( i - range_lower1)  = tempAMat1.row( m_intClass->blockArrays_ppp_hhh_Pjk(1,i) - range_lower1 );
        }

        tempAMat2 = tempAMat1 - Pik - Pjk;

        make3x3Block_inverse_I(tempAMat2, ku, 0,0,0,1,1,1, T3_elements_A_new, true);
    }
}

void MakeAmpMat::addElementsT3_T5e(){

    int ku;

    for (int channel = 0; channel<m_intClass->numOfKu3; channel++){
        ku = m_intClass->Vhhhppp_i[channel];

        int range_lower1 = m_intClass->boundsHolder_hhhppp_hhh(0,channel);
        int range_upper1 = m_intClass->boundsHolder_hhhppp_hhh(1,channel);
        int range_lower2 = m_intClass->boundsHolder_hhhppp_ppp(0,channel);
        int range_upper2 = m_intClass->boundsHolder_hhhppp_ppp(1,channel);

        Eigen::MatrixXi tempIMat = T3_directMat[channel]; //this holds indices for T3_elements_A (as well as _new and _temp)
        int rows = tempIMat.rows();
        int cols = tempIMat.cols();

        Eigen::MatrixXd tempAMat;
        tempAMat.conservativeResize( rows, cols );

        Eigen::MatrixXd tempAMat1;
        tempAMat1.conservativeResize( rows, cols );

        //should be easily parallizable?
        for (int col=0; col<cols; col++){
            for (int row=0; row<rows; row++){
                tempAMat(row, col) = T3_elements_A_temp[ tempIMat(row, col) ];
            }
        }

        //P(c/ab)
        Eigen::MatrixXd Pac(rows,cols);
        Eigen::MatrixXd Pbc(rows,cols);

        for (int i=range_lower2; i<range_upper2; i++){
            Pac.col( i-range_lower2 )   = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pac(1,i)-range_lower2 );
            Pbc.col( i-range_lower2 )   = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pbc(1,i)-range_lower2 );
        }

        tempAMat1 = tempAMat - Pac - Pbc;

        Eigen::MatrixXd tempAMat2;
        tempAMat2.conservativeResize( rows, cols );

        //P(i/jk)
        Eigen::MatrixXd Pij(rows,cols);
        Eigen::MatrixXd Pik(rows,cols);

        for (int i=range_lower1; i<range_upper1; i++){
            Pij.row( i - range_lower1)  = tempAMat1.row( m_intClass->blockArrays_ppp_hhh_Pij(1,i) - range_lower1 );
            Pik.row( i - range_lower1 ) = tempAMat1.row( m_intClass->blockArrays_ppp_hhh_Pik(1,i) - range_lower1 );
        }

        tempAMat2 = tempAMat1 - Pij - Pik;

        make3x3Block_inverse_I(tempAMat2, ku, 0,0,0,1,1,1, T3_elements_A_new, true);
    }
}

void MakeAmpMat::addElementsT3_T5f(){

    int ku;

    for (int channel = 0; channel<m_intClass->numOfKu3; channel++){
        ku = m_intClass->Vhhhppp_i[channel];

        int range_lower1 = m_intClass->boundsHolder_hhhppp_hhh(0,channel);
        int range_upper1 = m_intClass->boundsHolder_hhhppp_hhh(1,channel);

        Eigen::MatrixXi tempIMat = T3_directMat[channel]; //this holds indices for T3_elements_A (as well as _new and _temp)
        int rows = tempIMat.rows();
        int cols = tempIMat.cols();

        Eigen::MatrixXd tempAMat;
        tempAMat.conservativeResize( rows, cols );

        Eigen::MatrixXd tempAMat1;
        tempAMat1.conservativeResize( rows, cols );

        //should be easily parallizable?
        for (int col=0; col<cols; col++){
            for (int row=0; row<rows; row++){
                tempAMat(row, col) = T3_elements_A_temp[ tempIMat(row, col) ];
            }
        }

        //P(k/ij)
        Eigen::MatrixXd Pik(rows,cols);
        Eigen::MatrixXd Pjk(rows,cols);

        for (int i=range_lower1; i<range_upper1; i++){
            Pik.row( i - range_lower1 ) = tempAMat.row( m_intClass->blockArrays_ppp_hhh_Pik(1,i) - range_lower1 );
            Pjk.row( i - range_lower1)  = tempAMat.row( m_intClass->blockArrays_ppp_hhh_Pjk(1,i) - range_lower1 );
        }

        tempAMat1 = tempAMat - Pik - Pjk;

        make3x3Block_inverse_I(tempAMat1, ku, 0,0,0,1,1,1, T3_elements_A_new, true);
    }
}

void MakeAmpMat::addElementsT3_T5g(){

    int ku;

    for (int channel = 0; channel<m_intClass->numOfKu3; channel++){
        ku = m_intClass->Vhhhppp_i[channel];

        int range_lower2 = m_intClass->boundsHolder_hhhppp_ppp(0,channel);
        int range_upper2 = m_intClass->boundsHolder_hhhppp_ppp(1,channel);

        Eigen::MatrixXi tempIMat = T3_directMat[channel]; //this holds indices for T3_elements_A (as well as _new and _temp)
        int rows = tempIMat.rows();
        int cols = tempIMat.cols();

        Eigen::MatrixXd tempAMat;
        tempAMat.conservativeResize( rows, cols );

        Eigen::MatrixXd tempAMat1;
        tempAMat1.conservativeResize( rows, cols );

        //should be easily parallizable?
        for (int col=0; col<cols; col++){
            for (int row=0; row<rows; row++){
                tempAMat(row, col) = T3_elements_A_temp[ tempIMat(row, col) ];
            }
        }

        //P(c/ab)
        Eigen::MatrixXd Pac(rows,cols);
        Eigen::MatrixXd Pbc(rows,cols);

        for (int i=range_lower2; i<range_upper2; i++){
            Pac.col( i-range_lower2 )   = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pac(1,i)-range_lower2 );
            Pbc.col( i-range_lower2 )   = tempAMat.col( m_intClass->blockArrays_ppp_ppp_Pbc(1,i)-range_lower2 );
        }

        tempAMat1 = tempAMat - Pac - Pbc;

        make3x3Block_inverse_I(tempAMat1, ku, 0,0,0,1,1,1, T3_elements_A_new, true);
    }
}

Eigen::MatrixXd MakeAmpMat::make3x3Block_I(int ku, int i1, int i2, int i3, int i4, int i5, int i6, std::vector<double> &T_vec){

    bool cond_hhh1 = (i1 == 0 && i2 == 0 && i3==0);
    bool cond_pph1 = (i1 == 1 && i2 == 1 && i3==0);
    bool cond_hhp1 = (i1 == 0 && i2 == 0 && i3==1);
    //bool cond_ppp1 = (i1 == 1 && i2 == 1 && i3==1);

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
    unsigned long int id; int index;
    if (cond_hhp1 && cond_pph2){
        for (int i = range_lower1; i<range_upper1; i++){
            for (int j = range_lower2; j<range_upper2; j++){
                id = m_intClass->Identity_hhhppp((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays2_pointer)(3,j), (blockArrays2_pointer)(1,j), (blockArrays2_pointer)(2,j), (blockArrays1_pointer)(3,i));
                index = T3_elements_I.find(id)->second;
                returnMat(i-range_lower1, j-range_lower2) = T_vec[index];
            }
        }
    }
    else if (cond_hhh1 && cond_ppp2){
        for (int i = range_lower1; i<range_upper1; i++){
            for (int j = range_lower2; j<range_upper2; j++){
                id = m_intClass->Identity_hhhppp((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays1_pointer)(3,i), (blockArrays2_pointer)(1,j), (blockArrays2_pointer)(2,j), (blockArrays2_pointer)(3,j));
                index = T3_elements_I.find(id)->second;
                returnMat(i-range_lower1, j-range_lower2) = T_vec[index];
            }
        }
    }
    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::I1_makemat_1(int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_pp_hh(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pp_hh(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pp_pp(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pp_pp(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int i; int j;
    int c; int d;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            /*
            i = m_intClass->blockArrays_pp_hh(1,i1);
            j = m_intClass->blockArrays_pp_hh(2,i1);
            c = m_intClass->blockArrays_pp_pp(1,i2);
            d = m_intClass->blockArrays_pp_pp(2,i2);
            id = m_intClass->Identity_hhpp(i,j,c,d);
            */
            id = m_intClass->Identity_hhpp(m_intClass->blockArrays_pp_hh(1,i1),
                                           m_intClass->blockArrays_pp_hh(2,i1),
                                           m_intClass->blockArrays_pp_pp(1,i2),
                                           m_intClass->blockArrays_pp_pp(2,i2));
            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::I1_makemat_2(int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_pp_hh(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pp_hh(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pp_pp(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pp_pp(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int k; int l;
    int a; int b;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            /*
            k = m_intClass->blockArrays_pp_hh(1,i1);
            l = m_intClass->blockArrays_pp_hh(2,i1);
            a = m_intClass->blockArrays_pp_pp(1,i2);
            b = m_intClass->blockArrays_pp_pp(2,i2);
            id = m_intClass->Identity_hhpp(k,l,a,b);
            */

            id = m_intClass->Identity_hhpp(m_intClass->blockArrays_pp_hh(1,i1),
                                           m_intClass->blockArrays_pp_hh(2,i1),
                                           m_intClass->blockArrays_pp_pp(1,i2),
                                           m_intClass->blockArrays_pp_pp(2,i2));

            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::I2_makemat_1(int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_pm_hp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pm_hp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pm_ph(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pm_ph(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int j; int l;
    int d; int b;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            /*
            j = m_intClass->blockArrays_pm_hp(1,i1);
            b = m_intClass->blockArrays_pm_hp(2,i1);
            d = m_intClass->blockArrays_pm_ph(1,i2);
            l = m_intClass->blockArrays_pm_ph(2,i2);
            id = m_intClass->Identity_hhpp(j,l,d,b);
            */

            id = m_intClass->Identity_hhpp(m_intClass->blockArrays_pm_hp(1,i1),
                                           m_intClass->blockArrays_pm_ph(2,i2),
                                           m_intClass->blockArrays_pm_ph(1,i2),
                                           m_intClass->blockArrays_pm_hp(2,i1));

            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::I2_makemat_2(int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_pm_hp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pm_hp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pm_ph(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pm_ph(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int i; int k;
    int c; int a;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            /*
            k = m_intClass->blockArrays_pm_hp(1,i1);
            c = m_intClass->blockArrays_pm_hp(2,i1);
            a = m_intClass->blockArrays_pm_ph(1,i2);
            i = m_intClass->blockArrays_pm_ph(2,i2);
            id = m_intClass->Identity_hhpp(i,k,c,a);
            */

            id = m_intClass->Identity_hhpp(m_intClass->blockArrays_pm_ph(2,i2),
                                           m_intClass->blockArrays_pm_hp(1,i1),
                                           m_intClass->blockArrays_pm_hp(2,i1),
                                           m_intClass->blockArrays_pm_ph(1,i2));

            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::I3_makemat_1(int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_ppm_pph(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_pph(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_h(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_h(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int j; int l;
    int a; int b;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            /*
            a = m_intClass->blockArrays_ppm_pph(1,i1);
            b = m_intClass->blockArrays_ppm_pph(2,i1);
            j = m_intClass->blockArrays_ppm_pph(3,i1);
            l = m_intClass->blockArrays_p_h(1,i2);
            id = m_intClass->Identity_hhpp(j,l,a,b);
            */

            id = m_intClass->Identity_hhpp(m_intClass->blockArrays_ppm_pph(3,i1),
                                           m_intClass->blockArrays_p_h(1,i2),
                                           m_intClass->blockArrays_ppm_pph(1,i1),
                                           m_intClass->blockArrays_ppm_pph(2,i1));

            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::I3_makemat_2(int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_ppm_pph(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_pph(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_h(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_h(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int i; int k;
    int c; int d;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            /*
            c = m_intClass->blockArrays_ppm_pph(1,i1);
            d = m_intClass->blockArrays_ppm_pph(2,i1);
            k = m_intClass->blockArrays_ppm_pph(3,i1);
            i = m_intClass->blockArrays_p_h(1,i2);
            id = m_intClass->Identity_hhpp(i,k,c,d);
            */

            id = m_intClass->Identity_hhpp(m_intClass->blockArrays_p_h(1,i2),
                                           m_intClass->blockArrays_ppm_pph(3,i1),
                                           m_intClass->blockArrays_ppm_pph(1,i1),
                                           m_intClass->blockArrays_ppm_pph(2,i1));

            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::I4_makemat_1(int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_p(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_p(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int k; int l;
    int c; int a;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            /*
            k = m_intClass->blockArrays_ppm_hhp(1,i1);
            l = m_intClass->blockArrays_ppm_hhp(2,i1);
            c = m_intClass->blockArrays_ppm_hhp(3,i1);
            a = m_intClass->blockArrays_p_p(1,i2);
            id = m_intClass->Identity_hhpp(k,l,c,a);
            */

            id = m_intClass->Identity_hhpp(m_intClass->blockArrays_ppm_hhp(1,i1),
                                           m_intClass->blockArrays_ppm_hhp(2,i1),
                                           m_intClass->blockArrays_ppm_hhp(3,i1),
                                           m_intClass->blockArrays_p_p(1,i2));

            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::I4_makemat_2(int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_p(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_p(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int i; int j;
    int d; int b;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            /*
            i = m_intClass->blockArrays_ppm_hhp(1,i1);
            j = m_intClass->blockArrays_ppm_hhp(2,i1);
            b = m_intClass->blockArrays_ppm_hhp(3,i1);
            d = m_intClass->blockArrays_p_p(1,i2);
            id = m_intClass->Identity_hhpp(i,j,d,b);
            */

            id = m_intClass->Identity_hhpp(m_intClass->blockArrays_ppm_hhp(1,i1),
                                           m_intClass->blockArrays_ppm_hhp(2,i1),
                                           m_intClass->blockArrays_p_p(1,i2),
                                           m_intClass->blockArrays_ppm_hhp(3,i1));

            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

void MakeAmpMat::I1_inverse(Eigen::MatrixXd &inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_pp_hh(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pp_hh(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pp_pp(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pp_pp(1, channel2);


    int id;
    int i; int j;
    int a; int b;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            /*
            i = m_intClass->blockArrays_pp_hh(1,i1);
            j = m_intClass->blockArrays_pp_hh(2,i1);
            a = m_intClass->blockArrays_pp_pp(1,i2);
            b = m_intClass->blockArrays_pp_pp(2,i2);
            id = m_intClass->Identity_hhpp(i,j,a,b);
            */

            id = m_intClass->Identity_hhpp(m_intClass->blockArrays_pp_hh(1,i1),
                                           m_intClass->blockArrays_pp_hh(2,i1),
                                           m_intClass->blockArrays_pp_pp(1,i2),
                                           m_intClass->blockArrays_pp_pp(2,i2));

            T2_temp[id] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}

void MakeAmpMat::I2_inverse(Eigen::MatrixXd &inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_pm_hp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pm_hp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pm_ph(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pm_ph(1, channel2);


    int id;
    int i; int j;
    int a; int b;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            /*
            j = m_intClass->blockArrays_pm_hp(1,i1);
            b = m_intClass->blockArrays_pm_hp(2,i1);
            a = m_intClass->blockArrays_pm_ph(1,i2);
            i = m_intClass->blockArrays_pm_ph(2,i2);
            id = m_intClass->Identity_hhpp(i,j,a,b);
            */

            id = m_intClass->Identity_hhpp(m_intClass->blockArrays_pm_ph(2,i2),
                                           m_intClass->blockArrays_pm_hp(1,i1),
                                           m_intClass->blockArrays_pm_ph(1,i2),
                                           m_intClass->blockArrays_pm_hp(2,i1));

            T2_temp[id] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}

void MakeAmpMat::I3_inverse(Eigen::MatrixXd inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_ppm_pph(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_pph(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_h(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_h(1, channel2);


    int id;
    int i; int j;
    int a; int b;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            /*
            a = m_intClass->blockArrays_ppm_pph(1,i1);
            b = m_intClass->blockArrays_ppm_pph(2,i1);
            j = m_intClass->blockArrays_ppm_pph(3,i1);
            i = m_intClass->blockArrays_p_h(1,i2);
            id = m_intClass->Identity_hhpp(i,j,a,b);
            */

            id = m_intClass->Identity_hhpp(m_intClass->blockArrays_p_h(1,i2),
                                           m_intClass->blockArrays_ppm_pph(3,i1),
                                           m_intClass->blockArrays_ppm_pph(1,i1),
                                           m_intClass->blockArrays_ppm_pph(2,i1));

            T2_temp[id] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}

void MakeAmpMat::I4_inverse(Eigen::MatrixXd inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_p(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_p(1, channel2);


    int id;
    int i; int j;
    int a; int b;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            /*
            i = m_intClass->blockArrays_ppm_hhp(1,i1);
            j = m_intClass->blockArrays_ppm_hhp(2,i1);
            b = m_intClass->blockArrays_ppm_hhp(3,i1);
            a = m_intClass->blockArrays_p_p(1,i2);
            id = m_intClass->Identity_hhpp(i,j,a,b);
            */

            id = m_intClass->Identity_hhpp(m_intClass->blockArrays_ppm_hhp(1,i1),
                                           m_intClass->blockArrays_ppm_hhp(2,i1),
                                           m_intClass->blockArrays_p_p(1,i2),
                                           m_intClass->blockArrays_ppm_hhp(3,i1));

            T2_temp[id] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}

Eigen::MatrixXd MakeAmpMat::D10b_makemat(int channel1, int channel2){    //makes a 3x3 matrix

    int range_lower1 = m_intClass->indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_ppm_pph(0, channel2);
    int range_upper2 = m_intClass->indexHolder_ppm_pph(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    unsigned long int id; int index;
    int i; int j; int k;
    int c; int d; int b;

    for (int i2 = range_lower2; i2<range_upper2; i2++){
        c = m_intClass->blockArrays_ppm_pph(1,i2);
        d = m_intClass->blockArrays_ppm_pph(2,i2);
        k = m_intClass->blockArrays_ppm_pph(3,i2);
        for (int i1 = range_lower1; i1<range_upper1; i1++){
            i = m_intClass->blockArrays_ppm_hhp(1,i1);
            j = m_intClass->blockArrays_ppm_hhp(2,i1);
            b = m_intClass->blockArrays_ppm_hhp(3,i1);

            id = m_intClass->Identity_hhhppp(i,j,k,c,d,b);

            /*id = m_intClass->Identity_hhhppp(m_intClass->blockArrays_ppm_hhp(1,i1),
                                             m_intClass->blockArrays_ppm_hhp(2,i1),
                                             m_intClass->blockArrays_ppm_pph(3,i2),
                                             m_intClass->blockArrays_ppm_pph(1,i2),
                                             m_intClass->blockArrays_ppm_pph(2,i2),
                                             m_intClass->blockArrays_ppm_hhp(3,i1));*/

            index = T3_elements_I.find(id)->second;//T3_elements_IV[thread][id];
            returnMat(i1-range_lower1, i2-range_lower2) = T3_elements_A[index];
        }
    }

    return returnMat;
}

void MakeAmpMat::D10b_inverse(Eigen::MatrixXd &inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_p(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_p(1, channel2);


    int id;
    int i; int j;
    int a; int b;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        i = m_intClass->blockArrays_ppm_hhp(1,i1);
        j = m_intClass->blockArrays_ppm_hhp(2,i1);
        b = m_intClass->blockArrays_ppm_hhp(3,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            a = m_intClass->blockArrays_p_p(1,i2);

            id = m_intClass->Identity_hhpp(i,j,a,b);
            T2_temp[id] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}

Eigen::MatrixXd MakeAmpMat::D10c_makemat(int channel1, int channel2){    //makes a 3x3 matrix

    int range_lower1 = m_intClass->indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_ppm_pph(0, channel2);
    int range_upper2 = m_intClass->indexHolder_ppm_pph(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    unsigned long int id; int index;
    int i; int k; int l;
    int c; int a; int b;

    for (int i2 = range_lower2; i2<range_upper2; i2++){
        a = m_intClass->blockArrays_ppm_pph(1,i2);
        b = m_intClass->blockArrays_ppm_pph(2,i2);
        i = m_intClass->blockArrays_ppm_pph(3,i2);
        for (int i1 = range_lower1; i1<range_upper1; i1++){

            k = m_intClass->blockArrays_ppm_hhp(1,i1);
            l = m_intClass->blockArrays_ppm_hhp(2,i1);
            c = m_intClass->blockArrays_ppm_hhp(3,i1);

            id = m_intClass->Identity_hhhppp(i,k,l,c,a,b);

            /*id = m_intClass->Identity_hhhppp(m_intClass->blockArrays_ppm_pph(3,i2),
                                             m_intClass->blockArrays_ppm_hhp(1,i1),
                                             m_intClass->blockArrays_ppm_hhp(2,i1),
                                             m_intClass->blockArrays_ppm_hhp(3,i1),
                                             m_intClass->blockArrays_ppm_pph(1,i2),
                                             m_intClass->blockArrays_ppm_pph(2,i2));*/

            index = T3_elements_I.find(id)->second;//T3_elements_IV[thread][id];
            returnMat(i1-range_lower1, i2-range_lower2) = T3_elements_A[index];
        }
    }

    return returnMat;
}

void MakeAmpMat::D10c_inverse(Eigen::MatrixXd &inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_ppm_pph(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_pph(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_h(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_h(1, channel2);


    int id;
    int i; int j;
    int a; int b;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            a = m_intClass->blockArrays_ppm_pph(1,i1);
            b = m_intClass->blockArrays_ppm_pph(2,i1);
            i = m_intClass->blockArrays_ppm_pph(3,i1);
            j = m_intClass->blockArrays_p_h(1,i2);

            id = m_intClass->Identity_hhpp(i,j,a,b);
            T2_temp[id] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}

Eigen::MatrixXd MakeAmpMat::T1a_makemat(int channel1, int channel2){    //makes a 3x1 matrix

    int range_lower1 = m_intClass->indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_p(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_p(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int i; int j;
    int a; int d;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            i = m_intClass->blockArrays_ppm_hhp(1,i1);
            j = m_intClass->blockArrays_ppm_hhp(2,i1);
            a = m_intClass->blockArrays_ppm_hhp(3,i1);
            d = m_intClass->blockArrays_p_p(1,i2);

            id = m_intClass->Identity_hhpp(i,j,a,d);
            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T1b_makemat(int channel1, int channel2){    //makes a 3x1 matrix

    int range_lower1 = m_intClass->indexHolder_ppm_pph(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_pph(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_h(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_h(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int i; int l;
    int a; int b;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            a = m_intClass->blockArrays_ppm_pph(1,i1);
            b = m_intClass->blockArrays_ppm_pph(2,i1);
            i = m_intClass->blockArrays_ppm_pph(3,i1);
            l = m_intClass->blockArrays_p_h(1,i2);

            id = m_intClass->Identity_hhpp(i,l,a,b);
            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T2c_makemat(int channel1, int channel2){    //makes a 4x2 matrix
    //std::cout << "sup" << std::endl;
    int range_lower1 = m_intClass->indexHolder_pppm_hhhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pppm_hhhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pp_pp(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pp_pp(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);

    unsigned long int id; int index;
    int i; int j; int k;
    int d; int e; int c;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        i = m_intClass->blockArrays_pppm_hhhp(1,i1);
        j = m_intClass->blockArrays_pppm_hhhp(2,i1);
        k = m_intClass->blockArrays_pppm_hhhp(3,i1);
        c = m_intClass->blockArrays_pppm_hhhp(4,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            d = m_intClass->blockArrays_pp_pp(1,i2);
            e = m_intClass->blockArrays_pp_pp(2,i2);

            id = m_intClass->Identity_hhhppp(i,j,k,d,e,c);

            index = T3_elements_I.find(id)->second;
            returnMat(i1-range_lower1, i2-range_lower2) = T3_elements_A[index];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T2d_makemat(int channel1, int channel2){    //makes a 4x2 matrix

    int range_lower1 = m_intClass->indexHolder_pppm_ppph(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pppm_ppph(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pp_hh(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pp_hh(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    unsigned long int id; int index;
    int l; int m; int k;
    int a; int b; int c;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        a = m_intClass->blockArrays_pppm_ppph(1,i1);
        b = m_intClass->blockArrays_pppm_ppph(2,i1);
        c = m_intClass->blockArrays_pppm_ppph(3,i1);
        k = m_intClass->blockArrays_pppm_ppph(4,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            l = m_intClass->blockArrays_pp_hh(1,i2);
            m = m_intClass->blockArrays_pp_hh(2,i2);

            id = m_intClass->Identity_hhhppp(l,m,k,a,b,c);
            index = T3_elements_I.find(id)->second;
            returnMat(i1-range_lower1, i2-range_lower2) = T3_elements_A[index];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T2e_makemat(int channel1, int channel2){    //makes a 4x2 matrix

    int range_lower1 = m_intClass->indexHolder_ppmm_pphh(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppmm_pphh(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pm_hp(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pm_hp(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    unsigned long int id; int index;
    int l; int j; int k;
    int d; int b; int c;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        b = m_intClass->blockArrays_ppmm_pphh(1,i1);
        c = m_intClass->blockArrays_ppmm_pphh(2,i1);
        j = m_intClass->blockArrays_ppmm_pphh(3,i1);
        k = m_intClass->blockArrays_ppmm_pphh(4,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            l = m_intClass->blockArrays_pm_hp(1,i2);
            d = m_intClass->blockArrays_pm_hp(2,i2);

            //std::cout << l << " " << j << " " << k << " " << d << " " << b << " " << c << std::endl;
;
            id = m_intClass->Identity_hhhppp(l,j,k,d,b,c);
            index = T3_elements_I.find(id)->second;//T3_elements_I_um.find(id)->second;
            //std::cout << id << " " << T3_elements_I.find(id)->first <<std::endl;
            //std::cout << index << std::endl;
            returnMat(i1-range_lower1, i2-range_lower2) = T3_elements_A[index];
            //std::cout << "ARGH" << std::endl;
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T3b_makemat_1(int channel1, int channel2){    //makes a 2x2 matrix

    int range_lower1 = m_intClass->indexHolder_pm_hp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pm_hp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pm_ph(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pm_ph(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int i; int l;
    int a; int d;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        l = m_intClass->blockArrays_pm_hp(1,i1);
        d = m_intClass->blockArrays_pm_hp(2,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            a = m_intClass->blockArrays_pm_ph(1,i2);
            i = m_intClass->blockArrays_pm_ph(2,i2);

            id = m_intClass->Identity_hhpp(i,l,a,d);
            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T3b_makemat_2(int channel1, int channel2){    //makes a 3x1 matrix

    int range_lower1 = m_intClass->indexHolder_ppm_pph(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_pph(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_p(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_p(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    unsigned long int id;
    int e; int b;
    int a; int i;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        b = m_intClass->blockArrays_ppm_pph(1,i1);
        a = m_intClass->blockArrays_ppm_pph(2,i1);
        i = m_intClass->blockArrays_ppm_pph(3,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            e = m_intClass->blockArrays_p_p(1,i2);

            id = (unsigned long int)m_intClass->Identity_ppph(e,b,a,i);
            returnMat(i1-range_lower1, i2-range_lower2) = T3D_remap[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T3b_makemat_3(int channel1, int channel2){    //makes a 3x1 matrix

    int range_lower1 = m_intClass->indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_p(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_p(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int j; int k;
    int e; int c;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        j = m_intClass->blockArrays_ppm_hhp(1,i1);
        k = m_intClass->blockArrays_ppm_hhp(2,i1);
        c = m_intClass->blockArrays_ppm_hhp(3,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            e = m_intClass->blockArrays_p_p(1,i2);

            id = m_intClass->Identity_hhpp(j,k,e,c);
            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T3c_makemat_1(int channel1, int channel2){    //makes a 2x2 matrix

    int range_lower1 = m_intClass->indexHolder_pm_hp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pm_hp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pm_ph(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pm_ph(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int i; int l;
    int a; int d;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        i = m_intClass->blockArrays_pm_hp(1,i1);
        a = m_intClass->blockArrays_pm_hp(2,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            d = m_intClass->blockArrays_pm_ph(1,i2);
            l = m_intClass->blockArrays_pm_ph(2,i2);

            id = m_intClass->Identity_hhpp(i,l,a,d);
            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T3c_makemat_2(int channel1, int channel2){    //makes a 3x1 matrix

    int range_lower1 = m_intClass->indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_h(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_h(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    unsigned long int id;
    int i; int a;
    int m; int j;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        i = m_intClass->blockArrays_ppm_hhp(1,i1);
        j = m_intClass->blockArrays_ppm_hhp(2,i1);
        a = m_intClass->blockArrays_ppm_hhp(3,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            m = m_intClass->blockArrays_p_h(1,i2);

            id = (unsigned long int)m_intClass->Identity_hhhp(m,j,i,a);
            returnMat(i1-range_lower1, i2-range_lower2) = T3D_remap[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T3c_makemat_3(int channel1, int channel2){    //makes a 3x1 matrix

    int range_lower1 = m_intClass->indexHolder_ppm_pph(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_pph(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_h(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_h(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int m; int k;
    int b; int c;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        b = m_intClass->blockArrays_ppm_pph(1,i1);
        c = m_intClass->blockArrays_ppm_pph(2,i1);
        k = m_intClass->blockArrays_ppm_pph(3,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            m = m_intClass->blockArrays_p_h(1,i2);

            id = m_intClass->Identity_hhpp(m,k,b,c);
            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T3d_makemat_1(int channel1, int channel2){    //makes a 2x2 matrix

    int range_lower1 = m_intClass->indexHolder_pp_hh(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pp_hh(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pp_pp(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pp_pp(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int j; int k;
    int d; int e;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        j = m_intClass->blockArrays_pp_hh(1,i1);
        k = m_intClass->blockArrays_pp_hh(2,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            d = m_intClass->blockArrays_pp_pp(1,i2);
            e = m_intClass->blockArrays_pp_pp(2,i2);

            id = m_intClass->Identity_hhpp(j,k,d,e);
            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T3d_makemat_2(int channel1, int channel2){    //makes a 3x1 matrix

    int range_lower1 = m_intClass->indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_h(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_h(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    unsigned long int id;
    int j; int k;
    int c; int l;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        j = m_intClass->blockArrays_ppm_hhp(1,i1);
        k = m_intClass->blockArrays_ppm_hhp(2,i1);
        c = m_intClass->blockArrays_ppm_hhp(3,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            l = m_intClass->blockArrays_p_h(1,i2);

            id = (unsigned long int)m_intClass->Identity_hhhp(j,k,l,c);
            returnMat(i1-range_lower1, i2-range_lower2) = T3D_remap[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T3d_makemat_3(int channel1, int channel2){    //makes a 3x1 matrix

    int range_lower1 = m_intClass->indexHolder_ppm_pph(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_pph(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_h(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_h(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int i; int l;
    int a; int b;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        a = m_intClass->blockArrays_ppm_pph(1,i1);
        b = m_intClass->blockArrays_ppm_pph(2,i1);
        i = m_intClass->blockArrays_ppm_pph(3,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            l = m_intClass->blockArrays_p_h(1,i2);

            id = m_intClass->Identity_hhpp(i,l,a,b);
            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T3e_makemat_1(int channel1, int channel2){    //makes a 3x1 matrix

    int range_lower1 = m_intClass->indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_p(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_p(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int i; int j;
    int a; int d;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        i = m_intClass->blockArrays_ppm_hhp(1,i1);
        j = m_intClass->blockArrays_ppm_hhp(2,i1);
        a = m_intClass->blockArrays_ppm_hhp(3,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            d = m_intClass->blockArrays_p_p(1,i2);

            id = m_intClass->Identity_hhpp(i,j,a,d);
            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T3e_makemat_2(int channel1, int channel2){    //makes a 4x2 matrix

    int range_lower1 = m_intClass->indexHolder_pppm_hhhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pppm_hhhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pp_hh(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pp_hh(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    unsigned long int id;
    int i; int j; int a;
    int l; int m; int k;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        i = m_intClass->blockArrays_pppm_hhhp(1,i1);
        j = m_intClass->blockArrays_pppm_hhhp(2,i1);
        k = m_intClass->blockArrays_pppm_hhhp(3,i1);
        a = m_intClass->blockArrays_pppm_hhhp(4,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            l = m_intClass->blockArrays_pp_hh(1,i2);
            m = m_intClass->blockArrays_pp_hh(2,i2);

            id = m_intClass->Identity_hhhhhp(i,j,l,m,k,a);
            returnMat(i1-range_lower1, i2-range_lower2) = T3D_remap[id];
            //std::cout << T3D_remap[id] << std::endl;
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T3e_makemat_3(int channel1, int channel2){    //makes a 2x2 matrix

    int range_lower1 = m_intClass->indexHolder_pp_hh(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pp_hh(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pp_pp(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pp_pp(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int l; int m;
    int b; int c;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        l = m_intClass->blockArrays_pp_hh(1,i1);
        m = m_intClass->blockArrays_pp_hh(2,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            b = m_intClass->blockArrays_pp_pp(1,i2);
            c = m_intClass->blockArrays_pp_pp(2,i2);

            id = m_intClass->Identity_hhpp(l,m,b,c);
            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T5a_makemat_1(int channel1, int channel2){    //makes a 2x2 matrix

    int range_lower1 = m_intClass->indexHolder_pp_hp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pp_hp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pp_ph(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pp_ph(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int i; int l;
    int a; int d;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        l = m_intClass->blockArrays_pp_hp(1,i1);
        d = m_intClass->blockArrays_pp_hp(2,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            a = m_intClass->blockArrays_pp_ph(1,i2);
            i = m_intClass->blockArrays_pp_ph(2,i2);

            id = m_intClass->Identity_hhpp(i,l,a,d);
            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T5a_makemat_2(int channel1, int channel2){    //makes a 4x2 matrix

    int range_lower1 = m_intClass->indexHolder_ppmm_hhpp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppmm_hhpp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pm_ph(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pm_ph(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    unsigned long int id; int index;
    int m; int j; int k;
    int e; int b; int c;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        j = m_intClass->blockArrays_ppmm_hhpp(1,i1);
        k = m_intClass->blockArrays_ppmm_hhpp(2,i1);
        b = m_intClass->blockArrays_ppmm_hhpp(3,i1);
        c = m_intClass->blockArrays_ppmm_hhpp(4,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            e = m_intClass->blockArrays_pm_ph(1,i2);
            m = m_intClass->blockArrays_pm_ph(2,i2);

            id = m_intClass->Identity_hhhppp(m,j,k,e,b,c);
            index = T3_elements_I.find(id)->second;//T3_elements_IV[thread][id];
            returnMat(i1-range_lower1, i2-range_lower2) = T3_elements_A[index];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T5b_makemat_1(int channel1, int channel2){    //makes a 3x1 matrix

    int range_lower1 = m_intClass->indexHolder_ppm_pph(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_pph(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_h(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_h(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int l; int i;
    int d; int e;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        d = m_intClass->blockArrays_ppm_pph(1,i1);
        e = m_intClass->blockArrays_ppm_pph(2,i1);
        l = m_intClass->blockArrays_ppm_pph(3,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            i = m_intClass->blockArrays_p_h(1,i2);

            id = m_intClass->Identity_hhpp(l,i,d,e);
            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T5b_makemat_2(int channel1, int channel2){    //makes a 5x1 matrix

    int range_lower1 = m_intClass->indexHolder_pppmm_ppphh(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pppmm_ppphh(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_h(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_h(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    unsigned long int id; int index;
    int m; int j; int k;
    int a; int b; int c;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            //std::cout << "sup" << std::endl;

            a = m_intClass->blockArrays_pppmm_ppphh(1,i1);
            b = m_intClass->blockArrays_pppmm_ppphh(2,i1);
            c = m_intClass->blockArrays_pppmm_ppphh(3,i1);
            j = m_intClass->blockArrays_pppmm_ppphh(4,i1);
            k = m_intClass->blockArrays_pppmm_ppphh(5,i1);
            m = m_intClass->blockArrays_p_h(1,i2);

            id = m_intClass->Identity_hhhppp(m,j,k,a,b,c);
            index = T3_elements_I.find(id)->second;
            returnMat(i1-range_lower1, i2-range_lower2) = T3_elements_A[index];
        }
    }

    return returnMat;
}

Eigen::MatrixXi MakeAmpMat::T5b_makemat_2_I(int channel1, int channel2){    //makes a 5x1 matrix

    int range_lower1 = m_intClass->indexHolder_pppmm_ppphh(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pppmm_ppphh(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_h(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_h(1, channel2);

    Eigen::MatrixXi returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    unsigned long int id; int index;
    int m; int j; int k;
    int a; int b; int c;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            //std::cout << "sup" << std::endl;

            a = m_intClass->blockArrays_pppmm_ppphh(1,i1);
            b = m_intClass->blockArrays_pppmm_ppphh(2,i1);
            c = m_intClass->blockArrays_pppmm_ppphh(3,i1);
            j = m_intClass->blockArrays_pppmm_ppphh(4,i1);
            k = m_intClass->blockArrays_pppmm_ppphh(5,i1);
            m = m_intClass->blockArrays_p_h(1,i2);

            id = m_intClass->Identity_hhhppp(m,j,k,a,b,c);
            index = T3_elements_I.find(id)->second;
            returnMat(i1-range_lower1, i2-range_lower2) = index;
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T5c_makemat_1(int channel1, int channel2){    //makes a 3x1 matrix

    int range_lower1 = m_intClass->indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_p(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_p(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int l; int m;
    int d; int a;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        l = m_intClass->blockArrays_ppm_hhp(1,i1);
        m = m_intClass->blockArrays_ppm_hhp(2,i1);
        d = m_intClass->blockArrays_ppm_hhp(3,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            a = m_intClass->blockArrays_p_p(1,i2);

            id = m_intClass->Identity_hhpp(l,m,d,a);
            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T5c_makemat_2(int channel1, int channel2){    //makes a 5x1 matrix

    int range_lower1 = m_intClass->indexHolder_pppmm_hhhpp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pppmm_hhhpp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_p(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_p(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    unsigned long int id; int index;
    int i; int j; int k;
    int e; int b; int c;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            //std::cout << "sup" << std::endl;

            i = m_intClass->blockArrays_pppmm_hhhpp(1,i1);
            j = m_intClass->blockArrays_pppmm_hhhpp(2,i1);
            k = m_intClass->blockArrays_pppmm_hhhpp(3,i1);
            b = m_intClass->blockArrays_pppmm_hhhpp(4,i1);
            c = m_intClass->blockArrays_pppmm_hhhpp(5,i1);
            e = m_intClass->blockArrays_p_p(1,i2);

            id = m_intClass->Identity_hhhppp(i,j,k,e,b,c);
            index = T3_elements_I.find(id)->second;
            returnMat(i1-range_lower1, i2-range_lower2) = T3_elements_A[index];
        }
    }

    return returnMat;
}

Eigen::MatrixXi MakeAmpMat::T5c_makemat_2_I(int channel1, int channel2){    //makes a 5x1 matrix

    int range_lower1 = m_intClass->indexHolder_pppmm_hhhpp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pppmm_hhhpp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_p(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_p(1, channel2);

    Eigen::MatrixXi returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    unsigned long int id; int index;
    int i; int j; int k;
    int e; int b; int c;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            //std::cout << "sup" << std::endl;

            i = m_intClass->blockArrays_pppmm_hhhpp(1,i1);
            j = m_intClass->blockArrays_pppmm_hhhpp(2,i1);
            k = m_intClass->blockArrays_pppmm_hhhpp(3,i1);
            b = m_intClass->blockArrays_pppmm_hhhpp(4,i1);
            c = m_intClass->blockArrays_pppmm_hhhpp(5,i1);
            e = m_intClass->blockArrays_p_p(1,i2);

            id = m_intClass->Identity_hhhppp(i,j,k,e,b,c);
            index = T3_elements_I.find(id)->second;
            returnMat(i1-range_lower1, i2-range_lower2) = index;
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T5d_makemat_1(int channel1, int channel2){    //makes a 3x1 matrix

    int range_lower1 = m_intClass->indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_p(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_p(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int i; int j;
    int a; int d;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        i = m_intClass->blockArrays_ppm_hhp(1,i1);
        j = m_intClass->blockArrays_ppm_hhp(2,i1);
        a = m_intClass->blockArrays_ppm_hhp(3,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            d = m_intClass->blockArrays_p_p(1,i2);

            id = m_intClass->Identity_hhpp(i,j,a,d);
            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T5d_makemat_2(int channel1, int channel2){    //makes a 3x3 matrix

    int range_lower1 = m_intClass->indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_ppm_pph(0, channel2);
    int range_upper2 = m_intClass->indexHolder_ppm_pph(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    unsigned long int id; int index;
    int l; int m; int k;
    int b; int e; int c;

    for (int i2 = range_lower2; i2<range_upper2; i2++){
        b = m_intClass->blockArrays_ppm_pph(1,i2);
        c = m_intClass->blockArrays_ppm_pph(2,i2);
        k = m_intClass->blockArrays_ppm_pph(3,i2);
        for (int i1 = range_lower1; i1<range_upper1; i1++){
            l = m_intClass->blockArrays_ppm_hhp(1,i1);
            m = m_intClass->blockArrays_ppm_hhp(2,i1);
            e = m_intClass->blockArrays_ppm_hhp(3,i1);

            id = m_intClass->Identity_hhhppp(l,m,k,b,e,c);
            index = T3_elements_I.find(id)->second;//T3_elements_IV[thread][id];
            returnMat(i1-range_lower1, i2-range_lower2) = T3_elements_A[index];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T5e_makemat_1(int channel1, int channel2){    //makes a 3x1 matrix

    int range_lower1 = m_intClass->indexHolder_ppm_pph(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_pph(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_h(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_h(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int i; int l;
    int a; int b;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        a = m_intClass->blockArrays_ppm_pph(1,i1);
        b = m_intClass->blockArrays_ppm_pph(2,i1);
        i = m_intClass->blockArrays_ppm_pph(3,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            l = m_intClass->blockArrays_p_h(1,i2);

            id = m_intClass->Identity_hhpp(i,l,a,b);
            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T5e_makemat_2(int channel1, int channel2){    //makes a 3x3 matrix

    int range_lower1 = m_intClass->indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_ppm_pph(0, channel2);
    int range_upper2 = m_intClass->indexHolder_ppm_pph(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    unsigned long int id; int index;
    int j; int m; int k;
    int d; int e; int c;

    for (int i2 = range_lower2; i2<range_upper2; i2++){
        d = m_intClass->blockArrays_ppm_pph(1,i2);
        e = m_intClass->blockArrays_ppm_pph(2,i2);
        m = m_intClass->blockArrays_ppm_pph(3,i2);
        for (int i1 = range_lower1; i1<range_upper1; i1++){
            j = m_intClass->blockArrays_ppm_hhp(1,i1);
            k = m_intClass->blockArrays_ppm_hhp(2,i1);
            c = m_intClass->blockArrays_ppm_hhp(3,i1);

            id = m_intClass->Identity_hhhppp(j,m,k,d,e,c);
            index = T3_elements_I.find(id)->second;//T3_elements_IV[thread][id];
            returnMat(i1-range_lower1, i2-range_lower2) = T3_elements_A[index];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T5f_makemat_1(int channel1, int channel2){    //makes a 2x2 matrix

    int range_lower1 = m_intClass->indexHolder_pp_hh(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pp_hh(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pp_pp(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pp_pp(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int i; int j;
    int d; int e;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        i = m_intClass->blockArrays_pp_hh(1,i1);
        j = m_intClass->blockArrays_pp_hh(2,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            d = m_intClass->blockArrays_pp_pp(1,i2);
            e = m_intClass->blockArrays_pp_pp(2,i2);

            id = m_intClass->Identity_hhpp(i,j,d,e);
            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T5f_makemat_2(int channel1, int channel2){    //makes a 4x2 matrix

    int range_lower1 = m_intClass->indexHolder_pppm_ppph(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pppm_ppph(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pp_hh(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pp_hh(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    unsigned long int id; int index;
    int l; int m; int k;
    int a; int b; int c;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        a = m_intClass->blockArrays_pppm_ppph(1,i1);
        b = m_intClass->blockArrays_pppm_ppph(2,i1);
        c = m_intClass->blockArrays_pppm_ppph(3,i1);
        k = m_intClass->blockArrays_pppm_ppph(4,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            l = m_intClass->blockArrays_pp_hh(1,i2);
            m = m_intClass->blockArrays_pp_hh(2,i2);

            id = m_intClass->Identity_hhhppp(l,m,k,a,b,c);
            index = T3_elements_I.find(id)->second;//T3_elements_IV[thread][id];
            returnMat(i1-range_lower1, i2-range_lower2) = T3_elements_A[index];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T5g_makemat_1(int channel1, int channel2){    //makes a 2x2 matrix

    int range_lower1 = m_intClass->indexHolder_pp_hh(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pp_hh(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pp_pp(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pp_pp(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int l; int m;
    int a; int b;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        l = m_intClass->blockArrays_pp_hh(1,i1);
        m = m_intClass->blockArrays_pp_hh(2,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            a = m_intClass->blockArrays_pp_pp(1,i2);
            b = m_intClass->blockArrays_pp_pp(2,i2);

            id = m_intClass->Identity_hhpp(l,m,a,b);
            returnMat(i1-range_lower1, i2-range_lower2) = T2_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::T5g_makemat_2(int channel1, int channel2){    //makes a 4x2 matrix

    int range_lower1 = m_intClass->indexHolder_pppm_hhhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pppm_hhhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pp_pp(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pp_pp(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    unsigned long int id; int index;
    int i; int j; int k;
    int d; int e; int c;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        i = m_intClass->blockArrays_pppm_hhhp(1,i1);
        j = m_intClass->blockArrays_pppm_hhhp(2,i1);
        k = m_intClass->blockArrays_pppm_hhhp(3,i1);
        c = m_intClass->blockArrays_pppm_hhhp(4,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            d = m_intClass->blockArrays_pp_pp(1,i2);
            e = m_intClass->blockArrays_pp_pp(2,i2);

            id = m_intClass->Identity_hhhppp(i,j,k,d,e,c);
            index = T3_elements_I.find(id)->second;//T3_elements_IV[thread][id];
            returnMat(i1-range_lower1, i2-range_lower2) = T3_elements_A[index];
        }
    }

    return returnMat;
}

void MakeAmpMat::T3b_Inverse_temp(Eigen::MatrixXd &inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_pm_pp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pm_pp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pm_ph(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pm_ph(1, channel2);


    unsigned long int id;
    int e; int b;
    int a; int i;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        e = m_intClass->blockArrays_pm_pp(1,i1);
        b = m_intClass->blockArrays_pm_pp(2,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            a = m_intClass->blockArrays_pm_ph(1,i2);
            i = m_intClass->blockArrays_pm_ph(2,i2);

            id = (unsigned long int)m_intClass->Identity_ppph(e,b,a,i);
            T3D_remap[id] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}

void MakeAmpMat::T3c_Inverse_temp(Eigen::MatrixXd inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_pm_hh(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pm_hh(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pm_hp(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pm_hp(1, channel2);


    unsigned long int id;
    int i; int a;
    int m; int j;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        m = m_intClass->blockArrays_pm_hh(1,i1);
        j = m_intClass->blockArrays_pm_hh(2,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            i = m_intClass->blockArrays_pm_hp(1,i2);
            a = m_intClass->blockArrays_pm_hp(2,i2);

            id = (unsigned long int)m_intClass->Identity_hhhp(m,j,i,a);
            T3D_remap[id] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}

void MakeAmpMat::T3d_Inverse_temp(Eigen::MatrixXd &inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_pp_hh(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pp_hh(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pp_ph(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pp_ph(1, channel2);


    unsigned long int id;
    int j; int k;
    int c; int l;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        j = m_intClass->blockArrays_pp_hh(1,i1);
        k = m_intClass->blockArrays_pp_hh(2,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            c = m_intClass->blockArrays_pp_ph(1,i2);
            l = m_intClass->blockArrays_pp_ph(2,i2);

            id = (unsigned long int)m_intClass->Identity_hhhp(j,k,l,c);
            T3D_remap[id] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}

void MakeAmpMat::T3e_Inverse_temp(Eigen::MatrixXd &inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_ppm_hhh(0, channel2);
    int range_upper2 = m_intClass->indexHolder_ppm_hhh(1, channel2);


    unsigned long int id;
    int i; int j; int a;
    int l; int m; int k;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        i = m_intClass->blockArrays_ppm_hhp(1,i1);
        j = m_intClass->blockArrays_ppm_hhp(2,i1);
        a = m_intClass->blockArrays_ppm_hhp(3,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            l = m_intClass->blockArrays_ppm_hhh(1,i2);
            m = m_intClass->blockArrays_ppm_hhh(2,i2);
            k = m_intClass->blockArrays_ppm_hhh(3,i2);

            id = m_intClass->Identity_hhhhhp(i,j,l,m,k,a);
            T3D_remap[id] =  inMat(i1-range_lower1, i2-range_lower2);
            //std::cout << T3D_remap[id] << std::endl;
        }
    }
}

void MakeAmpMat::T1a_inverse(Eigen::MatrixXd &inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_ppm_pph(0, channel2);
    int range_upper2 = m_intClass->indexHolder_ppm_pph(1, channel2);


    unsigned long int id; int index;
    int i; int j; int k;
    int a; int b; int c;

    for (int i2 = range_lower2; i2<range_upper2; i2++){
        b = m_intClass->blockArrays_ppm_pph(1,i2);
        c = m_intClass->blockArrays_ppm_pph(2,i2);
        k = m_intClass->blockArrays_ppm_pph(3,i2);
        for (int i1 = range_lower1; i1<range_upper1; i1++){
            i = m_intClass->blockArrays_ppm_hhp(1,i1);
            j = m_intClass->blockArrays_ppm_hhp(2,i1);
            a = m_intClass->blockArrays_ppm_hhp(3,i1);

            id = m_intClass->Identity_hhhppp(i,j,k,a,b,c);

            /*id = m_intClass->Identity_hhhppp(m_intClass->blockArrays_ppm_hhp(1,i1),
                                             m_intClass->blockArrays_ppm_hhp(2,i1),
                                             m_intClass->blockArrays_ppm_pph(3,i2),
                                             m_intClass->blockArrays_ppm_hhp(3,i1),
                                             m_intClass->blockArrays_ppm_pph(1,i2),
                                             m_intClass->blockArrays_ppm_pph(2,i2));*/

            index = T3_elements_I.find(id)->second;//T3_elements_I.find(id)->second;//T3_elements_IV[thread][id];
            T3_elements_A_temp[index] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}


void MakeAmpMat::T1b_inverse(Eigen::MatrixXd &inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_ppm_pph(0, channel2);
    int range_upper2 = m_intClass->indexHolder_ppm_pph(1, channel2);


    unsigned long int id; int index;
    int i; int j; int k;
    int a; int b; int c;

    for (int i2 = range_lower2; i2<range_upper2; i2++){
        a = m_intClass->blockArrays_ppm_pph(1,i2);
        b = m_intClass->blockArrays_ppm_pph(2,i2);
        i = m_intClass->blockArrays_ppm_pph(3,i2);
        for (int i1 = range_lower1; i1<range_upper1; i1++){
            j = m_intClass->blockArrays_ppm_hhp(1,i1);
            k = m_intClass->blockArrays_ppm_hhp(2,i1);
            c = m_intClass->blockArrays_ppm_hhp(3,i1);

            id = m_intClass->Identity_hhhppp(i,j,k,a,b,c);

            /*id = m_intClass->Identity_hhhppp(m_intClass->blockArrays_ppm_pph(3,i2),
                                             m_intClass->blockArrays_ppm_hhp(1,i1),
                                             m_intClass->blockArrays_ppm_hhp(2,i1),
                                             m_intClass->blockArrays_ppm_pph(1,i2),
                                             m_intClass->blockArrays_ppm_pph(2,i2),
                                             m_intClass->blockArrays_ppm_hhp(3,i1));*/

            index = T3_elements_I.find(id)->second;//T3_elements_IV[thread][id];
            T3_elements_A_temp[index] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}

void MakeAmpMat::T2c_inverse(Eigen::MatrixXd &inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_pppm_hhhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pppm_hhhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pp_pp(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pp_pp(1, channel2);

    unsigned long int id; int index;
    int i; int j; int k;
    int a; int b; int c;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        i = m_intClass->blockArrays_pppm_hhhp(1,i1);
        j = m_intClass->blockArrays_pppm_hhhp(2,i1);
        k = m_intClass->blockArrays_pppm_hhhp(3,i1);
        c = m_intClass->blockArrays_pppm_hhhp(4,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            a = m_intClass->blockArrays_pp_pp(1,i2);
            b = m_intClass->blockArrays_pp_pp(2,i2);

            id = m_intClass->Identity_hhhppp(i,j,k,a,b,c);
            index = T3_elements_I.find(id)->second;//T3_elements_I.find(id)->second;//T3_elements_IV[thread][id];
            T3_elements_A_temp[index] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}

void MakeAmpMat::T2d_inverse(Eigen::MatrixXd &inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_pppm_ppph(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pppm_ppph(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pp_hh(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pp_hh(1, channel2);

    unsigned long int id; int index;
    int i; int j; int k;
    int a; int b; int c;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        a = m_intClass->blockArrays_pppm_ppph(1,i1);
        b = m_intClass->blockArrays_pppm_ppph(2,i1);
        c = m_intClass->blockArrays_pppm_ppph(3,i1);
        k = m_intClass->blockArrays_pppm_ppph(4,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            i = m_intClass->blockArrays_pp_hh(1,i2);
            j = m_intClass->blockArrays_pp_hh(2,i2);

            id = m_intClass->Identity_hhhppp(i,j,k,a,b,c);
            index = T3_elements_I.find(id)->second;//T3_elements_I.find(id)->second;//T3_elements_IV[thread][id];
            T3_elements_A_temp[index] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}

void MakeAmpMat::T2e_inverse(Eigen::MatrixXd &inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_ppmm_pphh(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppmm_pphh(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pm_hp(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pm_hp(1, channel2);

    unsigned long int id; int index;
    int i; int j; int k;
    int a; int b; int c;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        b = m_intClass->blockArrays_ppmm_pphh(1,i1);
        c = m_intClass->blockArrays_ppmm_pphh(2,i1);
        j = m_intClass->blockArrays_ppmm_pphh(3,i1);
        k = m_intClass->blockArrays_ppmm_pphh(4,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            i = m_intClass->blockArrays_pm_hp(1,i2);
            a = m_intClass->blockArrays_pm_hp(2,i2);

            id = m_intClass->Identity_hhhppp(i,j,k,a,b,c);
            index = T3_elements_I.find(id)->second;//T3_elements_I.find(id)->second;//T3_elements_IV[thread][id];
            T3_elements_A_temp[index] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}

void MakeAmpMat::T3b_inverse(Eigen::MatrixXd &inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_ppm_pph(0, channel2);
    int range_upper2 = m_intClass->indexHolder_ppm_pph(1, channel2);

    unsigned long int id; int index;
    int i; int j; int k;
    int a; int b; int c;

    int thread;
    #pragma omp parallel for num_threads(4) private(i,j,k,a,b,c,id,index,thread) firstprivate(range_lower1, range_upper1, range_lower2, range_upper2)
    for (int i2 = range_lower2; i2<range_upper2; i2++){
        thread = omp_get_thread_num();
        b = m_intClass->blockArrays_ppm_pph(1,i2);
        a = m_intClass->blockArrays_ppm_pph(2,i2);
        i = m_intClass->blockArrays_ppm_pph(3,i2);
        for (int i1 = range_lower1; i1<range_upper1; i1++){
            j = m_intClass->blockArrays_ppm_hhp(1,i1);
            k = m_intClass->blockArrays_ppm_hhp(2,i1);
            c = m_intClass->blockArrays_ppm_hhp(3,i1);

            id = m_intClass->Identity_hhhppp(i,j,k,a,b,c);
            index = T3_elements_I.find(id)->second;
            T3_elements_A_temp[index] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}

void MakeAmpMat::T3c_inverse(Eigen::MatrixXd &inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_ppm_pph(0, channel2);
    int range_upper2 = m_intClass->indexHolder_ppm_pph(1, channel2);

    unsigned long int id; int index;
    int i; int j; int k;
    int a; int b; int c;

    int thread;
    #pragma omp parallel for num_threads(4) private(i,j,k,a,b,c,id,index,thread) firstprivate(range_lower1, range_upper1, range_lower2, range_upper2)
    for (int i2 = range_lower2; i2<range_upper2; i2++){
        thread = omp_get_thread_num();
        b = m_intClass->blockArrays_ppm_pph(1,i2);
        c = m_intClass->blockArrays_ppm_pph(2,i2);
        k = m_intClass->blockArrays_ppm_pph(3,i2);
    for (int i1 = range_lower1; i1<range_upper1; i1++){
            i = m_intClass->blockArrays_ppm_hhp(1,i1);
            j = m_intClass->blockArrays_ppm_hhp(2,i1);
            a = m_intClass->blockArrays_ppm_hhp(3,i1);

            id = m_intClass->Identity_hhhppp(i,j,k,a,b,c);
            index = T3_elements_I.find(id)->second;
            T3_elements_A_temp[index] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}

void MakeAmpMat::T3d_inverse(Eigen::MatrixXd &inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_ppm_pph(0, channel2);
    int range_upper2 = m_intClass->indexHolder_ppm_pph(1, channel2);

    unsigned long int id; int index;
    int i; int j; int k;
    int a; int b; int c;

    int thread;
    #pragma omp parallel for num_threads(4) private(i,j,k,a,b,c,id,index,thread) firstprivate(range_lower1, range_upper1, range_lower2, range_upper2)
    for (int i2 = range_lower2; i2<range_upper2; i2++){
        thread = omp_get_thread_num();
        a = m_intClass->blockArrays_ppm_pph(1,i2);
        b = m_intClass->blockArrays_ppm_pph(2,i2);
        i = m_intClass->blockArrays_ppm_pph(3,i2);
        for (int i1 = range_lower1; i1<range_upper1; i1++){
            j = m_intClass->blockArrays_ppm_hhp(1,i1);
            k = m_intClass->blockArrays_ppm_hhp(2,i1);
            c = m_intClass->blockArrays_ppm_hhp(3,i1);

            id = m_intClass->Identity_hhhppp(i,j,k,a,b,c);
            index = T3_elements_I.find(id)->second;
            T3_elements_A_temp[index] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}

void MakeAmpMat::T3e_inverse(Eigen::MatrixXd &inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_pppm_hhhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pppm_hhhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pp_pp(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pp_pp(1, channel2);

    unsigned long int id; int index;
    int i; int j; int k;
    int a; int b; int c;

    int thread;
    #pragma omp parallel for num_threads(4) private(i,j,k,a,b,c,id,index,thread) firstprivate(range_lower1, range_upper1, range_lower2, range_upper2)
    for (int i1 = range_lower1; i1<range_upper1; i1++){
        thread = omp_get_thread_num();
        i = m_intClass->blockArrays_pppm_hhhp(1,i1);
        j = m_intClass->blockArrays_pppm_hhhp(2,i1);
        k = m_intClass->blockArrays_pppm_hhhp(3,i1);
        a = m_intClass->blockArrays_pppm_hhhp(4,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            b = m_intClass->blockArrays_pp_pp(1,i2);
            c = m_intClass->blockArrays_pp_pp(2,i2);

            id = m_intClass->Identity_hhhppp(i,j,k,a,b,c);
            index = T3_elements_I.find(id)->second;
            T3_elements_A_temp[index] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}

void MakeAmpMat::T5a_inverse(Eigen::MatrixXd &inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_ppmm_hhpp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppmm_hhpp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pm_ph(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pm_ph(1, channel2);

    unsigned long int id; int index;
    int i; int j; int k;
    int a; int b; int c;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        j = m_intClass->blockArrays_ppmm_hhpp(1,i1);
        k = m_intClass->blockArrays_ppmm_hhpp(2,i1);
        b = m_intClass->blockArrays_ppmm_hhpp(3,i1);
        c = m_intClass->blockArrays_ppmm_hhpp(4,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            a = m_intClass->blockArrays_pm_ph(1,i2);
            i = m_intClass->blockArrays_pm_ph(2,i2);

            id = m_intClass->Identity_hhhppp(i,j,k,a,b,c);
            index = T3_elements_I.find(id)->second;//T3_elements_IV[thread][id];
            T3_elements_A_temp[index] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}

void MakeAmpMat::T5b_inverse(Eigen::MatrixXd &inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_pppmm_ppphh(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pppmm_ppphh(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_h(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_h(1, channel2);

    unsigned long int id; int index;
    int i; int j; int k;
    int a; int b; int c;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            a = m_intClass->blockArrays_pppmm_ppphh(1,i1);
            b = m_intClass->blockArrays_pppmm_ppphh(2,i1);
            c = m_intClass->blockArrays_pppmm_ppphh(3,i1);
            j = m_intClass->blockArrays_pppmm_ppphh(4,i1);
            k = m_intClass->blockArrays_pppmm_ppphh(5,i1);
            i = m_intClass->blockArrays_p_h(1,i2);

            id = m_intClass->Identity_hhhppp(i,j,k,a,b,c);
            index = T3_elements_I.find(id)->second;
            T3_elements_A_temp[index] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}

void MakeAmpMat::T5c_inverse(Eigen::MatrixXd &inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_pppmm_hhhpp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pppmm_hhhpp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_p_p(0, channel2);
    int range_upper2 = m_intClass->indexHolder_p_p(1, channel2);

    unsigned long int id; int index;
    int i; int j; int k;
    int a; int b; int c;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            i = m_intClass->blockArrays_pppmm_hhhpp(1,i1);
            j = m_intClass->blockArrays_pppmm_hhhpp(2,i1);
            k = m_intClass->blockArrays_pppmm_hhhpp(3,i1);
            b = m_intClass->blockArrays_pppmm_hhhpp(4,i1);
            c = m_intClass->blockArrays_pppmm_hhhpp(5,i1);
            a = m_intClass->blockArrays_p_p(1,i2);

            id = m_intClass->Identity_hhhppp(i,j,k,a,b,c);
            index = T3_elements_I.find(id)->second;
            T3_elements_A_temp[index] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}

void MakeAmpMat::T5b_inverse_I(Eigen::MatrixXd &inMat, int channel){

    int rows = T3_T5b_indices[channel].rows();
    int cols = T3_T5b_indices[channel].cols();
    int index;

    #pragma omp parallel for num_threads(4) private(index) firstprivate(rows, cols, channel)
    for (int i1 = 0; i1<rows; i1++){
        for (int i2 = 0; i2<cols; i2++){
            index = T3_T5b_indices[channel](i1,i2);
            T3_elements_A_temp[index] =  inMat(i1,i2);
        }
    }
}

void MakeAmpMat::T5c_inverse_I(Eigen::MatrixXd &inMat, int channel){

    int rows = T3_T5c_indices[channel].rows();
    int cols = T3_T5c_indices[channel].cols();
    int index;

    #pragma omp parallel for num_threads(4) private(index) firstprivate(rows, cols, channel)
    for (int i1 = 0; i1<rows; i1++){
        for (int i2 = 0; i2<cols; i2++){
            index = T3_T5c_indices[channel](i1,i2);
            T3_elements_A_temp[index] =  inMat(i1,i2);
        }
    }
}

void MakeAmpMat::T5d_inverse(Eigen::MatrixXd &inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_ppm_pph(0, channel2);
    int range_upper2 = m_intClass->indexHolder_ppm_pph(1, channel2);

    unsigned long int id; int index;
    int i; int j; int k;
    int a; int b; int c;

    for (int i2 = range_lower2; i2<range_upper2; i2++){
        b = m_intClass->blockArrays_ppm_pph(1,i2);
        c = m_intClass->blockArrays_ppm_pph(2,i2);
        k = m_intClass->blockArrays_ppm_pph(3,i2);
        for (int i1 = range_lower1; i1<range_upper1; i1++){

            i = m_intClass->blockArrays_ppm_hhp(1,i1);
            j = m_intClass->blockArrays_ppm_hhp(2,i1);
            a = m_intClass->blockArrays_ppm_hhp(3,i1);

            id = m_intClass->Identity_hhhppp(i,j,k,a,b,c);
            index = T3_elements_I.find(id)->second;//T3_elements_IV[thread][id];
            T3_elements_A_temp[index] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}

void MakeAmpMat::T5e_inverse(Eigen::MatrixXd &inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_ppm_pph(0, channel2);
    int range_upper2 = m_intClass->indexHolder_ppm_pph(1, channel2);

    unsigned long int id; int index;
    int i; int j; int k;
    int a; int b; int c;

    for (int i2 = range_lower2; i2<range_upper2; i2++){
        a = m_intClass->blockArrays_ppm_pph(1,i2);
        b = m_intClass->blockArrays_ppm_pph(2,i2);
        i = m_intClass->blockArrays_ppm_pph(3,i2);
        for (int i1 = range_lower1; i1<range_upper1; i1++){
            j = m_intClass->blockArrays_ppm_hhp(1,i1);
            k = m_intClass->blockArrays_ppm_hhp(2,i1);
            c = m_intClass->blockArrays_ppm_hhp(3,i1);

            id = m_intClass->Identity_hhhppp(i,j,k,a,b,c);
            index = T3_elements_I.find(id)->second;//T3_elements_IV[thread][id];
            T3_elements_A_temp[index] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}

void MakeAmpMat::T5f_inverse(Eigen::MatrixXd inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_pppm_ppph(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pppm_ppph(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pp_hh(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pp_hh(1, channel2);

    unsigned long int id; int index;
    int i; int j; int k;
    int a; int b; int c;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        a = m_intClass->blockArrays_pppm_ppph(1,i1);
        b = m_intClass->blockArrays_pppm_ppph(2,i1);
        c = m_intClass->blockArrays_pppm_ppph(3,i1);
        k = m_intClass->blockArrays_pppm_ppph(4,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            i = m_intClass->blockArrays_pp_hh(1,i2);
            j = m_intClass->blockArrays_pp_hh(2,i2);

            id = m_intClass->Identity_hhhppp(i,j,k,a,b,c);
            index = T3_elements_I.find(id)->second;//T3_elements_IV[thread][id];
            T3_elements_A_temp[index] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}

void MakeAmpMat::T5g_inverse(Eigen::MatrixXd &inMat, int channel1, int channel2){

    int range_lower1 = m_intClass->indexHolder_pppm_hhhp(0, channel1);
    int range_upper1 = m_intClass->indexHolder_pppm_hhhp(1, channel1);
    int range_lower2 = m_intClass->indexHolder_pp_pp(0, channel2);
    int range_upper2 = m_intClass->indexHolder_pp_pp(1, channel2);

    unsigned long int id; int index;
    int i; int j; int k;
    int a; int b; int c;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        i = m_intClass->blockArrays_pppm_hhhp(1,i1);
        j = m_intClass->blockArrays_pppm_hhhp(2,i1);
        k = m_intClass->blockArrays_pppm_hhhp(3,i1);
        c = m_intClass->blockArrays_pppm_hhhp(4,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            a = m_intClass->blockArrays_pp_pp(1,i2);
            b = m_intClass->blockArrays_pp_pp(2,i2);

            id = m_intClass->Identity_hhhppp(i,j,k,a,b,c);
            index = T3_elements_I.find(id)->second;//T3_elements_IV[thread][id];
            T3_elements_A_temp[index] =  inMat(i1-range_lower1, i2-range_lower2);
        }
    }
}

Eigen::MatrixXd MakeAmpMat::make3x3Block_I_D10c(int ku, int i1, int i2, int i3, int i4, int i5, int i6, std::vector<double> &T_vec){

    bool cond_hhh1 = (i1 == 0 && i2 == 0 && i3==0);
    bool cond_pph1 = (i1 == 1 && i2 == 1 && i3==0);
    bool cond_hhp1 = (i1 == 0 && i2 == 0 && i3==1);
    //bool cond_ppp1 = (i1 == 1 && i2 == 1 && i3==1);

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
    unsigned long int id; int index;
    if (cond_hhp1 && cond_pph2){
        for (int i = range_lower1; i<range_upper1; i++){
            for (int j = range_lower2; j<range_upper2; j++){
                id = m_intClass->Identity_hhhppp((blockArrays2_pointer)(3,j), (blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays1_pointer)(3,i), (blockArrays2_pointer)(1,j), (blockArrays2_pointer)(2,j));
                index = T3_elements_I.find(id)->second;
                returnMat(i-range_lower1, j-range_lower2) = T_vec[index];
            }
        }
    }
    else if (cond_hhh1 && cond_ppp2){
        for (int i = range_lower1; i<range_upper1; i++){
            for (int j = range_lower2; j<range_upper2; j++){
                id = m_intClass->Identity_hhhppp((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays1_pointer)(3,i), (blockArrays2_pointer)(1,j), (blockArrays2_pointer)(2,j), (blockArrays2_pointer)(3,j));
                index = T3_elements_I.find(id)->second;
                returnMat(i-range_lower1, j-range_lower2) = T_vec[index];
            }
        }
    }
    return returnMat;
}

//returns a block matrix of dimensions 3x1, currently only made for Vhhpp
// i1,i2,i3,i4 specify whether there is a hole or particle (by a 0 or 1)  index at index ij, for j=1-4
Eigen::MatrixXd MakeAmpMat::make3x3Block(int ku, int i1, int i2, int i3, int i4, int i5, int i6, spp::sparse_hash_map<int, double> &T_list){

    bool cond_hhh1 = (i1 == 0 && i2 == 0 && i3==0);
    bool cond_pph1 = (i1 == 1 && i2 == 1 && i3==0);
    bool cond_hhp1 = (i1 == 0 && i2 == 0 && i3==1);
    //bool cond_ppp1 = (i1 == 1 && i2 == 1 && i3==1);

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
    unsigned long int id;
    if (cond_hhp1 && cond_pph2){
        for (int i = range_lower1; i<range_upper1; i++){
            for (int j = range_lower2; j<range_upper2; j++){
                id = m_intClass->Identity_hhhppp((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays2_pointer)(3,j), (blockArrays2_pointer)(1,j), (blockArrays2_pointer)(2,j), (blockArrays1_pointer)(3,i));
                //std::cout << id << std::endl;
                returnMat(i-range_lower1, j-range_lower2) = T_list[id];

                /*for(auto it = T_list.cbegin(); it != T_list.cend(); ++it)
                {
                    std::cout << it->second << "\n";
                }*/
            }
        }
    }
    else if (cond_hhh1 && cond_ppp2){
        for (int i = range_lower1; i<range_upper1; i++){
            for (int j = range_lower2; j<range_upper2; j++){
                id = m_intClass->Identity_hhhppp((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays1_pointer)(3,i), (blockArrays2_pointer)(1,j), (blockArrays2_pointer)(2,j), (blockArrays2_pointer)(3,j));
                returnMat(i-range_lower1, j-range_lower2) = T_list[id];
            }
        }
        //std::cout << "sup2" << std::endl;
    }
    //std::cout << "sup2" << std::endl;
    return returnMat;
}

Eigen::MatrixXd MakeAmpMat::make3x1Block(int ku, int i1, int i2, int i3, int i4, spp::sparse_hash_map<int, double> &T_list){

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
    //int id;
    //std::cout << "sup1" << std::endl;
    if (cond_hhp && cond_p){
        for (int i = range_lower1; i<range_upper1; i++){
            for (int j = range_lower2; j<range_upper2; j++){
                int id = m_intClass->Identity_hhpp((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays1_pointer)(3,i), (blockArrays2_pointer)(1,j));
                returnMat(i-range_lower1, j-range_lower2) = T_list[id];
            }
        }
    }
    else if (cond_pph && cond_h){
        for (int i = range_lower1; i<range_upper1; i++){
            for (int j = range_lower2; j<range_upper2; j++){
                int id = m_intClass->Identity_hhpp((blockArrays1_pointer)(3,i), (blockArrays2_pointer)(1,j), (blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i));
                returnMat(i-range_lower1, j-range_lower2) = T_list[id];
            }
        }
    }
    //std::cout << "sup2" << std::endl;
    return returnMat;
}


//returns a block matrix of dimensions 2x2, currently only made for Vhhpp
// i1,i2,i3,i4 specify whether there is a hole or particle (by a 0 or 1) index at index ij, for j=1-4
Eigen::MatrixXd MakeAmpMat::make2x2Block(int ku, int i1, int i2, int i3, int i4, spp::sparse_hash_map<int, double> &T_list){

    //std::cout << "hey" << std::endl;
    bool cond_hh1 = (i1 == 0 && i2 == 0);
    bool cond_hp1 = (i1 == 0 && i2 == 1);
    bool cond_ph1 = (i1 == 1 && i2 == 0);
    //bool cond_pp1 = (i1 == 1 && i2 == 1);

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
