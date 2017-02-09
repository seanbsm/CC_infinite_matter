#include "makeintmat.h"
#include "Systems/system.h"
#include "Systems/heg.h"
#include <eigen3/Eigen/Dense>
#include <iostream>

using namespace std;

MakeIntMat::MakeIntMat()
{
}

template <typename T> //credit: http://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes#12399290
vector<size_t> sort_indexes(const vector<T> &v) {

    // initialize original index locations
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
         [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return idx;
}









// ##################################################
// ##                                              ##
// ## Mapper functions                             ##
// ##                                              ##
// ##################################################

//i only store this in case my rewrite fails
//0 means hole, 1 means particle
/*void MakeIntMat::mapper_1(std::vector<int>& sortVec, Eigen::MatrixXi& blockArrays, int i1, int s1){
    Eigen::MatrixXi   blockArrays_temp;
    int index       = 0;
    int colSize     = 0;
    int range_upper = 0;
    int range_lower = 0;

    bool cond_h = (i1 == 0);
    bool cond_p = (i1 == 1);

    // 0
    if (cond_h){
        colSize = m_Nh;
        blockArrays_temp.conservativeResize(2, colSize);

        for (int i=0; i<m_Nh; i++){
            blockArrays_temp.col(index) << m_system->kUnique1(i, s1),i;
            index += 1;

        }
    }
    // 1
    else{
        colSize = (m_Ns-m_Nh);
        blockArrays_temp.conservativeResize(2, colSize);
        for (int a=m_Nh; a<m_Ns; a++){
            blockArrays_temp.col(index) << m_system->kUnique1(a, s1),a;
            index += 1;

        }
    }

    Eigen::MatrixXi blockArrays = blockArrays_temp;     //need this
    Eigen::VectorXi veryTempVec = blockArrays_temp.row(0);
    std::vector<int> sortVec(veryTempVec.data(), veryTempVec.row(0).data() + blockArrays_temp.cols());
    std::vector<int> tempVec = sortVec;

    index = 0;
    for (auto i: sort_indexes(tempVec)){
        blockArrays.col(index) = blockArrays_temp.col(i);
        sortVec[index]         = tempVec[i];
        index += 1;
    }

    sortVec.erase( unique( sortVec.begin(), sortVec.end() ), sortVec.end() ); //need this

    // 0
    if (cond_h && s1==1){
        blockArrays_h = blockArrays;
        sortVec_h     = sortVec;
    }
    // 1
    else if (cond_p && s1==1){
        blockArrays_p = blockArrays;
        sortVec_p     = sortVec;
    }
}*/

//0 means hole, 1 means particle
void MakeIntMat::mapper_1(std::vector<int>& sortVecIn, Eigen::MatrixXi& blockArraysIn, int i1, int s1){
    Eigen::MatrixXi   blockArrays_temp;
    int index       = 0;
    int colSize     = 0;
    int range_upper = 0;
    int range_lower = 0;

    bool cond_h = (i1 == 0);
    bool cond_p = (i1 == 1);

    // 0
    if (cond_h){
        colSize = m_Nh;
        blockArrays_temp.conservativeResize(2, colSize);

        for (int i=0; i<m_Nh; i++){
            blockArrays_temp.col(index) << m_system->kUnique1(i, s1),i;
            index += 1;

        }
    }
    // 1
    else{
        colSize = (m_Ns-m_Nh);
        blockArrays_temp.conservativeResize(2, colSize);
        for (int a=m_Nh; a<m_Ns; a++){
            blockArrays_temp.col(index) << m_system->kUnique1(a, s1),a;
            index += 1;

        }
    }

    Eigen::MatrixXi blockArrays = blockArrays_temp;     //need this
    Eigen::VectorXi veryTempVec = blockArrays_temp.row(0);
    std::vector<int> sortVec(veryTempVec.data(), veryTempVec.row(0).data() + blockArrays_temp.cols());
    std::vector<int> tempVec = sortVec;

    index = 0;
    for (auto i: sort_indexes(tempVec)){
        blockArrays.col(index) = blockArrays_temp.col(i);
        sortVec[index]         = tempVec[i];
        index += 1;
    }

    sortVec.erase( unique( sortVec.begin(), sortVec.end() ), sortVec.end() ); //need this

    blockArraysIn = blockArrays;
    sortVecIn     = sortVec;
}


void  MakeIntMat::mapper_hp(){

    int index       = 0;
    int colSize     = m_Nh*(m_Ns-m_Nh);
    int range_upper = 0;
    int range_lower = 0;
    Eigen::MatrixXi   blockArrays_temp(3, colSize);

    for (int i=0; i<m_Nh; i++){
        for (int a=m_Nh; a<m_Ns; a++){
            blockArrays_temp.col(index) << m_system->kUnique2(i,a,1,1),i,a;
            index += 1;
        }
    }

    Eigen::MatrixXi blockArrays = blockArrays_temp;     //need this
    Eigen::VectorXi veryTempVec = blockArrays_temp.row(0);
    std::vector<int> sortVec(veryTempVec.data(), veryTempVec.row(0).data() + blockArrays_temp.cols());
    std::vector<int> tempVec = sortVec;

    index = 0;
    for (auto i: sort_indexes(tempVec)){
        blockArrays.col(index) = blockArrays_temp.col(i);
        sortVec[index]         = tempVec[i];
        index += 1;
    }

    sortVec.erase( unique( sortVec.begin(), sortVec.end() ), sortVec.end() ); //need this

    blockArrays_hp_s = blockArrays;
    sortVec_hp_s     = sortVec;
}



void  MakeIntMat::mapper_2(std::vector<int>& sortVecIn, Eigen::MatrixXi& blockArraysIn, int i1, int i2, int s1, int s2){

    Eigen::MatrixXi   blockArrays_temp;
    int index       = 0;
    int colSize     = 0;
    int range_upper = 0;
    int range_lower = 0;

    bool cond_hh = (i1 == 0 && i2 == 0);
    bool cond_hp = (i1 == 0 && i2 == 1);
    bool cond_ph = (i1 == 1 && i2 == 0);
    bool cond_pp = (i1 == 1 && i2 == 1);

    // 0 0
    if (cond_hh){
        colSize = m_Nh*m_Nh;
        blockArrays_temp.conservativeResize(3, colSize);

        for (int i=0; i<m_Nh; i++){
            for (int j=0; j<m_Nh; j++){
                blockArrays_temp.col(index) << m_system->kUnique2(i,j,s1,s2),i,j;
                index += 1;
            }
        }
    }
    // 0 1
    else if (cond_hp){
        colSize = (m_Ns-m_Nh)*m_Nh;
        blockArrays_temp.conservativeResize(3, colSize);

        for (int i=0; i<m_Nh; i++){
            for (int a=m_Nh; a<m_Ns; a++){
                blockArrays_temp.col(index) << m_system->kUnique2(i,a,s1,s2),i,a;
                index += 1;
            }
        }
    }
    // 1 0
    else if (cond_ph){
        colSize = (m_Ns-m_Nh)*m_Nh;
        blockArrays_temp.conservativeResize(3, colSize);

        for (int a=m_Nh; a<m_Ns; a++){
            for (int i=0; i<m_Nh; i++){
                blockArrays_temp.col(index) << m_system->kUnique2(a,i,s1,s2),a,i;
                index += 1;
            }
        }
    }
    // 1 1
    else{
        colSize = (m_Ns-m_Nh)*(m_Ns-m_Nh);
        blockArrays_temp.conservativeResize(3, colSize);
        for (int a=m_Nh; a<m_Ns; a++){
            for (int b=m_Nh; b<m_Ns; b++){
                blockArrays_temp.col(index) << m_system->kUnique2(a,b,s1,s2),a,b;
                index += 1;
            }
        }
    }

    Eigen::MatrixXi blockArrays = blockArrays_temp;     //need this
    Eigen::VectorXi veryTempVec = blockArrays_temp.row(0);
    std::vector<int> sortVec(veryTempVec.data(), veryTempVec.row(0).data() + blockArrays_temp.cols());
    std::vector<int> tempVec = sortVec;

    index = 0;
    for (auto i: sort_indexes(tempVec)){
        blockArrays.col(index) = blockArrays_temp.col(i);
        sortVec[index]         = tempVec[i];
        index += 1;
    }

    sortVec.erase( unique( sortVec.begin(), sortVec.end() ), sortVec.end() ); //need this

    blockArraysIn = blockArrays;
    sortVecIn     = sortVec;
}

void MakeIntMat::mapper_3(std::vector<int>& sortVecIn, Eigen::MatrixXi& blockArraysIn, int i1, int i2, int i3, int s1, int s2, int s3){

    Eigen::MatrixXi   blockArrays_temp;
    int index       = 0;
    int colSize     = 0;
    int range_upper = 0;
    int range_lower = 0;

    bool cond_hhh = (i1 == 0 && i2 == 0 && i3==0);
    bool cond_hhp = (i1 == 0 && i2 == 0 && i3==1);
    bool cond_pph = (i1 == 1 && i2 == 1 && i3==0);
    bool cond_ppp = (i1 == 1 && i2 == 1 && i3==1);

    // 0 0 0
    if (cond_hhh){
        colSize = (m_Ns-m_Nh)*m_Nh*m_Nh;
        blockArrays_temp.conservativeResize(4, colSize);

        for (int i=0; i<m_Nh; i++){
            for (int j=0; j<m_Nh; j++){
                for (int k=0; k<m_Nh; k++){
                    blockArrays_temp.col(index) << m_system->kUnique3(i,j,k,s1,s2,s3),i,j,k;
                    index += 1;
                }
            }
        }
    }
    // 0 0 1
    if (cond_hhp){
        colSize = (m_Ns-m_Nh)*m_Nh*m_Nh;
        blockArrays_temp.conservativeResize(4, colSize);

        for (int i=0; i<m_Nh; i++){
            for (int j=0; j<m_Nh; j++){
                for (int a=m_Nh; a<m_Ns; a++){
                    blockArrays_temp.col(index) << m_system->kUnique3(i,j,a,s1,s2,s3),i,j,a;
                    index += 1;
                }
            }
        }
    }
    // 1 1 0
    else if (cond_pph){
        colSize = (m_Ns-m_Nh)*(m_Ns-m_Nh)*m_Nh;
        blockArrays_temp.conservativeResize(4, colSize);

        for (int a=m_Nh; a<m_Ns; a++){
            for (int b=m_Nh; b<m_Ns; b++){
                for (int i=0; i<m_Nh; i++){
                    blockArrays_temp.col(index) << m_system->kUnique3(a,b,i,s1,s2,s3),a,b,i;
                    index += 1;
                }
            }
        }
    }
    // 1 1 1
    else if (cond_ppp){
        colSize = (m_Ns-m_Nh)*(m_Ns-m_Nh)*(m_Ns-m_Nh);
        blockArrays_temp.conservativeResize(4, colSize);

        for (int a=m_Nh; a<m_Ns; a++){
            for (int b=m_Nh; b<m_Ns; b++){
                for (int c=m_Nh; c<m_Ns; c++){
                    blockArrays_temp.col(index) << m_system->kUnique3(a,b,c,s1,s2,s3),a,b,c;
                    index += 1;
                }
            }
        }
    }


    Eigen::MatrixXi blockArrays = blockArrays_temp;     //need this
    Eigen::VectorXi veryTempVec = blockArrays_temp.row(0);
    std::vector<int> sortVec(veryTempVec.data(), veryTempVec.row(0).data() + blockArrays_temp.cols());
    std::vector<int> tempVec = sortVec;

    index = 0;
    for (auto i: sort_indexes(tempVec)){
        blockArrays.col(index) = blockArrays_temp.col(i);
        sortVec[index]         = tempVec[i];
        index += 1;
    }

    sortVec.erase( unique( sortVec.begin(), sortVec.end() ), sortVec.end() ); //need this

    blockArraysIn = blockArrays;
    sortVecIn     = sortVec;
}


// ##################################################
// ##                                              ##
// ## Make blocks                                  ##
// ##                                              ##
// ##################################################

//returns a block matrix of dimensions 3x1, currently only made for Vhhpp
// i1,i2,i3,i4 specify whether there is a hole or particle (by a 0 or 1)  index at index ij, for j=1-4
Eigen::MatrixXd MakeIntMat::make3x1Block(int ku, int i1, int i2, int i3, int i4){

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
        blockArrays1_pointer = blockArrays_hhp;
        sortVec1             = sortVec_hhp;
        indexHolder1_pointer = indexHolder_hhp;
    }
    // 1 1 0
    else if (cond_pph){
        blockArrays1_pointer = blockArrays_pph;
        sortVec1             = sortVec_pph;
        indexHolder1_pointer = indexHolder_pph;
    }

    // 0
    if (cond_h){
        blockArrays2_pointer = blockArrays_h;
        sortVec2             = sortVec_h;
        indexHolder2_pointer = indexHolder_h;
    }
    // 1
    else if (cond_p){
        blockArrays2_pointer = blockArrays_p;
        sortVec2             = sortVec_p;
        indexHolder2_pointer = indexHolder_p;
    }

    int length1 = indexHolder1_pointer.cols();
    int length2 = indexHolder2_pointer.cols();

    int index1; int index2;

    Eigen::MatrixXd returnMat;

    auto it1 = std::find(sortVec1.begin(), sortVec1.end(), ku);
    if (it1 == sortVec1.end()){
        returnMat.conservativeResize(1,1);
        returnMat(0,0) = 0;
      return returnMat;
    }
    else{
      index1 = distance(sortVec1.begin(), it1);
    }

    auto it2 = std::find(sortVec2.begin(), sortVec2.end(), ku);
    if (it2 == sortVec2.end()){
        returnMat.conservativeResize(1,1);
        returnMat(0,0) = 0;
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
                id = Identity_hhpp((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays1_pointer)(3,i), (blockArrays2_pointer)(1,j));
                returnMat(i-range_lower1, j-range_lower2) = Vhhpp_elements[id];
            }
        }
    }
    else if (cond_pph && cond_h){
        for (int i = range_lower1; i<range_upper1; i++){
            for (int j = range_lower2; j<range_upper2; j++){
                id = Identity_hhpp((blockArrays1_pointer)(3,i), (blockArrays2_pointer)(1,j), (blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i));
                returnMat(i-range_lower1, j-range_lower2) = Vhhpp_elements[id];
            }
        }
    }
    return returnMat;
}

Eigen::MatrixXd MakeIntMat::make2x2Block_alt(int channel){
    int dim1 = V_hh_pp[channel].rows();
    int dim2 = V_hh_pp[channel].cols();

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(dim1, dim2);
    for (int i = 0; i<dim1; i++){
        for (int j = 0; j<dim2; j++){
            //returnMat(i-range_lower1, j-range_lower2) = m_system->assym((array1)(1,i), (array1)(2,i), (array2)(1,j), (array2)(2,j));
            returnMat(i, j) = Vhhpp_vector[V_hh_pp[channel](i, j)];
        }
    }
    return returnMat;

}

//returns a block matrix of dimensions 2x2, currently only made for Vhhpp
// i1,i2,i3,i4 specify whether there is a hole or particle (by a 0 or 1) index at index ij, for j=1-4
Eigen::MatrixXd MakeIntMat::make2x2Block(int ku, int i1, int i2, int i3, int i4){

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
        blockArrays1_pointer = blockArrays_hh;
        sortVec1             = sortVec_hh;
        indexHolder1_pointer = indexHolder_hh;
    }
    // 0 1
    else if (cond_hp1){
        blockArrays1_pointer = blockArrays_hp;
        sortVec1             = sortVec_hp;
        indexHolder1_pointer = indexHolder_hp;
    }
    // 1 0
    else if (cond_ph1){
        blockArrays1_pointer = blockArrays_ph;
        sortVec1             = sortVec_ph;
        indexHolder1_pointer = indexHolder_ph;
    }
    // 1 1
    else {
        blockArrays1_pointer = blockArrays_pp;
        sortVec1             = sortVec_pp;
        indexHolder1_pointer = indexHolder_pp;
    }

    // 0 0
    if (cond_hh2){
        blockArrays2_pointer = blockArrays_hh;
        sortVec2             = sortVec_hh;
        indexHolder2_pointer = indexHolder_hh;
    }
    // 0 1
    else if (cond_hp2){
        blockArrays2_pointer = blockArrays_hp;
        sortVec2             = sortVec_hp;
        indexHolder2_pointer = indexHolder_hp;
    }
    // 1 0
    else if (cond_ph2){
        blockArrays2_pointer = blockArrays_ph;
        sortVec2             = sortVec_ph;
        indexHolder2_pointer = indexHolder_ph;
    }
    // 1 1
    else {
        blockArrays2_pointer = blockArrays_pp;
        sortVec2             = sortVec_pp;
        indexHolder2_pointer = indexHolder_pp;
    }

    int length1 = indexHolder1_pointer.cols();
    int length2 = indexHolder2_pointer.cols();

    int index1; int index2;

    Eigen::MatrixXd returnMat;

    auto it1 = std::find(sortVec1.begin(), sortVec1.end(), ku);
    if (it1 == sortVec1.end()){
        returnMat.conservativeResize(1,1);
        returnMat(0,0) = 0;
        std::cout << "make2x2Block in MakeIntMat, kUnique not found for rows" << std::endl;
      return returnMat;
    }
    else{
      index1 = distance(sortVec1.begin(), it1);
    }

    auto it2 = std::find(sortVec2.begin(), sortVec2.end(), ku);
    if (it2 == sortVec2.end()){
        returnMat.conservativeResize(1,1);
        returnMat(0,0) = 0;
        std::cout << "make2x2Block in MakeIntMat, kUnique not found for columns" << std::endl;
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
        int count1 = 0;
        int count2 = 0;
        for (int i = range_lower1; i<range_upper1; i++){
            for (int j = range_lower2; j<range_upper2; j++){
                count1 ++;
                returnMat(i-range_lower1, j-range_lower2) = Vhhpp_elements[Identity_hhpp((blockArrays1_pointer)(1,i), (blockArrays2_pointer)(2,j), (blockArrays1_pointer)(2,i), (blockArrays2_pointer)(1,j))];
                /*int ii = (blockArrays1_pointer)(1,i);
                int jj = (blockArrays2_pointer)(2,j);
                int aa = (blockArrays1_pointer)(2,i);
                int bb = (blockArrays2_pointer)(1,j);
                if (m_system->kUnique2(ii,jj,1,1) == m_system->kUnique2(aa,bb,1,1) ){
                    count2 ++;
                }*/
            }
        }
        /*if (count1 != count2){
            std::cout << count1 << " " << count2 << std::endl;
        }*/
    }
    if (cond_ph1 && cond_ph2){
        int count1 = 0;
        int count2 = 0;
        for (int i = range_lower1; i<range_upper1; i++){
            for (int j = range_lower2; j<range_upper2; j++){
                returnMat(i-range_lower1, j-range_lower2) = Vhhpp_elements[Identity_hhpp((blockArrays1_pointer)(2,i), (blockArrays2_pointer)(2,j), (blockArrays1_pointer)(1,i), (blockArrays2_pointer)(1,j))];
                /*int ii = (blockArrays1_pointer)(2,i);
                int jj = (blockArrays2_pointer)(2,j);
                int aa = (blockArrays1_pointer)(1,i);
                int bb = (blockArrays2_pointer)(1,j);
                if (m_system->kUnique2(ii,jj,1,1) == m_system->kUnique2(aa,bb,1,1) ){
                    count2 ++;
                }*/
            }
        }
        /*if (count1 != count2){
            std::cout << count1 << " " << count2 << std::endl;
        }*/
    }
    else if (cond_hh1 && cond_pp2){
        for (int i = range_lower1; i<range_upper1; i++){
            for (int j = range_lower2; j<range_upper2; j++){
                returnMat(i-range_lower1, j-range_lower2) = Vhhpp_elements[Identity_hhpp((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays2_pointer)(1,j), (blockArrays2_pointer)(2,j))];
            }
        }
    }
    return returnMat;
}









// ##################################################
// ##                                              ##
// ## Make interaction matrix elements             ##
// ##                                              ##
// ##################################################

void MakeIntMat::makeBlockMat(System* system, int Nh, int Ns){

    m_system = system;
    m_Nh = Nh;
    m_Ns = Ns;

    //The mapper functions take ascending values, i.e. they won't accept (1,0,1) as argument, but (0,1,1) works
    //I did this for several reasons, none of which I see the point in jotting down here

    mapper_1(sortVec_h, blockArrays_h, 0, 1);        //h
    mapper_1(sortVec_p, blockArrays_p, 1, 1);        //p

    mapper_2(sortVec_hh, blockArrays_hh, 0,0, 1,1);      //hh
    mapper_2(sortVec_hp, blockArrays_hp, 0,1, 1,-1);      //hp
    mapper_2(sortVec_ph, blockArrays_ph, 1,0, 1,-1);      //ph
    mapper_2(sortVec_pp, blockArrays_pp, 1,1, 1,1);      //pp
    mapper_hp();        //hp for Vhphp

    mapper_3(sortVec_hhp, blockArrays_hhp, 0,0,1, 1,1,-1);    //hhp
    mapper_3(sortVec_pph, blockArrays_pph, 1,1,0, 1,1,-1);    //pph

    cout << "made all blockArrays" << endl;

    int counter         = 0;

    int range_lower     = 0;
    int range_upper     = 0;

    int range_lower_hh  = 0;
    int range_upper_hh  = 0;

    int range_lower_hp  = 0;
    int range_upper_hp  = 0;

    int range_lower_pp  = 0;
    int range_upper_pp  = 0;

    //set indexHolders
    //I refrain from doing it for Vpppp for now, since I see no purpose as of yet for doing so.
    boundsHolder_hhpp_hh.conservativeResize(2, Eigen::NoChange);
    boundsHolder_hhpp_pp.conservativeResize(2, Eigen::NoChange);

    indexHolder_h.conservativeResize(2,sortVec_h.size());
    for (int h=0; h<sortVec_h.size(); h++){
        int val_h = sortVec_h[h];
        for (int bA_h=0; bA_h<m_Nh; bA_h++){
            if (val_h == blockArrays_h(0,bA_h)){
                range_upper = bA_h+1;
                counter += 1;
            }
        }
        range_lower = range_upper - counter;
        counter = 0;
        indexHolder_h.col(h) << range_lower, range_upper; //this now has same indices as sortVec
    }

    counter     = 0;
    range_lower = 0;
    range_upper = 0;
    indexHolder_p.conservativeResize(2,sortVec_p.size());
    for (int p=0; p<sortVec_p.size(); p++){
        int val_p = sortVec_p[p];
        for (int bA_p=0; bA_p<(m_Ns-m_Nh); bA_p++){
            if (val_p == blockArrays_p(0,bA_p)){
                range_upper = bA_p+1;
                counter += 1;
            }
        }
        range_lower = range_upper - counter;
        counter = 0;
        indexHolder_p.col(p) << range_lower, range_upper; //this now has same indices as sortVec
    }

    cout << "made one-particle indexHolders" << endl;

    counter     = 0;
    range_lower = 0;
    range_upper = 0;
    indexHolder_hh.conservativeResize(2,sortVec_hh.size());
    for (int hh=0; hh<sortVec_hh.size(); hh++){
        int val_hh = sortVec_hh[hh];
        for (int bA_hh=0; bA_hh<m_Nh*m_Nh; bA_hh++){
            if (val_hh == blockArrays_hh(0,bA_hh)){
                range_upper = bA_hh+1;
                counter += 1;
            }
        }
        range_lower = range_upper - counter;
        counter = 0;
        indexHolder_hh.col(hh) << range_lower, range_upper; //this now has same indices as sortVec
    }

    counter     = 0;
    range_lower = 0;
    range_upper = 0;
    indexHolder_hp.conservativeResize(2,sortVec_hp.size());
    for (int hp=0; hp<sortVec_hp.size(); hp++){
        int val_hp = sortVec_hp[hp];
        auto it = std::find(sortVec_ph.begin(), sortVec_ph.end(), val_hp);
        if (it != sortVec_ph.end()){
            for (int bA_hp=0; bA_hp<m_Nh*(m_Ns-m_Nh); bA_hp++){
                if (val_hp == blockArrays_hp(0,bA_hp)){
                    range_upper = bA_hp+1;
                    counter += 1;
                }
            }
        }
        range_lower = range_upper - counter;
        counter = 0;
        indexHolder_hp.col(hp) << range_lower, range_upper; //this now has same indices as sortVec
    }

    counter     = 0;
    range_lower = 0;
    range_upper = 0;
    indexHolder_ph.conservativeResize(2,sortVec_ph.size());
    for (int ph=0; ph<sortVec_ph.size(); ph++){
        int val_ph = sortVec_ph[ph];
        auto it = std::find(sortVec_hp.begin(), sortVec_hp.end(), val_ph);
        if (it != sortVec_hp.end()){
            for (int bA_ph=0; bA_ph<m_Nh*(m_Ns-m_Nh); bA_ph++){
                if (val_ph == blockArrays_ph(0,bA_ph)){
                    range_upper = bA_ph+1;
                    counter += 1;
                }
            }
        }
        range_lower = range_upper - counter;
        counter = 0;
        indexHolder_ph.col(ph) << range_lower, range_upper; //this now has same indices as sortVec
    }

    counter     = 0;
    range_lower = 0;
    range_upper = 0;
    indexHolder_hp_s.conservativeResize(2,sortVec_hp_s.size());
    for (int hp=0; hp<sortVec_hp_s.size(); hp++){
        int val_hp = sortVec_hp_s[hp];
        for (int bA_hp=0; bA_hp<m_Nh*(m_Ns-m_Nh); bA_hp++){
            if (val_hp == blockArrays_hp_s(0,bA_hp)){
                range_upper = bA_hp+1;
                counter += 1;
            }
        }
        range_lower = range_upper - counter;
        counter = 0;
        indexHolder_hp_s.col(hp) << range_lower, range_upper; //this now has same indices as sortVec
    }

    counter     = 0;
    range_lower = 0;
    range_upper = 0;
    indexHolder_pp.conservativeResize(2,sortVec_pp.size());
    for (int pp=0; pp<sortVec_pp.size(); pp++){
        int val_pp = sortVec_pp[pp];
        auto it = std::find(sortVec_hh.begin(), sortVec_hh.end(), val_pp);
        if (it != sortVec_hh.end()){
            for (int bA_pp=0; bA_pp<(m_Ns-m_Nh)*(m_Ns-m_Nh); bA_pp++){
                if (val_pp == blockArrays_pp(0,bA_pp)){
                    range_upper = bA_pp+1;
                    counter += 1;
                }
            }
        }
        range_lower = range_upper - counter;
        counter = 0;
        indexHolder_pp.col(pp) << range_lower, range_upper; //this now has same indices as sortVec
    }

    cout << "made two-particle indexHolders" << endl;

    counter     = 0;
    range_lower = 0;
    range_upper = 0;
    indexHolder_hhp.conservativeResize(2,sortVec_hhp.size());
    for (int hhp=0; hhp<sortVec_hhp.size(); hhp++){
        int val_hhp = sortVec_hhp[hhp];
        auto it = std::find(sortVec_p.begin(), sortVec_p.end(), val_hhp);
        if (it != sortVec_p.end()){
            for (int bA_hhp=0; bA_hhp<m_Nh*m_Nh*(m_Ns-m_Nh); bA_hhp++){
                if (val_hhp == blockArrays_hhp(0,bA_hhp)){
                    range_upper = bA_hhp+1;
                    counter += 1;
                }
            }
        }
        range_lower = range_upper - counter;
        counter = 0;
        indexHolder_hhp.col(hhp) << range_lower, range_upper; //this now has same indices as sortVec
    }

    counter     = 0;
    range_lower = 0;
    range_upper = 0;
    indexHolder_pph.conservativeResize(2,sortVec_pph.size());
    for (int pph=0; pph<sortVec_pph.size(); pph++){
        int val_pph = sortVec_pph[pph];
        auto it = std::find(sortVec_h.begin(), sortVec_h.end(), val_pph);
        if (it != sortVec_h.end()){
            for (int bA_pph=0; bA_pph<m_Nh*(m_Ns-m_Nh)*(m_Ns-m_Nh); bA_pph++){
                if (val_pph == blockArrays_pph(0,bA_pph)){
                    range_upper = bA_pph+1;
                    counter += 1;
                }
            }
        }
        range_lower = range_upper - counter;
        counter = 0;
        indexHolder_pph.col(pph) << range_lower, range_upper; //this now has same indices as sortVec
    }

    cout << "made three-particle indexHolders" << endl;

    counter     = 0;
    range_lower = 0;
    range_upper = 0;


    for (int h=0; h<sortVec_hh.size(); h++){
        for (int p=0; p<sortVec_pp.size(); p++){
            int val_hh = sortVec_hh[h];
            int val_pp = sortVec_pp[p];
            if (val_hh == val_pp){      //ensures I only work on cases where hh and pp have equal kunique
                for (int hh=0; hh<m_Nh*m_Nh; hh++){
                    if ( val_hh == blockArrays_hh(0,hh) ){
                        range_upper_hh = hh+1;
                        counter += 1;
                    }
                }
                range_lower_hh = range_upper_hh - counter;
                counter = 0;
                for (int pp=0; pp<(m_Ns-m_Nh)*(m_Ns-m_Nh); pp++){
                    if ( val_hh == blockArrays_pp(0,pp) ){
                        range_upper_pp = pp+1;
                        counter += 1;
                    }
                }
                range_lower_pp = range_upper_pp - counter;
                counter = 0;
                boundsHolder_hhpp_hh.conservativeResize(Eigen::NoChange, boundsHolder_hhpp_hh.cols()+1);
                boundsHolder_hhpp_pp.conservativeResize(Eigen::NoChange, boundsHolder_hhpp_pp.cols()+1);
                boundsHolder_hhpp_hh.col(boundsHolder_hhpp_hh.cols()-1) << range_lower_hh, range_upper_hh;
                boundsHolder_hhpp_pp.col(boundsHolder_hhpp_pp.cols()-1) << range_lower_pp, range_upper_pp;
                Vhhpp_i.push_back(val_hh);
            }
        }
    }

    numOfKu = boundsHolder_hhpp_hh.cols();
    //cout << numOfKu <<" "<< sortVec_hh.size() << endl;

    cout << "made indexHolders" << endl;

    /*for (int h=0; h<boundsHolder_hhpp_hh.cols(); h++){
        range_lower_hh = boundsHolder_hhpp_hh(0,h);
        range_upper_hh = boundsHolder_hhpp_hh(1,h);
        range_lower_pp = boundsHolder_hhpp_pp(0,h);
        range_upper_pp = boundsHolder_hhpp_pp(1,h);
        Vhhpp.push_back( makeRektBlock(blockArrays_hh, blockArrays_pp,range_lower_hh, range_upper_hh, range_lower_pp, range_upper_pp) );
        Vhhpp_i.push_back( sortVec_hh[h] );
    }*/

    /*for (int h=0; h<boundsHolder_hhpp_hh.cols(); h++){
        range_lower_hh = boundsHolder_hhpp_hh(0,h);
        range_upper_hh = boundsHolder_hhpp_hh(1,h);
        range_lower_pp = boundsHolder_hhpp_pp(0,h);
        range_upper_pp = boundsHolder_hhpp_pp(1,h);

        //makeMatMap_hhpp makes ALL non-zero elements of Vhhpp, there are no others
        //i.e. V_hhp_p can't have any elements outside those in Vhhpp
        makeMatMap_hhpp(blockArrays_hh, blockArrays_pp,range_lower_hh, range_upper_hh, range_lower_pp, range_upper_pp);

        //These two lines should no longer be needed
        Vhhpp.push_back( makeRektBlock(blockArrays_hh, blockArrays_pp,range_lower_hh, range_upper_hh, range_lower_pp, range_upper_pp) );
        //Vhhpp_i.push_back( sortVec_hh[h] );
    }*/


    //make Vhhpp map
    for (int h=0; h<boundsHolder_hhpp_hh.cols(); h++){
        range_lower_hh = boundsHolder_hhpp_hh(0,h);
        range_upper_hh = boundsHolder_hhpp_hh(1,h);
        range_lower_pp = boundsHolder_hhpp_pp(0,h);
        range_upper_pp = boundsHolder_hhpp_pp(1,h);
        makeMatMap_hhpp(blockArrays_hh, blockArrays_pp,range_lower_hh, range_upper_hh, range_lower_pp, range_upper_pp);
        makeMatVec(blockArrays_hh, blockArrays_pp,range_lower_hh, range_upper_hh, range_lower_pp, range_upper_pp);
    }

    //currently this is riddled with errors, like ph has convention a-i instead of a+i as it should be
    //and makeMatMap_hhpp uses Vhhpp, not Vppph
    /*for (int pp=0; pp<sortVec_pp.size(); pp++){
        for (int ph=0; ph<sortVec_ph.size(); ph++){
            if (sortVec_pp[pp] == sortVec_ph[ph] ){
                range_lower_pp = indexHolder_pp(0,pp);
                range_upper_pp = indexHolder_pp(1,pp);
                range_lower_ph = indexHolder_ph(0,ph);
                range_upper_ph = indexHolder_ph(1,ph);
                makeMatMap_hhpp(blockArrays_pp, blockArrays_ph,range_lower_pp, range_upper_pp, range_lower_ph, range_upper_ph);
            }
        }
    }*/


    cout << "made Vhhpp" << endl;










    //Below we make Vpppp, Vhhhh, and Vhphp, which are distinct from Vhhpp --> don't need any remapping of elements
    for (int h=0; h<boundsHolder_hhpp_hh.cols(); h++){
        range_lower_hh = boundsHolder_hhpp_hh(0,h);
        range_upper_hh = boundsHolder_hhpp_hh(1,h);
        Vhhhh.push_back( makeSquareBlock(blockArrays_hh, range_lower_hh, range_upper_hh) );
    }

    int prev = 0;
    int length = boundsHolder_hhpp_hh.cols();

    for (int hp=0; hp<indexHolder_hp.cols(); hp++){
        range_lower_hp = indexHolder_hp(0,hp);
        range_upper_hp = indexHolder_hp(1,hp);
        Vhphp.push_back( makeSquareBlock_s(blockArrays_hp, range_lower_hp, range_upper_hp) );
    }

    for (int p=0; p<boundsHolder_hhpp_hh.cols(); p++){
        range_lower_pp = boundsHolder_hhpp_pp(0,p);
        range_upper_pp = boundsHolder_hhpp_pp(1,p);
        //Vpppp.push_back( makeRektBlock(blockArrays_pp,blockArrays_pp, range_lower_pp, range_upper_pp,range_lower_pp, range_upper_pp) );
        Vpppp.push_back( makeSquareBlock(blockArrays_pp, range_lower_pp, range_upper_pp) );

        //this printout is handy for large systems
        /*if (prev < p/(double)length){
            cout << 100*p/(double)length << "%  " << range_upper_pp-range_lower_pp << endl;
            prev = p/(double)length;
        }*/
        //Vpppp_i.push_back( sortVec_pp[p] );
    }

    //for (int i = 0; i<Vhhpp_i.size(); i++){
    //    cout << blockArrays_pp(0,boundsHolder_hhpp_pp(0,i)) << endl;
    //}

    //Vpppp
    /*for (int l=0; l<sortVec_pp.size(); l++){
        int val = sortVec_pp[l];
        for (int h=0; h<(m_Ns-m_Nh)*(m_Ns-m_Nh); h++){ //goes through the entire blockArrays (which is sorted, so there'll be no interrupted repetitions)
            if ( val == blockArrays_pp(0,h) ){
                range_upper = h+1;
            }
        }
        Eigen::MatrixXd newElement = makeSquareBlock(blockArrays_pp, range_lower, range_upper);

        Vpppp.push_back( newElement );
        Vpppp_i.push_back( val );
        range_lower = range_upper; //doing this reset here is not a problem, unlike for that of Vhhpp
    }*/
    cout << "made Vpppp" << endl;
}










// ##################################################
// ##                                              ##
// ## Make channels for interaction matris         ##
// ##                                              ##
// ##################################################

//I get the same result using this as I get using Rektblock
Eigen::MatrixXd MakeIntMat::makeSquareBlock(Eigen::MatrixXi& array, int range_lower, int range_upper){
    int dim = range_upper - range_lower;
    //cout << dim << endl;
    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(dim, dim);
    for (int i = range_lower; i<range_upper; i++){
        for (int j = i; j<range_upper; j++){
            returnMat(i-range_lower,j-range_lower) = m_system->assym((array)(1,i), (array)(2,i), (array)(1,j), (array)(2,j));
        }
    }
    returnMat = returnMat.selfadjointView<Eigen::Upper>();  //Since Vpppp is symmetric, we reflect the upper half onto the lower half
    return returnMat;
}

Eigen::MatrixXd MakeIntMat::makeSquareBlock_s(Eigen::MatrixXi& array, int range_lower, int range_upper){
    int dim = range_upper - range_lower;
    //cout << dim << endl;
    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(dim, dim);
    for (int i = range_lower; i<range_upper; i++){
        for (int j = range_lower; j<range_upper; j++){
            returnMat(i-range_lower,j-range_lower) = m_system->assym((array)(1,j), (array)(2,i), (array)(1,i), (array)(2,j));
        }
    }
    //returnMat = returnMat.selfadjointView<Eigen::Upper>();  //Since Vpppp is symmetric, we reflect the upper half onto the lower half
    return returnMat;
}

Eigen::MatrixXd MakeIntMat::makeRektBlock(Eigen::MatrixXi& array1, Eigen::MatrixXi& array2, int range_lower1, int range_upper1, int range_lower2, int range_upper2){
    int dim1 = range_upper1 - range_lower1;
    int dim2 = range_upper2 - range_lower2;
    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(dim1, dim2);
    for (int i = range_lower1; i<range_upper1; i++){
        for (int j = range_lower2; j<range_upper2; j++){
            //returnMat(i-range_lower1, j-range_lower2) = m_system->assym((array1)(1,i), (array1)(2,i), (array2)(1,j), (array2)(2,j));
            returnMat(i-range_lower1, j-range_lower2) = Vhhpp_elements[Identity_hhpp((array1)(1,i), (array1)(2,i), (array2)(1,j), (array2)(2,j))];
        }
    }
    return returnMat;
}


// ##################################################
// ##                                              ##
// ## Make maps with interaction matrix elements   ##
// ##                                              ##
// ##################################################

void MakeIntMat::makeMatMap_hhhp(Eigen::MatrixXi& array1, Eigen::MatrixXi& array2, int range_lower1, int range_upper1, int range_lower2, int range_upper2){
    for (int i = range_lower1; i<range_upper1; i++){
        for (int j = range_lower2; j<range_upper2; j++){
            double val = m_system->assym((array1)(1,i), (array1)(2,i), (array2)(1,j), (array2)(2,j));
            if (val != 0){
                Vhhhp_elements[Identity_hhhp((array1)(1,i), (array1)(2,i), (array2)(1,j), (array2)(2,j))] = val;
            }

        }
    }
}

void MakeIntMat::makeMatMap_hhpp(Eigen::MatrixXi& array1, Eigen::MatrixXi& array2, int range_lower1, int range_upper1, int range_lower2, int range_upper2){
    for (int i = range_lower1; i<range_upper1; i++){
        for (int j = range_lower2; j<range_upper2; j++){
            double val = m_system->assym((array1)(1,i), (array1)(2,i), (array2)(1,j), (array2)(2,j));
            if (val != 0){
                Vhhpp_elements[Identity_hhpp((array1)(1,i), (array1)(2,i), (array2)(1,j), (array2)(2,j))] = val;
            }

        }
    }
}

void MakeIntMat::makeMatMap_hppp(Eigen::MatrixXi& array1, Eigen::MatrixXi& array2, int range_lower1, int range_upper1, int range_lower2, int range_upper2){
    for (int i = range_lower1; i<range_upper1; i++){
        for (int j = range_lower2; j<range_upper2; j++){
            double val = m_system->assym((array1)(1,i), (array1)(2,i), (array2)(1,j), (array2)(2,j));
            if (val != 0){
                Vhppp_elements[Identity_hppp((array1)(1,i), (array1)(2,i), (array2)(1,j), (array2)(2,j))] = val;
            }

        }
    }
}

void MakeIntMat::makeMatVec(Eigen::MatrixXi& array1, Eigen::MatrixXi& array2, int range_lower1, int range_upper1, int range_lower2, int range_upper2){
    //std::vector<double> tempVec;

    int dim1 = range_upper1 - range_lower1;
    int dim2 = range_upper2 - range_lower2;
    Eigen::MatrixXi insertMat;

    insertMat.conservativeResize(dim1, dim2);
    for (int i = range_lower1; i<range_upper1; i++){
        for (int j = range_lower2; j<range_upper2; j++){
            double val = m_system->assym((array1)(1,i), (array1)(2,i), (array2)(1,j), (array2)(2,j));
            if (val != 0){
                //tempVec.push_back(val);
                Vhhpp_vector.push_back(val);
                insertMat(i-range_lower1, j-range_lower2) = Vhhpp_counter;
                Vhhpp_counter += 1; //class counter, keeps value across function calls
            }
            else{
                insertMat(i-range_lower1, j-range_lower2) = -1;
            }
        }
    }

   // These lines are required (elsewhere than this function, I jsut haven't moved them) to have an Eigen vector
   // rather than an std vector
   // double* ptr = &tempVec[0];
   // Eigen::Map<Eigen::VectorXd> tempVec2(ptr, tempVec.size());
   // Vhhpp_vector = tempVec2;

    V_hh_pp.push_back(insertMat);
}

// ##################################################
// ##                                              ##
// ## Identity functions                           ##
// ##                                              ##
// ##################################################

int MakeIntMat::Identity_hhhp(int h1, int h2, int h3, int p1){
    return h1 + h2*m_Nh + h3*m_Nh*m_Nh + p1*m_Nh*m_Nh*(m_Ns-m_Nh);
}

int MakeIntMat::Identity_hhpp(int h1, int h2, int p1, int p2){
    return h1 + h2*m_Nh + p1*m_Nh*(m_Ns-m_Nh) + p2*m_Nh*(m_Ns-m_Nh)*(m_Ns-m_Nh);
}

int MakeIntMat::Identity_hppp(int h1, int p1, int p2, int p3){
    return h1 + p1*(m_Ns-m_Nh) + p2*(m_Ns-m_Nh)*(m_Ns-m_Nh) + p3*(m_Ns-m_Nh)*(m_Ns-m_Nh)*(m_Ns-m_Nh);
}
