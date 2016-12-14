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

//0 means hole, 1 means particle
void  MakeIntMat::mapper(int i1, int i2){

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
                blockArrays_temp.col(index) << m_system->kUnique2(i,j),i,j;
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
                blockArrays_temp.col(index) << m_system->kUnique2(i,a),i,a;
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
                blockArrays_temp.col(index) << m_system->kUnique2(a,i),a,i;
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
                blockArrays_temp.col(index) << m_system->kUnique2(a,b),a,b;
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

    // 0 0
    if (cond_hh){
        blockArrays_hh = blockArrays;
        sortVec_hh     = sortVec;
    }
    // 0 1
    else if (cond_hp){
        blockArrays_hp = blockArrays;
        sortVec_hp     = sortVec;
    }
    // 1 0
    else if (cond_ph){
        blockArrays_ph = blockArrays;
        sortVec_ph     = sortVec;
    }
    // 1 1
    else{
        blockArrays_pp = blockArrays;
        sortVec_pp     = sortVec;
    }
}

void MakeIntMat::makeBlockMat(System* system, int Nh, int Ns){

    m_system = system;
    m_Nh = Nh;
    m_Ns = Ns;

    mapper(0,0); //hh
    cout << "made blockArrays and sortVec for hh" << endl;

    mapper(1,1); //pp
    cout << "made blockArrays and sortVec for pp" << endl;

    int counter         = 0;

    int range_lower     = 0;
    int range_upper     = 0;

    int range_lower_hh  = 0;
    int range_upper_hh  = 0;

    int range_lower_pp  = 0;
    int range_upper_pp  = 0;

    //DON'T DELETE
    //Vhhhh
    /*for (int l=0; l<sortVec_hh.size(); l++){
        int val = sortVec[l];
        for (int h=0; h<m_Nh*m_Nh; h++){ //goes through the entire blockArrays (which is sorted, so there'll be no interrupted repetitions)
            if ( val == blockArrays_hh(0,h) ){
                range_upper += 1;
            }
        }
        Vhhhh.push_back( makeSquareBlock(range_lower, range_upper) );//NEED MATRIX GENERATOR HERE
        range_lower = range_upper; //should i set range_upper + 1 ?
    }

    //Vhphp
    for (int l=0; l<sortVec_hh.size(); l++){
        int val = sortVec[l];
        for (int h=0; h<m_Nh*m_Nh; h++){ //goes through the entire blockArrays (which is sorted, so there'll be no interrupted repetitions)
            if ( val == blockArrays_hh(0,h) ){
                range_upper += 1;
            }
        }
        Vhhhh.push_back( makeSquareBlock(range_lower, range_upper) );//NEED MATRIX GENERATOR HERE
        range_lower = range_upper; //should i set range_upper + 1 ?
    }*/

    //Vhhpp
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
                Vhhpp.push_back( makeRektBlock(blockArrays_hh, blockArrays_pp,range_lower_hh, range_upper_hh, range_lower_pp, range_upper_pp) );
                Vhhpp_i.push_back( val_hh );
            }
        }
    }
    cout << "made Vhhpp" << endl;
    //Vpppp
    for (int l=0; l<sortVec_pp.size(); l++){
        int val = sortVec_pp[l];
        for (int h=0; h<(m_Ns-m_Nh)*(m_Ns-m_Nh); h++){ //goes through the entire blockArrays (which is sorted, so there'll be no interrupted repetitions)
            if ( val == blockArrays_pp(0,h) ){
                range_upper = h+1;
            }
        }
        Eigen::MatrixXf newElement = makeSquareBlock(blockArrays_pp, range_lower, range_upper);

        Vpppp.push_back( newElement );
        Vpppp_i.push_back( val );
        range_lower = range_upper; //doing this reset here is not a problem, unlike for that of Vhhpp
    }
    cout << "make Vpppp" << endl;
}



Eigen::MatrixXf MakeIntMat::makeSquareBlock(Eigen::MatrixXi& array, int range_lower, int range_upper){
    int dim = range_upper - range_lower;
    Eigen::MatrixXf returnMat;
    returnMat.conservativeResize(dim, dim);
    for (int i = range_lower; i<range_upper; i++){
        for (int j = range_lower; j<range_upper; j++){
            returnMat(i-range_lower,j-range_lower) = m_system->assym((array)(1,i), (array)(2,i), (array)(1,j), (array)(2,j));
        }
    }
    return returnMat;
}

Eigen::MatrixXf MakeIntMat::makeRektBlock(Eigen::MatrixXi& array1, Eigen::MatrixXi& array2, int range_lower1, int range_upper1, int range_lower2, int range_upper2){
    int dim1 = range_upper1 - range_lower1;
    int dim2 = range_upper2 - range_lower2;
    Eigen::MatrixXf returnMat;
    returnMat.conservativeResize(dim1, dim2);
    for (int i = range_lower1; i<range_upper1; i++){
        for (int j = range_lower2; j<range_upper2; j++){
            returnMat(i-range_lower1, j-range_lower2) = m_system->assym((array1)(1,i), (array1)(2,i), (array2)(1,j), (array2)(2,j));
        }
    }
    return returnMat;
}
