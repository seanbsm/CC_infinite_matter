#include "makeintmat.h"
#include "Systems/system.h"
#include "Systems/heg.h"
#include "eigen3/Eigen/Dense"
#include <omp.h>
#include <mpi.h>
#include <iostream>
#include <iomanip>

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

void MakeIntMat::setThreads(int numthreads){
    m_numThreads = numthreads;
}

void MakeIntMat::makePermutations(){

    int Nh  = m_Nh;
    int Nh2 = Nh*Nh;
    int Np  = m_Ns;//m_Ns-m_Nh;
    int Np2 = Np*Np;

    int cols_h = blockArrays_ppp_hhh.cols();
    int cols_p = blockArrays_ppp_ppp.cols();

    //note some of these elements will be zero since they also hold channels not present in each other
    //however these will never be called, and thus not cause problems (in theory -.-' )
    blockArrays_ppp_hhh_Pij.conservativeResize(2,cols_h);
    blockArrays_ppp_hhh_Pik.conservativeResize(2,cols_h);
    blockArrays_ppp_hhh_Pjk.conservativeResize(2,cols_h);
    blockArrays_ppp_hhh_Pijik.conservativeResize(2,cols_h);
    blockArrays_ppp_hhh_Pijjk.conservativeResize(2,cols_h);

    blockArrays_ppp_ppp_Pab.conservativeResize(2,cols_p);
    blockArrays_ppp_ppp_Pac.conservativeResize(2,cols_p);
    blockArrays_ppp_ppp_Pbc.conservativeResize(2,cols_p);
    blockArrays_ppp_ppp_Pabac.conservativeResize(2,cols_p);
    blockArrays_ppp_ppp_Pabbc.conservativeResize(2,cols_p);

    int range_lower_h;
    int range_upper_h;
    int range_lower_p;
    int range_upper_p;

    int i1; int j1; int k1;
    int i2; int j2; int k2;
    int val_ij2; int val_ik2; int val_jk2; int val_ijik2; int val_ijjk2;

    int a1; int b1; int c1;
    int a2; int b2; int c2;
    int val_ab2; int val_ac2; int val_bc2; int val_abac2;; int val_abbc2;

    int val1;

    for (int channel = 0; channel<numOfKu3; channel++){
        range_lower_h = boundsHolder_hhhppp_hhh(0,channel);
        range_upper_h = boundsHolder_hhhppp_hhh(1,channel);
        range_lower_p = boundsHolder_hhhppp_ppp(0,channel);
        range_upper_p = boundsHolder_hhhppp_ppp(1,channel);

        //hhh permutations
        for (int it1 = range_lower_h; it1<range_upper_h; it1++){
            i1 = blockArrays_ppp_hhh(1,it1);
            j1 = blockArrays_ppp_hhh(2,it1);
            k1 = blockArrays_ppp_hhh(3,it1);

            val1 = i1 + j1*Nh + k1*Nh2;

            int counter = 0;

            /*if (i1 == j1 && i1 == k1){
                blockArrays_ppp_hhh_Pijk.col(it1) << it1, it1;
                counter ++;
            }
            if (i1 == j1){
                blockArrays_ppp_hhh_Pij.col(it1) << it1, it1;
                counter ++;
            }
            if (i1 == k1){
                blockArrays_ppp_hhh_Pik.col(it1) << it1, it1;
                counter ++;
            }
            if (j1 == k1){
                blockArrays_ppp_hhh_Pjk.col(it1) << it1, it1;
                counter ++;
            }*/

            //int n = omp_get_max_threads();
            #pragma omp parallel for num_threads(m_numThreads) private(i2,j2,k2, val_ij2,val_ik2,val_jk2,val_ijik2,val_ijjk2)
            for (int it2 = range_lower_h; it2<range_upper_h; it2++){
                i2 = blockArrays_ppp_hhh(1,it2);
                j2 = blockArrays_ppp_hhh(2,it2);
                k2 = blockArrays_ppp_hhh(3,it2);

                val_ij2  = j2 + i2*Nh + k2*Nh2;
                val_ik2  = k2 + j2*Nh + i2*Nh2;
                val_jk2  = i2 + k2*Nh + j2*Nh2;
                val_ijik2 = k2 + i2*Nh + j2*Nh2; // P(ki)P(kj): (i,j,k) -> (k,i,j)
                val_ijjk2 = j2 + k2*Nh + i2*Nh2; // P(ki)P(kj): (i,j,k) -> (j,k,i)

                if (val1 == val_ij2){
                    blockArrays_ppp_hhh_Pij.col(it1) << it1, it2;
                    //blockArrays_ppp_hhh_Pij.col(it2) << it2, it1;
                    //counter ++;
                }
                if (val1 == val_ik2){
                    blockArrays_ppp_hhh_Pik.col(it1) << it1, it2;
                    //blockArrays_ppp_hhh_Pik.col(it2) << it2, it1;
                    //counter ++;
                }
                if (val1 == val_jk2){
                    blockArrays_ppp_hhh_Pjk.col(it1) << it1, it2;
                    //blockArrays_ppp_hhh_Pjk.col(it2) << it2, it1;
                    //counter ++;
                }
                if (val1 == val_ijik2){
                    blockArrays_ppp_hhh_Pijik.col(it1) << it1, it2;
                    //blockArrays_ppp_hhh_Pijik.col(it2) << it2, it1;
                    //counter ++;
                }
                if (val1 == val_ijjk2){
                    blockArrays_ppp_hhh_Pijjk.col(it1) << it1, it2;
                    //blockArrays_ppp_hhh_Pijjk.col(it2) << it2, it1;
                    //counter ++;
                }
            }
            //if (counter != 5){ std::cout << "not all hole permutations found" << std::endl;}
        }

        //ppp permutations
        for (int it1 = range_lower_p; it1<range_upper_p; it1++){
            a1 = blockArrays_ppp_ppp(1,it1);
            b1 = blockArrays_ppp_ppp(2,it1);
            c1 = blockArrays_ppp_ppp(3,it1);

            val1 = a1 + b1*Np + c1*Np2;

            int counter = 0;
            //if indices are the same in the case of permutation
            /*if (a1 == b1 && a1 == c1){
                blockArrays_ppp_ppp_Pabc.col(it1) << it1, it1;
                counter ++;
            }
            if (a1 == b1){
                blockArrays_ppp_ppp_Pab.col(it1) << it1, it1;
                counter ++;
            }
            if (a1 == c1){
                blockArrays_ppp_ppp_Pac.col(it1) << it1, it1;
                counter ++;
            }
            if (b1 == c1){
                blockArrays_ppp_ppp_Pbc.col(it1) << it1, it1;
                counter ++;
            }*/

            //int n = omp_get_max_threads();
            #pragma omp parallel for num_threads(m_numThreads) private(a2,b2,c2, val_ab2,val_ac2,val_bc2,val_abac2,val_abbc2)
            for (int it2 = range_lower_p; it2<range_upper_p; it2++){  //starting at it1+1 means I'll never do the same permutation twice
                a2 = blockArrays_ppp_ppp(1,it2);
                b2 = blockArrays_ppp_ppp(2,it2);
                c2 = blockArrays_ppp_ppp(3,it2);

                val_ab2  = b2 + a2*Np + c2*Np2;
                val_ac2  = c2 + b2*Np + a2*Np2;
                val_bc2  = a2 + c2*Np + b2*Np2;
                val_abac2 = c2 + a2*Np + b2*Np2; // P(ab)P(ac): (a,b,c) -> (c,a,b)
                val_abbc2 = b2 + c2*Np + a2*Np2; // P(ab)P(bc): (a,b,c) -> (b,c,a)

                if (val1 == val_ab2){
                    blockArrays_ppp_ppp_Pab.col(it1) << it1, it2;
                    //blockArrays_ppp_ppp_Pab.col(it2) << it2, it1;
                    //if (it2-range_lower_p < 0){std::cout << it2-range_lower_p << std::endl;}
                    //counter ++;
                }
                if (val1 == val_ac2){
                    blockArrays_ppp_ppp_Pac.col(it1) << it1, it2;
                    //blockArrays_ppp_ppp_Pac.col(it2) << it2, it1;
                    //counter ++;
                }
                if (val1 == val_bc2){
                    blockArrays_ppp_ppp_Pbc.col(it1) << it1, it2;
                    //blockArrays_ppp_ppp_Pbc.col(it2) << it2, it1;
                    //counter ++;
                }
                if (val1 == val_abac2){
                    blockArrays_ppp_ppp_Pabac.col(it1) << it1, it2;
                    //blockArrays_ppp_ppp_Pabac.col(it2) << it2, it1;
                    //counter ++;
                }
                if (val1 == val_abbc2){
                    blockArrays_ppp_ppp_Pabbc.col(it1) << it1, it2;
                    //blockArrays_ppp_ppp_Pabbc.col(it2) << it2, it1;
                    //counter ++;
                }
            }
            //if (counter != 5){ std::cout << "not all particle permutations found" << std::endl;;}
        }
    }

    //std::cout << blockArrays_ppp_ppp_Pabbc.cols() << std::endl;

    /*for (int channel = 0; channel<5; channel++){
        range_lower_p = boundsHolder_hhhppp_ppp(0,channel);
        range_upper_p = boundsHolder_hhhppp_ppp(1,channel);
        for (int col=range_lower_p; col<range_upper_p; col++){
            if (blockArrays_ppp_ppp_Pab(0,col)-range_lower_p < 0){std::cout << blockArrays_ppp_ppp_Pab(0,col)-range_lower_p << std::endl;}
        }
    }*/

    //std::cout << blockArrays_ppp_ppp_Pab(1, 0) << std::endl;
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

void MakeIntMat::setTriples(bool argument){
    m_triplesOn = argument;
}

//0 means hole, 1 means particle
void MakeIntMat::mapper_1(std::vector<int>& sortVecIn, Eigen::MatrixXi& blockArraysIn, int i1, int s1){
    Eigen::MatrixXi   blockArrays_temp;
    int index       = 0;
    int colSize     = 0;

    bool cond_h = (i1 == 0);
    //bool cond_p = (i1 == 1);

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

    bool cond_hh = (i1 == 0 && i2 == 0);
    bool cond_hp = (i1 == 0 && i2 == 1);
    bool cond_ph = (i1 == 1 && i2 == 0);
    //bool cond_pp = (i1 == 1 && i2 == 1);

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

    bool cond_hhh = (i1 == 0 && i2 == 0 && i3==0);
    bool cond_hhp = (i1 == 0 && i2 == 0 && i3==1);
    bool cond_pph = (i1 == 1 && i2 == 1 && i3==0);
    bool cond_ppp = (i1 == 1 && i2 == 1 && i3==1);

    // 0 0 0
    if (cond_hhh){
        colSize = m_Nh*m_Nh*m_Nh;
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
        colSize = 10000;//colSize = (m_Ns-m_Nh)*(m_Ns-m_Nh)*m_Nh;
        //blockArrays_temp.conservativeResize(4, colSize);

        Eigen::MatrixXi   blockArrays_temp2;
        blockArrays_temp2.conservativeResize(4, colSize);

        int ku;

        for (int a=m_Nh; a<m_Ns; a++){
            for (int b=m_Nh; b<m_Ns; b++){
                for (int i=0; i<m_Nh; i++){
                    ku = m_system->kUnique3(a,b,i,s1,s2,s3);
                    auto it1 = std::find(sortVec_ppm_hhp.begin(), sortVec_ppm_hhp.end(), ku);
                    auto it2 = std::find(sortVec_p_h.begin(), sortVec_p_h.end(), ku);
                    auto it3 = std::find(sortVec_p_p.begin(), sortVec_p_p.end(), ku);
                    if (it1 != sortVec_ppp_hhh.end() || it2 != sortVec_p_h.end() || it3 != sortVec_p_p.end()){
                        blockArrays_temp2.col(index) << ku,a,b,i;
                        index += 1;

                        if (index >= colSize){
                            colSize += 10000;
                            blockArrays_temp2.conservativeResize(4, colSize);
                        }
                    }
                }
            }
        }

        std::cout << "made blockArrays_ppm_pph" << std::endl;

        blockArrays_temp.conservativeResize(4, index);
        for (int c=0; c<index; c++){
            blockArrays_temp.col(c) = blockArrays_temp2.col(c);
        }

        blockArrays_temp2.resize(0,0);

        /*for (int a=m_Nh; a<m_Ns; a++){
            for (int b=m_Nh; b<m_Ns; b++){
                for (int i=0; i<m_Nh; i++){
                    blockArrays_temp.col(index) << m_system->kUnique3(a,b,i,s1,s2,s3),a,b,i;
                    index += 1;
                }
            }
        }*/
    }
    // 1 1 1
    else if (cond_ppp){
        colSize = 10000;//(m_Ns-m_Nh)*(m_Ns-m_Nh)*(m_Ns-m_Nh);
        //blockArrays_temp.conservativeResize(4, colSize);

        Eigen::MatrixXi   blockArrays_temp2;
        blockArrays_temp2.conservativeResize(4, colSize);

        int ku;

        for (int a=m_Nh; a<m_Ns; a++){
            for (int b=m_Nh; b<m_Ns; b++){
                for (int c=m_Nh; c<m_Ns; c++){
                    ku = m_system->kUnique3(a,b,c,s1,s2,s3);
                    auto it = std::find(sortVec_ppp_hhh.begin(), sortVec_ppp_hhh.end(), ku);
                    if (it != sortVec_ppp_hhh.end()){
                        blockArrays_temp2.col(index) << ku,a,b,c;
                        index += 1;

                        if (index >= colSize){
                            colSize += 10000;
                            blockArrays_temp2.conservativeResize(4, colSize);
                        }
                    }
                }
            }
        }
        std::cout << "made blockArrays_ppp_ppp" << std::endl;

        blockArrays_temp.conservativeResize(4, index);
        for (int c=0; c<index; c++){
            blockArrays_temp.col(c) = blockArrays_temp2.col(c);
        }

        blockArrays_temp2.resize(0,0);

        /*for (int a=m_Nh; a<m_Ns; a++){
            for (int b=m_Nh; b<m_Ns; b++){
                for (int c=m_Nh; c<m_Ns; c++){
                    blockArrays_temp.col(index) << m_system->kUnique3(a,b,c,s1,s2,s3),a,b,c;
                    index += 1;
                }
            }
        }*/
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

void MakeIntMat::mapper_4(std::vector<int>& sortVecIn, Eigen::MatrixXi& blockArraysIn, int i1, int i2, int i3, int i4, int s1, int s2, int s3, int s4){

    Eigen::MatrixXi   blockArrays_temp1;
    Eigen::MatrixXi   blockArrays_temp2;
    int index       = 0;
    int colSize     = 0;

    bool cond_hhpp = (i1 == 0 && i2 == 0 && i3==1 && i4==1);
    bool cond_pphh = (i1 == 1 && i2 == 1 && i3==0 && i4==0);
    bool cond_hhhp = (i1 == 0 && i2 == 0 && i3==0 && i4==1);
    bool cond_ppph = (i1 == 1 && i2 == 1 && i3==1 && i4==0);

    //std::cout << "mapper 4 running with ";

    // 0 0 0 1
    if (cond_hhhp){
        //std::cout << "hhhp" << std::endl;
        colSize = 10000;//m_Nh*m_Nh*m_Nh*(m_Ns-m_Nh);
        blockArrays_temp1.conservativeResize(5, colSize);

        int ku;

        for (int i=0; i<m_Nh; i++){
            for (int j=0; j<m_Nh; j++){
                for (int k=0; k<m_Nh; k++){
                    for (int a=m_Nh; a<m_Ns; a++){
                        ku = m_system->kUnique4(i,j,k,a,s1,s2,s3,s4);
                        auto it = std::find(sortVec_pp_pp.begin(), sortVec_pp_pp.end(), ku);
                        if (it != sortVec_pp_pp.end()){
                            blockArrays_temp1.col(index) << ku,i,j,k,a;
                            index += 1;

                            if (index >= colSize){
                                colSize += 10000;
                                blockArrays_temp1.conservativeResize(5, colSize);
                            }
                        }
                    }
                }
            }
        }


        blockArrays_temp2.conservativeResize(5, index);
        for (int c=0; c<index; c++){
            blockArrays_temp2.col(c) = blockArrays_temp1.col(c);
        }

        blockArrays_temp1.resize(0,0);

        /*for (int i=0; i<m_Nh; i++){
            for (int j=0; j<m_Nh; j++){
                for (int k=0; k<m_Nh; k++){
                    for (int a=m_Nh; a<m_Ns; a++){
                        blockArrays_temp.col(index) << m_system->kUnique4(i,j,k,a,s1,s2,s3,s4),i,j,k,a;
                        index += 1;
                    }
                }
            }
        }*/
    }
    // 0 0 1 1
    if (cond_hhpp){ //this isn't used anymore after fix of T5a
        std::cout << "hhpp" << std::endl;
        colSize = 10000;//m_Nh*m_Nh*(m_Ns-m_Nh)*(m_Ns-m_Nh);
        blockArrays_temp1.conservativeResize(5, colSize);

        int ku;

        for (int i=0; i<m_Nh; i++){
            for (int j=0; j<m_Nh; j++){
                for (int a=m_Nh; a<m_Ns; a++){
                    for (int b=m_Nh; b<m_Ns; b++){
                        ku = m_system->kUnique4(i,j,a,b,s1,s2,s3,s4);
                        auto it = std::find(sortVec_pm_ph.begin(), sortVec_pm_ph.end(), ku);
                        if (it != sortVec_pm_ph.end()){
                            blockArrays_temp1.col(index) << ku,i,j,a,b;
                            index += 1;

                            if (index >= colSize){
                                colSize += 10000;
                                blockArrays_temp1.conservativeResize(5, colSize);
                            }
                        }
                    }
                }
            }
        }

        std::cout << "made blockArrays_ppmm_hhpp" << std::endl;

        blockArrays_temp2.conservativeResize(5, index);
        for (int c=0; c<index; c++){
            blockArrays_temp2.col(c) = blockArrays_temp1.col(c);
        }

        blockArrays_temp1.resize(0,0);

        /*for (int i=0; i<m_Nh; i++){
            for (int j=0; j<m_Nh; j++){
                for (int a=m_Nh; a<m_Ns; a++){
                    for (int b=m_Nh; b<m_Ns; b++){
                        blockArrays_temp.col(index) << m_system->kUnique4(i,j,a,b,s1,s2,s3,s4),i,j,a,b;
                        index += 1;
                    }
                }
            }
        }*/
    }
    // 1 1 0 0
    if (cond_pphh){
        std::cout << "pphh" << std::endl;
        colSize = 10000;//m_Nh*m_Nh*(m_Ns-m_Nh)*(m_Ns-m_Nh);
        blockArrays_temp1.conservativeResize(5, colSize);

        int ku;

        for (int a=m_Nh; a<m_Ns; a++){
            for (int b=m_Nh; b<m_Ns; b++){
                for (int i=0; i<m_Nh; i++){
                    for (int j=0; j<m_Nh; j++){
                        ku = m_system->kUnique4(a,b,i,j,s1,s2,s3,s4);
                        auto it = std::find(sortVec_pm_hp.begin(), sortVec_pm_hp.end(), ku);
                        if (it != sortVec_pm_hp.end()){
                            blockArrays_temp1.col(index) << ku,a,b,i,j;
                            index += 1;

                            if (index >= colSize){
                                colSize += 10000;
                                blockArrays_temp1.conservativeResize(5, colSize);
                            }
                        }
                    }
                }
            }
        }

        std::cout << "made blockArrays_ppmm_pphh" << std::endl;

        blockArrays_temp2.conservativeResize(5, index);
        for (int c=0; c<index; c++){
            blockArrays_temp2.col(c) = blockArrays_temp1.col(c);
        }

        blockArrays_temp1.resize(0,0);

        /*for (int i=0; i<m_Nh; i++){
            for (int j=0; j<m_Nh; j++){
                for (int a=m_Nh; a<m_Ns; a++){
                    for (int b=m_Nh; b<m_Ns; b++){
                        blockArrays_temp.col(index) << m_system->kUnique4(i,j,a,b,s1,s2,s3,s4),i,j,a,b;
                        index += 1;
                    }
                }
            }
        }*/
    }
    // 1 1 1 0
    else if (cond_ppph){
        std::cout << "ppph" << std::endl;
        colSize = 10000;//(m_Ns-m_Nh)*(m_Ns-m_Nh)*(m_Ns-m_Nh)*m_Nh;
        blockArrays_temp1.conservativeResize(5, colSize);

        int ku;

        for (int a=m_Nh; a<m_Ns; a++){
            for (int b=m_Nh; b<m_Ns; b++){
                for (int c=m_Nh; c<m_Ns; c++){
                    for (int i=0; i<m_Nh; i++){
                        ku = m_system->kUnique4(a,b,c,i,s1,s2,s3,s4);
                        auto it = std::find(sortVec_pp_hh.begin(), sortVec_pp_hh.end(), ku);
                        if (it != sortVec_pp_hh.end()){
                            blockArrays_temp1.col(index) << ku,a,b,c,i;
                            index += 1;

                            if (index >= colSize){
                                colSize += 10000;
                                blockArrays_temp1.conservativeResize(5, colSize);
                            }
                        }
                    }
                }
            }
        }


        blockArrays_temp2.conservativeResize(5, index);
        for (int c=0; c<index; c++){
            blockArrays_temp2.col(c) = blockArrays_temp1.col(c);
        }

        blockArrays_temp1.resize(0,0);

        std::cout << "made blockArrays_pppm_ppph" << std::endl;

        /*for (int a=m_Nh; a<m_Ns; a++){
            for (int b=m_Nh; b<m_Ns; b++){
                for (int c=m_Nh; c<m_Ns; c++){
                    for (int i=0; i<m_Nh; i++){
                        blockArrays_temp.col(index) << m_system->kUnique4(a,b,c,i,s1,s2,s3,s4),a,b,c,i;
                        index += 1;
                    }
                }
            }
        }*/
    }

    Eigen::MatrixXi blockArrays = blockArrays_temp2;     //need this
    Eigen::VectorXi veryTempVec = blockArrays_temp2.row(0);
    std::vector<int> sortVec(veryTempVec.data(), veryTempVec.row(0).data() + blockArrays_temp2.cols());
    std::vector<int> tempVec = sortVec;

    index = 0;
    for (auto i: sort_indexes(tempVec)){
        blockArrays.col(index) = blockArrays_temp2.col(i);
        sortVec[index]         = tempVec[i];
        index += 1;
    }

    sortVec.erase( unique( sortVec.begin(), sortVec.end() ), sortVec.end() ); //need this

    blockArraysIn = blockArrays;
    sortVecIn     = sortVec;
}


//this mapper is distinct form the previous because it otherwise completely hogs the RAM
void MakeIntMat::mapper_5(std::vector<int>& sortVecIn, Eigen::MatrixXi& blockArraysIn, int i1, int i2, int i3, int i4, int i5, int s1, int s2, int s3, int s4, int s5){

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Status status;

    Eigen::MatrixXi   blockArrays_temp1;
    Eigen::MatrixXi   blockArrays_temp2;
    int index       = 0;
    int colSize     = 0;

    bool cond_hhhpp = (i1 == 0 && i2 == 0 && i3 == 0 && i4==1 && i5==1);
    bool cond_ppphh = (i1 == 1 && i2 == 1 && i3 == 1 && i4==0 && i5==0);

    //std::cout << world_rank << std::endl;
    //std::cout << "mapper 5 running with ";

    // 0 0 0 1 1
    if (cond_hhhpp){
        if (world_rank==0){std::cout << "hhhpp" << std::endl;}
        colSize = 10000;//m_Nh*m_Nh*m_Nh*(m_Ns-m_Nh)*(m_Ns-m_Nh);
        blockArrays_temp1.conservativeResize(6, colSize);

        int ku;

        int dist = (m_Ns-m_Nh)/world_size;
        int remains = (m_Ns-m_Nh) % world_size;
        bool test = (world_rank == world_size-1);

        for (int i=0; i<m_Nh; i++){
            for (int j=0; j<m_Nh; j++){
                for (int k=0; k<m_Nh; k++){
                    for (int a=m_Nh + dist*world_rank; a<m_Nh + dist*(world_rank+1) + remains*test; a++){
                        for (int b=m_Nh; b<m_Ns; b++){
                            ku = m_system->kUnique5(i,j,k,a,b,s1,s2,s3,s4,s5);
                            auto it = std::find(sortVec_p_p.begin(), sortVec_p_p.end(), ku);
                            if (it != sortVec_p_p.end()){
                                blockArrays_temp1.col(index) << ku,i,j,k,a,b;
                                index += 1;

                                if (index >= colSize){
                                    colSize += 10000;
                                    blockArrays_temp1.conservativeResize(6, colSize);
                                }
                            }
                        }
                    }
                }
            }
        }
        //std::cout << "sup" << std::endl;

        blockArrays_temp2.conservativeResize(6, index);
        for (int c=0; c<index; c++){
            blockArrays_temp2.col(c) = blockArrays_temp1.col(c);
        }

        blockArrays_temp1.resize(0,0);

        //std::cout << "made blockArrays_pppmm_hhhpp" << std::endl;

        /*for (int i=0; i<m_Nh; i++){
            for (int j=0; j<m_Nh; j++){
                for (int k=0; k<m_Nh; k++){
                    for (int a=m_Nh; a<m_Ns; a++){
                        for (int b=m_Nh; b<m_Ns; b++){
                            blockArrays_temp.col(index) << m_system->kUnique5(i,j,k,a,b,s1,s2,s3,s4,s5),i,j,k,a,b;
                            index += 1;
                        }
                    }
                }
            }
        }*/
    }
    // 1 1 1 0 0
    else if (cond_ppphh){
        if (world_rank==0){std::cout << "ppphh" << std::endl;}
        colSize = 10000;//(m_Ns-m_Nh)*(m_Ns-m_Nh)*(m_Ns-m_Nh)*m_Nh*m_Nh;
        blockArrays_temp1.conservativeResize(6, colSize);

        int ku;

        //std::cout << "ppphh iterator start" << std::endl;

        int dist = (m_Ns-m_Nh)/world_size;
        int remains = (m_Ns-m_Nh) % world_size;
        bool test = (world_rank == world_size-1);

        //int sumA = 0;

        for (int a=m_Nh + dist*world_rank; a<m_Nh + dist*(world_rank+1) + remains*test; a++){
            //sumA ++;
            for (int b=m_Nh; b<m_Ns; b++){
                for (int c=m_Nh; c<m_Ns; c++){
                    for (int i=0; i<m_Nh; i++){
                        for (int j=0; j<m_Nh; j++){
                            ku = m_system->kUnique5(a,b,c,i,j,s1,s2,s3,s4,s5);
                            auto it = std::find(sortVec_p_h.begin(), sortVec_p_h.end(), ku);
                            if (it != sortVec_p_h.end()){
                                //std::cout << ku << std::endl;
                                blockArrays_temp1.col(index) << ku,a,b,c,i,j;
                                index += 1;

                                if (index >= colSize){
                                    colSize += 10000;
                                    blockArrays_temp1.conservativeResize(6, colSize);
                                }
                            }
                        }
                    }
                }
            }
        }

        //std::cout << sumA << std::endl;

        //std::cout << "made blockArrays_pppmm_ppphh" << std::endl;

        //std::cout << "ppphh iterator end" << std::endl;

        blockArrays_temp2.conservativeResize(6, index);
        for (int c=0; c<index; c++){
            blockArrays_temp2.col(c) = blockArrays_temp1.col(c);
        }

        blockArrays_temp1.resize(0,0);

        /*for (int a=m_Nh; a<m_Ns; a++){
            for (int b=m_Nh; b<m_Ns; b++){
                for (int c=m_Nh; c<m_Ns; c++){
                    for (int i=0; i<m_Nh; i++){
                        for (int j=0; j<m_Nh; j++){
                            blockArrays_temp.col(index) << m_system->kUnique5(a,b,c,i,j,s1,s2,s3,s4,s5),a,b,c,i,j;
                            index += 1;
                        }
                    }
                }
            }
        }*/
    }

    int cols = blockArrays_temp2.cols();
    int total_columns;
    MPI_Reduce(&cols, &total_columns, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    int displs[world_size];
    int recvCount[world_size];

    int sendCount = blockArrays_temp2.cols()*blockArrays_temp2.rows();


    MPI_Gather(&sendCount, 1, MPI_INT, recvCount, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (world_rank==0){
        displs[0] = 0;
        for (int rank = 1; rank < world_size; rank++){
            displs[rank] = displs[rank-1]+recvCount[rank-1];
        }
    }

    Eigen::MatrixXi blockArrays_recv;
    if (world_rank==0){
        blockArrays_recv.conservativeResize(6, total_columns);
    }

    MPI_Gatherv(blockArrays_temp2.data(), sendCount, MPI_INT, blockArrays_recv.data(), recvCount, displs, MPI_INT, 0, MPI_COMM_WORLD);


    if(world_rank == 0){
        Eigen::MatrixXi blockArrays = blockArrays_recv;

        Eigen::VectorXi veryTempVec = blockArrays_recv.row(0);
        std::vector<int> sortVec(veryTempVec.data(), veryTempVec.row(0).data() + blockArrays_recv.cols());
        std::vector<int> tempVec = sortVec;

        //std::cout << "making sortVec" << std::endl;
        index = 0;
        for (auto i: sort_indexes(tempVec)){
            blockArrays.col(index) = blockArrays_recv.col(i);
            sortVec[index]         = tempVec[i];
            index += 1;
        }

        sortVec.erase( unique( sortVec.begin(), sortVec.end() ), sortVec.end() ); //need this

        blockArraysIn = blockArrays;
        sortVecIn     = sortVec;
    }

    //blockArrays.conservativeResize(0,0);

    /*int size; if(world_rank==0){size = total_columns*5;}
    MPI_Bcast( size, 1, MPI_INT, 0, MPI_COMM_WORLD );
    MPI_Bcast( blockArrays.data(), size, MPI_INT, 0, MPI_COMM_WORLD );*/
}


// ##################################################
// ##                                              ##
// ## Make blocks                                  ##
// ##                                              ##
// ##################################################

Eigen::MatrixXd MakeIntMat::I1_makemat(int channel1, int channel2){    //makes a 2x2 matrix

    int range_lower1 = indexHolder_pp_hh(0, channel1);
    int range_upper1 = indexHolder_pp_hh(1, channel1);
    int range_lower2 = indexHolder_pp_pp(0, channel2);
    int range_upper2 = indexHolder_pp_pp(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int k; int l;
    int c; int d;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            /*
            k = blockArrays_pp_hh(1,i1);
            l = blockArrays_pp_hh(2,i1);
            c = blockArrays_pp_pp(1,i2);
            d = blockArrays_pp_pp(2,i2);
            id = Identity_hhpp(k,l,c,d);
            */

            id = Identity_hhpp(blockArrays_pp_hh(1,i1),
                               blockArrays_pp_hh(2,i1),
                               blockArrays_pp_pp(1,i2),
                               blockArrays_pp_pp(2,i2));
            returnMat(i1-range_lower1, i2-range_lower2) = Vhhpp_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeIntMat::I2_makemat(int channel1, int channel2){    //makes a 2x2 matrix

    int range_lower1 = indexHolder_pm_hp(0, channel1);
    int range_upper1 = indexHolder_pm_hp(1, channel1);
    int range_lower2 = indexHolder_pm_ph(0, channel2);
    int range_upper2 = indexHolder_pm_ph(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int k; int l;
    int c; int d;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            /*
            k = blockArrays_pm_hp(1,i1);
            c = blockArrays_pm_hp(2,i1);
            d = blockArrays_pm_ph(1,i2);
            l = blockArrays_pm_ph(2,i2);
            id = Identity_hhpp(k,l,c,d);
            */

            id = Identity_hhpp(blockArrays_pm_hp(1,i1),
                               blockArrays_pm_ph(2,i2),
                               blockArrays_pm_hp(2,i1),
                               blockArrays_pm_ph(1,i2));

            returnMat(i1-range_lower1, i2-range_lower2) = Vhhpp_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeIntMat::I3_makemat(int channel1, int channel2){    //makes a 3x1 matrix

    int range_lower1 = indexHolder_ppm_pph(0, channel1);
    int range_upper1 = indexHolder_ppm_pph(1, channel1);
    int range_lower2 = indexHolder_p_h(0, channel2);
    int range_upper2 = indexHolder_p_h(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int k; int l;
    int c; int d;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            /*
            c = blockArrays_ppm_pph(1,i1);
            d = blockArrays_ppm_pph(2,i1);
            k = blockArrays_ppm_pph(3,i1);
            l = blockArrays_p_h(1,i2);
            id = Identity_hhpp(k,l,c,d);
            */

            id = Identity_hhpp(blockArrays_ppm_pph(3,i1),
                               blockArrays_p_h(1,i2),
                               blockArrays_ppm_pph(1,i1),
                               blockArrays_ppm_pph(2,i1));

            returnMat(i1-range_lower1, i2-range_lower2) = Vhhpp_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeIntMat::I4_makemat(int channel1, int channel2){    //makes a 3x1 matrix

    int range_lower1 = indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = indexHolder_p_p(0, channel2);
    int range_upper2 = indexHolder_p_p(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int k; int l;
    int c; int d;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            /*
            k = blockArrays_ppm_hhp(1,i1);
            l = blockArrays_ppm_hhp(2,i1);
            c = blockArrays_ppm_hhp(3,i1);
            d = blockArrays_p_p(1,i2);
            id = Identity_hhpp(k,l,c,d);
            */

            id = Identity_hhpp(blockArrays_ppm_hhp(1,i1),
                               blockArrays_ppm_hhp(2,i1),
                               blockArrays_ppm_hhp(3,i1),
                               blockArrays_p_p(1,i2));

            returnMat(i1-range_lower1, i2-range_lower2) = Vhhpp_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeIntMat::D10b_makemat(int channel1, int channel2){    //makes a 3x1 matrix

    int range_lower1 = indexHolder_ppm_pph(0, channel1);
    int range_upper1 = indexHolder_ppm_pph(1, channel1);
    int range_lower2 = indexHolder_p_p(0, channel2);
    int range_upper2 = indexHolder_p_p(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int a; int k;
    int c; int d;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            c = blockArrays_ppm_pph(1,i1);
            d = blockArrays_ppm_pph(2,i1);
            k = blockArrays_ppm_pph(3,i1);
            a = blockArrays_p_p(1,i2);

            id = Identity_ppph(c,d,a,k);
            returnMat(i1-range_lower1, i2-range_lower2) = Vppph_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeIntMat::D10c_makemat(int channel1, int channel2){

    int range_lower1 = indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = indexHolder_p_h(0, channel2);
    int range_upper2 = indexHolder_p_h(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int k; int l;
    int c; int j;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            k = blockArrays_ppm_hhp(1,i1);
            l = blockArrays_ppm_hhp(2,i1);
            c = blockArrays_ppm_hhp(3,i1);
            j = blockArrays_p_h(1,i2);

            id = Identity_hhhp(k,l,j,c);
            returnMat(i1-range_lower1, i2-range_lower2) = Vhhhp_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeIntMat::T1a_makemat(int channel1, int channel2){    //makes a 3x1 matrix

    int range_lower1 = indexHolder_ppm_pph(0, channel1);
    int range_upper1 = indexHolder_ppm_pph(1, channel1);
    int range_lower2 = indexHolder_p_p(0, channel2);
    int range_upper2 = indexHolder_p_p(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int b; int c;
    int d; int k;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            b = blockArrays_ppm_pph(1,i1);
            c = blockArrays_ppm_pph(2,i1);
            k = blockArrays_ppm_pph(3,i1);
            d = blockArrays_p_p(1,i2);

            id = Identity_ppph(b,c,d,k);
            returnMat(i1-range_lower1, i2-range_lower2) = Vppph_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeIntMat::T1b_makemat(int channel1, int channel2){    //makes a 3x1 matrix

    int range_lower1 = indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = indexHolder_p_h(0, channel2);
    int range_upper2 = indexHolder_p_h(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int j; int k;
    int l; int c;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        for (int i2 = range_lower2; i2<range_upper2; i2++){

            j = blockArrays_ppm_hhp(1,i1);
            k = blockArrays_ppm_hhp(2,i1);
            c = blockArrays_ppm_hhp(3,i1);
            l = blockArrays_p_h(1,i2);

            id = Identity_hhhp(j,k,l,c);
            returnMat(i1-range_lower1, i2-range_lower2) = Vhhhp_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeIntMat::T3b_makemat(int channel1, int channel2){    //makes a 2x2 matrix

    int range_lower1 = indexHolder_pm_pp(0, channel1);
    int range_upper1 = indexHolder_pm_pp(1, channel1);
    int range_lower2 = indexHolder_pm_hp(0, channel2);
    int range_upper2 = indexHolder_pm_hp(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int d; int e;
    int b; int l;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        e = blockArrays_pm_pp(1,i1);
        b = blockArrays_pm_pp(2,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            l = blockArrays_pm_hp(1,i2);
            d = blockArrays_pm_hp(2,i2);

            id = Identity_ppph(d,e,b,l);
            returnMat(i1-range_lower1, i2-range_lower2) = Vppph_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeIntMat::T3c_makemat(int channel1, int channel2){    //makes a 2x2 matrix

    int range_lower1 = indexHolder_pm_hh(0, channel1);
    int range_upper1 = indexHolder_pm_hh(1, channel1);
    int range_lower2 = indexHolder_pm_ph(0, channel2);
    int range_upper2 = indexHolder_pm_ph(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int l; int m;
    int j; int d;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        m = blockArrays_pm_hh(1,i1);
        j = blockArrays_pm_hh(2,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            d = blockArrays_pm_ph(1,i2);
            l = blockArrays_pm_ph(2,i2);

            id = Identity_hhhp(l,m,j,d);
            returnMat(i1-range_lower1, i2-range_lower2) = Vhhhp_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeIntMat::T3d_makemat(int channel1, int channel2){    //makes a 2x2 matrix

    int range_lower1 = indexHolder_pp_pp(0, channel1);
    int range_upper1 = indexHolder_pp_pp(1, channel1);
    int range_lower2 = indexHolder_pp_ph(0, channel2);
    int range_upper2 = indexHolder_pp_ph(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int d; int e;
    int c; int l;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        d = blockArrays_pp_pp(1,i1);
        e = blockArrays_pp_pp(2,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            c = blockArrays_pp_ph(1,i2);
            l = blockArrays_pp_ph(2,i2);

            id = Identity_ppph(d,e,c,l);
            returnMat(i1-range_lower1, i2-range_lower2) = Vppph_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeIntMat::T3e_makemat(int channel1, int channel2){    //makes a 3x1 matrix

    int range_lower1 = indexHolder_ppm_hhh(0, channel1);
    int range_upper1 = indexHolder_ppm_hhh(1, channel1);
    int range_lower2 = indexHolder_p_p(0, channel2);
    int range_upper2 = indexHolder_p_p(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int l; int m;
    int k; int d;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        l = blockArrays_ppm_hhh(1,i1);
        m = blockArrays_ppm_hhh(2,i1);
        k = blockArrays_ppm_hhh(3,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            d = blockArrays_p_p(1,i2);

            id = Identity_hhhp(l,m,k,d);
            returnMat(i1-range_lower1, i2-range_lower2) = Vhhhp_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeIntMat::T5a_makemat(int channel1, int channel2){    //makes a 2x2 matrix

    int range_lower1 = indexHolder_pm_hp(0, channel1);
    int range_upper1 = indexHolder_pm_hp(1, channel1);
    int range_lower2 = indexHolder_pm_ph(0, channel2);
    int range_upper2 = indexHolder_pm_ph(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int l; int m;
    int d; int e;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        m = blockArrays_pm_hp(1,i1);
        e = blockArrays_pm_hp(2,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            d = blockArrays_pm_ph(1,i2);
            l = blockArrays_pm_ph(2,i2);

            id = Identity_hhpp(l,m,d,e);
            returnMat(i1-range_lower1, i2-range_lower2) = Vhhpp_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeIntMat::T5b_makemat(int channel1, int channel2){    //makes a 3x1 matrix

    int range_lower1 = indexHolder_ppm_pph(0, channel1);
    int range_upper1 = indexHolder_ppm_pph(1, channel1);
    int range_lower2 = indexHolder_p_h(0, channel2);
    int range_upper2 = indexHolder_p_h(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int l; int m;
    int d; int e;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        d = blockArrays_ppm_pph(1,i1);
        e = blockArrays_ppm_pph(2,i1);
        l = blockArrays_ppm_pph(3,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            m = blockArrays_p_h(1,i2);

            id = Identity_hhpp(l,m,d,e);
            returnMat(i1-range_lower1, i2-range_lower2) = Vhhpp_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeIntMat::T5c_makemat(int channel1, int channel2){    //makes a 3x1 matrix

    int range_lower1 = indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = indexHolder_p_p(0, channel2);
    int range_upper2 = indexHolder_p_p(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int l; int m;
    int d; int e;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        l = blockArrays_ppm_hhp(1,i1);
        m = blockArrays_ppm_hhp(2,i1);
        d = blockArrays_ppm_hhp(3,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            e = blockArrays_p_p(1,i2);

            id = Identity_hhpp(l,m,d,e);
            returnMat(i1-range_lower1, i2-range_lower2) = Vhhpp_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeIntMat::T5d_makemat(int channel1, int channel2){    //makes a 3x1 matrix

    int range_lower1 = indexHolder_ppm_hhp(0, channel1);
    int range_upper1 = indexHolder_ppm_hhp(1, channel1);
    int range_lower2 = indexHolder_p_p(0, channel2);
    int range_upper2 = indexHolder_p_p(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int l; int m;
    int d; int e;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        l = blockArrays_ppm_hhp(1,i1);
        m = blockArrays_ppm_hhp(2,i1);
        e = blockArrays_ppm_hhp(3,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            d = blockArrays_p_p(1,i2);

            id = Identity_hhpp(l,m,d,e);
            returnMat(i1-range_lower1, i2-range_lower2) = Vhhpp_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeIntMat::T5e_makemat(int channel1, int channel2){    //makes a 3x1 matrix

    int range_lower1 = indexHolder_ppm_pph(0, channel1);
    int range_upper1 = indexHolder_ppm_pph(1, channel1);
    int range_lower2 = indexHolder_p_h(0, channel2);
    int range_upper2 = indexHolder_p_h(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int l; int m;
    int d; int e;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        d = blockArrays_ppm_pph(1,i1);
        e = blockArrays_ppm_pph(2,i1);
        m = blockArrays_ppm_pph(3,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            l = blockArrays_p_h(1,i2);

            id = Identity_hhpp(l,m,d,e);
            returnMat(i1-range_lower1, i2-range_lower2) = Vhhpp_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeIntMat::T5f_makemat(int channel1, int channel2){    //makes a 2x2 matrix

    int range_lower1 = indexHolder_pp_hh(0, channel1);
    int range_upper1 = indexHolder_pp_hh(1, channel1);
    int range_lower2 = indexHolder_pp_pp(0, channel2);
    int range_upper2 = indexHolder_pp_pp(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int l; int m;
    int d; int e;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        l = blockArrays_pp_hh(1,i1);
        m = blockArrays_pp_hh(2,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            d = blockArrays_pp_pp(1,i2);
            e = blockArrays_pp_pp(2,i2);

            id = Identity_hhpp(l,m,d,e);
            returnMat(i1-range_lower1, i2-range_lower2) = Vhhpp_elements[id];
        }
    }

    return returnMat;
}

Eigen::MatrixXd MakeIntMat::T5g_makemat(int channel1, int channel2){    //makes a 2x2 matrix

    int range_lower1 = indexHolder_pp_hh(0, channel1);
    int range_upper1 = indexHolder_pp_hh(1, channel1);
    int range_lower2 = indexHolder_pp_pp(0, channel2);
    int range_upper2 = indexHolder_pp_pp(1, channel2);

    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(range_upper1 - range_lower1, range_upper2 - range_lower2);


    int id;
    int l; int m;
    int d; int e;

    for (int i1 = range_lower1; i1<range_upper1; i1++){
        l = blockArrays_pp_hh(1,i1);
        m = blockArrays_pp_hh(2,i1);
        for (int i2 = range_lower2; i2<range_upper2; i2++){
            d = blockArrays_pp_pp(1,i2);
            e = blockArrays_pp_pp(2,i2);

            id = Identity_hhpp(l,m,d,e);
            returnMat(i1-range_lower1, i2-range_lower2) = Vhhpp_elements[id];
        }
    }

    return returnMat;
}



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
        blockArrays1_pointer = blockArrays_ppm_hhp;
        sortVec1             = sortVec_ppm_hhp;
        indexHolder1_pointer = indexHolder_ppm_hhp;
    }
    // 1 1 0
    else if (cond_pph){
        blockArrays1_pointer = blockArrays_ppm_pph;
        sortVec1             = sortVec_ppm_pph;
        indexHolder1_pointer = indexHolder_ppm_pph;
    }

    // 0
    if (cond_h){
        blockArrays2_pointer = blockArrays_p_h;
        sortVec2             = sortVec_p_h;
        indexHolder2_pointer = indexHolder_p_h;
    }
    // 1
    else if (cond_p){
        blockArrays2_pointer = blockArrays_p_p;
        sortVec2             = sortVec_p_p;
        indexHolder2_pointer = indexHolder_p_p;
    }

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

    //cout << indexHolder1_pointer.cols() << endl;

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
    else if (cond_hhp && cond_h){
        for (int i = range_lower1; i<range_upper1; i++){
            for (int j = range_lower2; j<range_upper2; j++){
                id = Identity_hhhp((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays2_pointer)(1,j), (blockArrays1_pointer)(3,i));
                returnMat(i-range_lower1, j-range_lower2) = Vhhhp_elements[id];
            }
        }
    }
    else if (cond_pph && cond_p){
        for (int i = range_lower1; i<range_upper1; i++){
            for (int j = range_lower2; j<range_upper2; j++){
                //cout << "sup" << endl;
                id = Identity_ppph((blockArrays1_pointer)(1,i), (blockArrays1_pointer)(2,i), (blockArrays2_pointer)(1,j), (blockArrays1_pointer)(3,i));
                returnMat(i-range_lower1, j-range_lower2) = Vppph_elements[id];
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
        blockArrays1_pointer = blockArrays_pp_hh;
        sortVec1             = sortVec_pp_hh;
        indexHolder1_pointer = indexHolder_pp_hh;
    }
    // 0 1
    else if (cond_hp1){
        blockArrays1_pointer = blockArrays_pm_hp;
        sortVec1             = sortVec_pm_hp;
        indexHolder1_pointer = indexHolder_pm_hp;
    }
    // 1 0
    else if (cond_ph1){
        blockArrays1_pointer = blockArrays_pm_ph;
        sortVec1             = sortVec_pm_ph;
        indexHolder1_pointer = indexHolder_pm_ph;
    }
    // 1 1
    else {
        blockArrays1_pointer = blockArrays_pp_pp;
        sortVec1             = sortVec_pp_pp;
        indexHolder1_pointer = indexHolder_pp_pp;
    }

    // 0 0
    if (cond_hh2){
        blockArrays2_pointer = blockArrays_pp_hh;
        sortVec2             = sortVec_pp_hh;
        indexHolder2_pointer = indexHolder_pp_hh;
    }
    // 0 1
    else if (cond_hp2){
        blockArrays2_pointer = blockArrays_pm_hp;
        sortVec2             = sortVec_pm_hp;
        indexHolder2_pointer = indexHolder_pm_hp;
    }
    // 1 0
    else if (cond_ph2){
        blockArrays2_pointer = blockArrays_pm_ph;
        sortVec2             = sortVec_pm_ph;
        indexHolder2_pointer = indexHolder_pm_ph;
    }
    // 1 1
    else {
        blockArrays2_pointer = blockArrays_pp_pp;
        sortVec2             = sortVec_pp_pp;
        indexHolder2_pointer = indexHolder_pp_pp;
    }

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
        //int count1 = 0;
        //int count2 = 0;
        for (int i = range_lower1; i<range_upper1; i++){
            for (int j = range_lower2; j<range_upper2; j++){
                //count1 ++;
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
        //int count1 = 0;
        //int count2 = 0;
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

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    m_system = system;
    m_Nh     = Nh;
    m_Ns     = Ns;
    m_Np     = Ns-Nh;

    //we use these in the indentity functions (less calculations)
    m_Nh2    = Nh*Nh;
    m_Nh3    = Nh*Nh*Nh;
    m_Nh4    = Nh*Nh*Nh*Nh;
    m_Nh5    = Nh*Nh*Nh*Nh*Nh;
    m_NhNs2  = Nh*Ns*Ns;
    m_Nh2Ns  = Nh*Nh*Ns;
    m_Nh3Ns  = Nh*Nh*Nh*Ns;
    m_Nh3Ns2 = Nh*Nh*Nh*Ns*Ns;

    //I only bothered to let the mapper functions take on the values shown here
    //Any other blockArrays you want to make will have to be implemented by yourself (it's not hard)

    mapper_1(sortVec_p_h, blockArrays_p_h, 0, 1);        //h
    mapper_1(sortVec_p_p, blockArrays_p_p, 1, 1);        //p

    mapper_2(sortVec_pp_hh, blockArrays_pp_hh, 0,0, 1,1);      //hh
    mapper_2(sortVec_pm_hp, blockArrays_pm_hp, 0,1, 1,-1);     //hp
    mapper_2(sortVec_pm_ph, blockArrays_pm_ph, 1,0, 1,-1);     //ph
    mapper_2(sortVec_pp_pp, blockArrays_pp_pp, 1,1, 1,1);      //pp
    mapper_hp();        //hp for Vhphp

    mapper_3(sortVec_ppm_hhp, blockArrays_ppm_hhp, 0,0,1, 1,1,-1);    //hhp
    mapper_3(sortVec_ppm_pph, blockArrays_ppm_pph, 1,1,0, 1,1,-1);    //pph

    if (m_triplesOn){
        mapper_2(sortVec_pp_ph, blockArrays_pp_ph, 1,0, 1,1);             //ph
        mapper_2(sortVec_pp_hp, blockArrays_pp_hp, 0,1, 1,1);             //hp
        mapper_2(sortVec_pm_hh, blockArrays_pm_hh, 0,0, 1,-1);            //ph
        mapper_2(sortVec_pm_pp, blockArrays_pm_pp, 1,1, 1,-1);            //hp
        mapper_3(sortVec_ppp_hhh, blockArrays_ppp_hhh, 0,0,0, 1,1,1);     //hhh
        mapper_3(sortVec_ppm_hhh, blockArrays_ppm_hhh, 0,0,0, 1,1,-1);    //hhh
        mapper_3(sortVec_ppp_ppp, blockArrays_ppp_ppp, 1,1,1, 1,1,1);     //ppp
        if (world_rank==0){"finished ppp";}
        mapper_4(sortVec_pppm_hhhp, blockArrays_pppm_hhhp, 0,0,0,1, 1,1,1,-1);     //hhhp
        if (world_rank==0){"finished hhhp";}
        //mapper_4(sortVec_ppmm_hhpp, blockArrays_ppmm_hhpp, 0,0,1,1, 1,1,-1,-1);    //hhpp
        //if (world_rank==0){"finished hhpp";}
        mapper_4(sortVec_ppmm_pphh, blockArrays_ppmm_pphh, 1,1,0,0, 1,1,-1,-1);    //pphh
        if (world_rank==0){"finished pphh";}
        mapper_4(sortVec_pppm_ppph, blockArrays_pppm_ppph, 1,1,1,0, 1,1,1,-1);     //ppph
        if (world_rank==0){"finished ppph";}

        mapper_5(sortVec_pppmm_hhhpp, blockArrays_pppmm_hhhpp, 0,0,0,1,1, 1,1,1,-1,-1);     //hhhpp
        if (world_rank==0){"finished hhppp";}
        mapper_5(sortVec_pppmm_ppphh, blockArrays_pppmm_ppphh, 1,1,1,0,0, 1,1,1,-1,-1);     //ppphh
        if (world_rank==0){"finished ppphh";}

    }

    cout << "made all blockArrays" << endl;

    int counter         = 0;

    int range_lower     = 0;
    int range_upper     = 0;

    int range_lower_hh  = 0;
    int range_upper_hh  = 0;

    int range_lower_hp  = 0;
    int range_upper_hp  = 0;

    int range_lower_ph = 0;
    int range_upper_ph = 0;

    int range_lower_pp  = 0;
    int range_upper_pp  = 0;

    int range_lower_hhh  = 0;
    int range_upper_hhh  = 0;

    int range_lower_ppp  = 0;
    int range_upper_ppp  = 0;


    //set indexHolders
    //I refrain from doing it for Vpppp for now, since I see no purpose as of yet for doing so.
    boundsHolder_hhpp_hh.conservativeResize(2, Eigen::NoChange);
    boundsHolder_hhpp_pp.conservativeResize(2, Eigen::NoChange);

    if (m_triplesOn){
        boundsHolder_hhhppp_hhh.conservativeResize(2, Eigen::NoChange);
        boundsHolder_hhhppp_ppp.conservativeResize(2, Eigen::NoChange);
    }

    indexHolder_p_h.conservativeResize(2,sortVec_p_h.size());
    for (int h=0; h<sortVec_p_h.size(); h++){
        int val_h = sortVec_p_h[h];
        for (int bA_h=0; bA_h<m_Nh; bA_h++){
            if (val_h == blockArrays_p_h(0,bA_h)){
                range_upper = bA_h+1;
                counter += 1;
            }
        }
        range_lower = range_upper - counter;
        counter = 0;
        indexHolder_p_h.col(h) << range_lower, range_upper; //this now has same indices as sortVec
    }

    counter     = 0;
    range_lower = 0;
    range_upper = 0;
    indexHolder_p_p.conservativeResize(2,sortVec_p_p.size());
    for (int p=0; p<sortVec_p_p.size(); p++){
        int val_p = sortVec_p_p[p];
        for (int bA_p=0; bA_p<(m_Ns-m_Nh); bA_p++){
            if (val_p == blockArrays_p_p(0,bA_p)){
                range_upper = bA_p+1;
                counter += 1;
            }
        }
        range_lower = range_upper - counter;
        counter = 0;
        indexHolder_p_p.col(p) << range_lower, range_upper; //this now has same indices as sortVec
    }

    if (world_rank==0){cout << "made one-particle indexHolders" << endl;}

    counter     = 0;
    range_lower = 0;
    range_upper = 0;
    indexHolder_pp_hh.conservativeResize(2,sortVec_pp_hh.size());
    for (int hh=0; hh<sortVec_pp_hh.size(); hh++){
        int val_hh = sortVec_pp_hh[hh];
        for (int bA_hh=0; bA_hh<m_Nh*m_Nh; bA_hh++){
            if (val_hh == blockArrays_pp_hh(0,bA_hh)){
                range_upper = bA_hh+1;
                counter += 1;
            }
        }
        range_lower = range_upper - counter;
        counter = 0;
        indexHolder_pp_hh.col(hh) << range_lower, range_upper; //this now has same indices as sortVec
    }

    counter     = 0;
    range_lower = 0;
    range_upper = 0;
    indexHolder_pm_hp.conservativeResize(2,sortVec_pm_hp.size());
    for (int hp=0; hp<sortVec_pm_hp.size(); hp++){
        int val_hp = sortVec_pm_hp[hp];
        auto it = std::find(sortVec_pm_ph.begin(), sortVec_pm_ph.end(), val_hp);
        if (it != sortVec_pm_ph.end()){
            for (int bA_hp=0; bA_hp<m_Nh*(m_Ns-m_Nh); bA_hp++){
                if (val_hp == blockArrays_pm_hp(0,bA_hp)){
                    range_upper = bA_hp+1;
                    counter += 1;
                }
            }
        }
        range_lower = range_upper - counter;
        counter = 0;
        indexHolder_pm_hp.col(hp) << range_lower, range_upper; //this now has same indices as sortVec
    }

    counter     = 0;
    range_lower = 0;
    range_upper = 0;
    indexHolder_pm_ph.conservativeResize(2,sortVec_pm_ph.size());
    for (int ph=0; ph<sortVec_pm_ph.size(); ph++){
        int val_ph = sortVec_pm_ph[ph];
        auto it1 = std::find(sortVec_pm_hp.begin(), sortVec_pm_hp.end(), val_ph);
        auto it2 = std::find(sortVec_ppmm_hhpp.begin(), sortVec_ppmm_hhpp.end(), val_ph);
        if (it1 != sortVec_pm_hp.end() || it2 != sortVec_ppmm_hhpp.end()){
            for (int bA_ph=0; bA_ph<m_Nh*(m_Ns-m_Nh); bA_ph++){
                if (val_ph == blockArrays_pm_ph(0,bA_ph)){
                    range_upper = bA_ph+1;
                    counter += 1;
                }
            }
        }
        range_lower = range_upper - counter;
        counter = 0;
        indexHolder_pm_ph.col(ph) << range_lower, range_upper; //this now has same indices as sortVec
    }

    if (m_triplesOn){
        counter     = 0;
        range_lower = 0;
        range_upper = 0;
        indexHolder_pm_hh.conservativeResize(2,sortVec_pm_hh.size());
        for (int pp=0; pp<sortVec_pm_hh.size(); pp++){
            int val_hh = sortVec_pm_hh[pp];
            auto it1 = std::find(sortVec_pm_ph.begin(), sortVec_pm_ph.end(), val_hh);
            if (it1 != sortVec_pm_ph.end()){
                for (int bA_hh=0; bA_hh<m_Nh*m_Nh; bA_hh++){
                    if (val_hh == blockArrays_pm_hh(0,bA_hh)){
                        range_upper = bA_hh+1;
                        counter += 1;
                    }
                }
            }
            range_lower = range_upper - counter;
            counter = 0;
            indexHolder_pm_hh.col(pp) << range_lower, range_upper; //this now has same indices as sortVec
        }


        counter     = 0;
        range_lower = 0;
        range_upper = 0;
        indexHolder_pp_hp.conservativeResize(2,sortVec_pp_hp.size());
        for (int hp=0; hp<sortVec_pp_hp.size(); hp++){
            int val_hp = sortVec_pp_hp[hp];
            auto it1 = std::find(sortVec_pp_ph.begin(), sortVec_pp_ph.end(), val_hp);
            auto it2 = std::find(sortVec_pp_hh.begin(), sortVec_pp_hh.end(), val_hp);
            if (it1 != sortVec_pp_ph.end() || it2 != sortVec_pp_hh.end()){
                for (int bA_hp=0; bA_hp<m_Nh*(m_Ns-m_Nh); bA_hp++){
                    if (val_hp == blockArrays_pp_hp(0,bA_hp)){
                        range_upper = bA_hp+1;
                        counter += 1;
                    }
                }
            }
            range_lower = range_upper - counter;
            counter = 0;
            indexHolder_pp_hp.col(hp) << range_lower, range_upper; //this now has same indices as sortVec
        }


        counter     = 0;
        range_lower = 0;
        range_upper = 0;
        indexHolder_pp_ph.conservativeResize(2,sortVec_pp_ph.size());
        for (int ph=0; ph<sortVec_pp_ph.size(); ph++){
            int val_ph = sortVec_pp_ph[ph];
            auto it1 = std::find(sortVec_pp_hp.begin(), sortVec_pp_hp.end(), val_ph);
            auto it2 = std::find(sortVec_pp_pp.begin(), sortVec_pp_pp.end(), val_ph);
            if (it1 != sortVec_pp_hp.end() || it2 != sortVec_pp_pp.end()){
                for (int bA_ph=0; bA_ph<m_Nh*(m_Ns-m_Nh); bA_ph++){
                    if (val_ph == blockArrays_pp_ph(0,bA_ph)){
                        range_upper = bA_ph+1;
                        counter += 1;
                    }
                }
            }
            range_lower = range_upper - counter;
            counter = 0;
            indexHolder_pp_ph.col(ph) << range_lower, range_upper; //this now has same indices as sortVec
        }

        counter     = 0;
        range_lower = 0;
        range_upper = 0;
        indexHolder_pm_pp.conservativeResize(2,sortVec_pm_pp.size());
        for (int pp=0; pp<sortVec_pm_pp.size(); pp++){
            int val_pp = sortVec_pm_pp[pp];
            auto it1 = std::find(sortVec_pm_ph.begin(), sortVec_pm_ph.end(), val_pp);
            if (it1 != sortVec_pm_ph.end()){
                for (int bA_pp=0; bA_pp<(m_Ns-m_Nh)*(m_Ns-m_Nh); bA_pp++){
                    if (val_pp == blockArrays_pm_pp(0,bA_pp)){
                        range_upper = bA_pp+1;
                        counter += 1;
                    }
                }
            }
            range_lower = range_upper - counter;
            counter = 0;
            indexHolder_pm_pp.col(pp) << range_lower, range_upper; //this now has same indices as sortVec
        }
    }

    counter     = 0;
    range_lower = 0;
    range_upper = 0;
    indexHolder_hp_s.conservativeResize(2,sortVec_hp_s.size());
    for (int hp=0; hp<sortVec_hp_s.size(); hp++){
        int val_hp = sortVec_hp_s[hp];
        for (int bA_hp=range_upper; bA_hp<m_Nh*(m_Ns-m_Nh); bA_hp++){
            if (val_hp == blockArrays_hp_s(0,bA_hp)){
                range_upper = bA_hp+1;
                counter += 1;
            }
            else if (val_hp < blockArrays_hp_s(0,bA_hp)){
                break;
            }
        }
        range_lower = range_upper - counter;
        counter = 0;
        indexHolder_hp_s.col(hp) << range_lower, range_upper; //this now has same indices as sortVec
    }

    counter     = 0;
    range_lower = 0;
    range_upper = 0;
    indexHolder_pp_pp.conservativeResize(2,sortVec_pp_pp.size());
    for (int pp=0; pp<sortVec_pp_pp.size(); pp++){
        int val_pp = sortVec_pp_pp[pp];
        auto it1 = std::find(sortVec_pp_hh.begin(), sortVec_pp_hh.end(), val_pp);
        auto it2 = std::find(sortVec_pp_ph.begin(), sortVec_pp_ph.end(), val_pp);
        auto it3 = std::find(sortVec_pppm_hhhp.begin(), sortVec_pppm_hhhp.end(), val_pp);
        if (it1 != sortVec_pp_hh.end() || it2 != sortVec_pp_ph.end() || it3 != sortVec_pppm_hhhp.end()){
            for (int bA_pp=range_upper; bA_pp<(m_Ns-m_Nh)*(m_Ns-m_Nh); bA_pp++){
                if (val_pp == blockArrays_pp_pp(0,bA_pp)){
                    range_upper = bA_pp+1;
                    counter += 1;
                }
                else if(val_pp < blockArrays_pp_pp(0,bA_pp)){
                    break;
                }
            }
        }
        range_lower = range_upper - counter;
        counter = 0;
        indexHolder_pp_pp.col(pp) << range_lower, range_upper; //this now has same indices as sortVec
    }

    if (world_rank==0){cout << "made two-particle indexHolders" << endl;}

    counter     = 0;
    range_lower = 0;
    range_upper = 0;
    indexHolder_ppm_hhp.conservativeResize(2,sortVec_ppm_hhp.size());
    for (int hhp=0; hhp<sortVec_ppm_hhp.size(); hhp++){
        int val_hhp = sortVec_ppm_hhp[hhp];
        auto it1 = std::find(sortVec_p_p.begin(), sortVec_p_p.end(), val_hhp);
        auto it2 = std::find(sortVec_ppm_pph.begin(), sortVec_ppm_pph.end(), val_hhp);
        if (it1 != sortVec_p_p.end() || it2 != sortVec_ppm_pph.end()){
            for (int bA_hhp=range_upper; bA_hhp<m_Nh*m_Nh*(m_Ns-m_Nh); bA_hhp++){
                if (val_hhp == blockArrays_ppm_hhp(0,bA_hhp)){
                    range_upper = bA_hhp+1;
                    counter += 1;
                }
                else if(val_hhp < blockArrays_ppm_hhp(0,bA_hhp)){
                    break;
                }
            }
        }
        range_lower = range_upper - counter;
        counter = 0;
        indexHolder_ppm_hhp.col(hhp) << range_lower, range_upper; //this now has same indices as sortVec
    }

    counter     = 0;
    range_lower = 0;
    range_upper = 0;
    indexHolder_ppm_pph.conservativeResize(2,sortVec_ppm_pph.size());
    for (int pph=0; pph<sortVec_ppm_pph.size(); pph++){
        int val_pph = sortVec_ppm_pph[pph];
        auto it1 = std::find(sortVec_p_h.begin(), sortVec_p_h.end(), val_pph);
        auto it2 = std::find(sortVec_ppm_hhp.begin(), sortVec_ppm_hhp.end(), val_pph);
        if (it1 != sortVec_p_h.end() || it2 != sortVec_ppm_hhp.end()){
            for (int bA_pph=range_upper; bA_pph<m_Nh*(m_Ns-m_Nh)*(m_Ns-m_Nh); bA_pph++){
                if (val_pph == blockArrays_ppm_pph(0,bA_pph)){
                    range_upper = bA_pph+1;
                    counter += 1;
                }
                else if(val_pph < blockArrays_ppm_pph(0,bA_pph)){
                    break;
                }
            }
        }
        range_lower = range_upper - counter;
        counter = 0;
        indexHolder_ppm_pph.col(pph) << range_lower, range_upper; //this now has same indices as sortVec
    }

    if (m_triplesOn){
        counter     = 0;
        range_lower = 0;
        range_upper = 0;
        indexHolder_ppp_hhh.conservativeResize(2,sortVec_ppp_hhh.size());
        for (int hhh=0; hhh<sortVec_ppp_hhh.size(); hhh++){
            int val_hhh = sortVec_ppp_hhh[hhh];
            auto it = std::find(sortVec_ppp_ppp.begin(), sortVec_ppp_ppp.end(), val_hhh);
            if (it != sortVec_ppp_ppp.end()){
                for (int bA_hhh=range_upper; bA_hhh<m_Nh*m_Nh*m_Nh; bA_hhh++){
                    if (val_hhh == blockArrays_ppp_hhh(0,bA_hhh)){
                        range_upper = bA_hhh+1;
                        counter += 1;
                    }
                    else if (val_hhh < blockArrays_ppp_hhh(0,bA_hhh)){
                        break;
                    }
                }
            }
            range_lower = range_upper - counter;
            counter = 0;
            indexHolder_ppp_hhh.col(hhh) << range_lower, range_upper; //this now has same indices as sortVec
        }

        counter     = 0;
        range_lower = 0;
        range_upper = 0;
        indexHolder_ppm_hhh.conservativeResize(2,sortVec_ppm_hhh.size());
        for (int hhh=0; hhh<sortVec_ppm_hhh.size(); hhh++){
            int val_hhh = sortVec_ppm_hhh[hhh];
            auto it = std::find(sortVec_p_p.begin(), sortVec_p_p.end(), val_hhh);
            if (it != sortVec_p_p.end()){
                for (int bA_hhh=range_upper; bA_hhh<m_Nh*m_Nh*m_Nh; bA_hhh++){
                    if (val_hhh == blockArrays_ppm_hhh(0,bA_hhh)){
                        range_upper = bA_hhh+1;
                        counter += 1;
                    }
                    else if (val_hhh < blockArrays_ppm_hhh(0,bA_hhh)){
                        break;
                    }
                }
            }
            range_lower = range_upper - counter;
            counter = 0;
            indexHolder_ppm_hhh.col(hhh) << range_lower, range_upper; //this now has same indices as sortVec
        }

        counter     = 0;
        range_lower = 0;
        range_upper = 0;
        indexHolder_ppp_ppp.conservativeResize(2,sortVec_ppp_ppp.size());
        for (int ppp=0; ppp<sortVec_ppp_ppp.size(); ppp++){
            int val_ppp = sortVec_ppp_ppp[ppp];
            auto it = std::find(sortVec_ppp_hhh.begin(), sortVec_ppp_hhh.end(), val_ppp);
            if (it != sortVec_ppp_hhh.end()){
                for (int bA_ppp=range_upper; bA_ppp<blockArrays_ppp_ppp.cols(); bA_ppp++){
                    if (val_ppp == blockArrays_ppp_ppp(0,bA_ppp)){
                        range_upper = bA_ppp+1;
                        counter += 1;
                    }
                    else if (val_ppp < blockArrays_ppp_ppp(0,bA_ppp)){
                        break;
                    }
                }
            }
            range_lower = range_upper - counter;
            counter = 0;
            indexHolder_ppp_ppp.col(ppp) << range_lower, range_upper; //this now has same indices as sortVec
        }
    }

    if (world_rank==0){cout << "made three-particle indexHolders" << endl;}

    if (m_triplesOn){
        counter     = 0;
        range_lower = 0;
        range_upper = 0;
        indexHolder_pppm_hhhp.conservativeResize(2,sortVec_pppm_hhhp.size());
        for (int hhhp=0; hhhp<sortVec_pppm_hhhp.size(); hhhp++){
            int val_hhhp = sortVec_pppm_hhhp[hhhp];
            auto it = std::find(sortVec_pp_pp.begin(), sortVec_pp_pp.end(), val_hhhp);
            if (it != sortVec_pp_pp.end()){
                for (int bA_hhhp=range_upper; bA_hhhp<blockArrays_pppm_hhhp.cols(); bA_hhhp++){
                    if (val_hhhp == blockArrays_pppm_hhhp(0,bA_hhhp)){
                        range_upper = bA_hhhp+1;
                        counter += 1;
                    }
                    else if (val_hhhp < blockArrays_pppm_hhhp(0,bA_hhhp)){
                        break;
                    }
                }
            }
            range_lower = range_upper - counter;
            counter = 0;
            indexHolder_pppm_hhhp.col(hhhp) << range_lower, range_upper; //this now has same indices as sortVec
        }

        //this isn't used anymore after T5a fix
        /*counter     = 0;
        range_lower = 0;
        range_upper = 0;
        indexHolder_ppmm_hhpp.conservativeResize(2,sortVec_ppmm_hhpp.size());
        for (int hhpp=0; hhpp<sortVec_ppmm_hhpp.size(); hhpp++){
            int val_hhpp = sortVec_ppmm_hhpp[hhpp];
            auto it = std::find(sortVec_pm_ph.begin(), sortVec_pm_ph.end(), val_hhpp);
            if (it != sortVec_pm_ph.end()){
                for (int bA_hhpp=range_upper; bA_hhpp<blockArrays_ppmm_hhpp.cols(); bA_hhpp++){
                    if (val_hhpp == blockArrays_ppmm_hhpp(0,bA_hhpp)){
                        range_upper = bA_hhpp+1;
                        counter += 1;
                    }
                    else if (val_hhpp < blockArrays_ppmm_hhpp(0,bA_hhpp)){
                        break;
                    }
                }
            }
            range_lower = range_upper - counter;
            counter = 0;
            indexHolder_ppmm_hhpp.col(hhpp) << range_lower, range_upper; //this now has same indices as sortVec
        }*/

            counter     = 0;
            range_lower = 0;
            range_upper = 0;
            indexHolder_ppmm_pphh.conservativeResize(2,sortVec_ppmm_pphh.size());
            for (int pphh=0; pphh<sortVec_ppmm_pphh.size(); pphh++){
                int val_pphh = sortVec_ppmm_pphh[pphh];
                auto it = std::find(sortVec_pm_hp.begin(), sortVec_pm_hp.end(), val_pphh);
                if (it != sortVec_pm_hp.end()){
                    for (int bA_pphh=range_upper; bA_pphh<blockArrays_ppmm_pphh.cols(); bA_pphh++){
                        if (val_pphh == blockArrays_ppmm_pphh(0,bA_pphh)){
                            range_upper = bA_pphh+1;
                            counter += 1;
                        }
                        else if (val_pphh < blockArrays_ppmm_pphh(0,bA_pphh)){
                            break;
                        }
                    }
                }
                range_lower = range_upper - counter;
                counter = 0;
                indexHolder_ppmm_pphh.col(pphh) << range_lower, range_upper; //this now has same indices as sortVec
        }

        counter     = 0;
        range_lower = 0;
        range_upper = 0;
        indexHolder_pppm_ppph.conservativeResize(2,sortVec_pppm_ppph.size());
        for (int ppph=0; ppph<sortVec_pppm_ppph.size(); ppph++){
            int val_ppph = sortVec_pppm_ppph[ppph];
            auto it = std::find(sortVec_pp_hh.begin(), sortVec_pp_hh.end(), val_ppph);
            if (it != sortVec_pp_hh.end()){
                for (int bA_ppph=range_upper; bA_ppph<blockArrays_pppm_ppph.cols(); bA_ppph++){
                    if (val_ppph == blockArrays_pppm_ppph(0,bA_ppph)){
                        range_upper = bA_ppph+1;
                        counter += 1;
                    }
                    else if (val_ppph < blockArrays_pppm_ppph(0,bA_ppph)){
                        break;
                    }
                }
            }
            range_lower = range_upper - counter;
            counter = 0;
            indexHolder_pppm_ppph.col(ppph) << range_lower, range_upper; //this now has same indices as sortVec
        }
    }

    if (world_rank==0){cout << "made four-particle indexHolders" << endl;}

    if (m_triplesOn){
        counter     = 0;
        range_lower = 0;
        range_upper = 0;
        indexHolder_pppmm_hhhpp.conservativeResize(2,sortVec_pppmm_hhhpp.size());
        for (int hhhpp=0; hhhpp<sortVec_pppmm_hhhpp.size(); hhhpp++){
            int val_hhhpp = sortVec_pppmm_hhhpp[hhhpp];
            auto it = std::find(sortVec_p_p.begin(), sortVec_p_p.end(), val_hhhpp);
            if (it != sortVec_p_p.end()){
                for (int bA_hhhpp=range_upper; bA_hhhpp<blockArrays_pppmm_hhhpp.cols(); bA_hhhpp++){
                    if (val_hhhpp == blockArrays_pppmm_hhhpp(0,bA_hhhpp)){
                        range_upper = bA_hhhpp+1;
                        counter += 1;
                    }
                    else if (val_hhhpp < blockArrays_pppmm_hhhpp(0,bA_hhhpp)){
                        break;
                    }
                }
            }
            range_lower = range_upper - counter;
            counter = 0;
            indexHolder_pppmm_hhhpp.col(hhhpp) << range_lower, range_upper; //this now has same indices as sortVec
        }

        counter     = 0;
        range_lower = 0;
        range_upper = 0;
        indexHolder_pppmm_ppphh.conservativeResize(2,sortVec_pppmm_ppphh.size());
        for (int ppphh=0; ppphh<sortVec_pppmm_ppphh.size(); ppphh++){
            int val_ppphh = sortVec_pppmm_ppphh[ppphh];
            auto it = std::find(sortVec_p_h.begin(), sortVec_p_h.end(), val_ppphh);
            if (it != sortVec_p_h.end()){
                for (int bA_ppphh=range_upper; bA_ppphh<blockArrays_pppmm_ppphh.cols(); bA_ppphh++){
                    if (val_ppphh == blockArrays_pppmm_ppphh(0,bA_ppphh)){
                        range_upper = bA_ppphh+1;
                        counter += 1;
                    }
                    else if (val_ppphh < blockArrays_pppmm_ppphh(0,bA_ppphh)){
                        break;
                    }
                }
            }
            range_lower = range_upper - counter;
            counter = 0;
            indexHolder_pppmm_ppphh.col(ppphh) << range_lower, range_upper; //this now has same indices as sortVec
        }
    }

    if (world_rank==0){cout << "made five-particle indexHolders" << endl;}

    counter     = 0;
    range_lower = 0;
    range_upper = 0;
    for (int h=0; h<sortVec_pp_hh.size(); h++){
        for (int p=0; p<sortVec_pp_pp.size(); p++){
            int val_hh = sortVec_pp_hh[h];
            int val_pp = sortVec_pp_pp[p];
            if (val_hh == val_pp){      //ensures I only work on cases where hh and pp have equal kunique
                for (int hh=0; hh<m_Nh*m_Nh; hh++){
                    if ( val_hh == blockArrays_pp_hh(0,hh) ){
                        range_upper_hh = hh+1;
                        counter += 1;
                    }
                }
                range_lower_hh = range_upper_hh - counter;
                counter = 0;
                for (int pp=0; pp<(m_Ns-m_Nh)*(m_Ns-m_Nh); pp++){
                    if ( val_hh == blockArrays_pp_pp(0,pp) ){
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
    //cout << numOfKu <<" "<< sortVec_pp_hh.size() << endl;


    /*std::vector<int> temps;
    for (int i1=0; i1<sortVec_p_p.size(); i1++){
        for (int i2=0; i2<sortVec_ppm_pph.size(); i2++){
            for (int i3=0; i3<sortVec_ppm_hhp.size(); i3++){
                if ( sortVec_p_p[i1] == sortVec_ppm_pph[i2] && sortVec_p_p[i1] == sortVec_ppm_hhp[i3]){
                    if (std::find(temps.begin(), temps.end(), sortVec_p_p[i1]) == temps.end()){
                        temps.push_back( sortVec_p_p[i1] );
                    }
                }
            }
        }
    }
    for (int i1=0; i1<sortVec_p_h.size(); i1++){
        for (int i2=0; i2<sortVec_ppm_hhp.size(); i2++){
            for (int i3=0; i3<sortVec_ppm_pph.size(); i3++){
                if ( sortVec_p_h[i1] == sortVec_ppm_hhp[i2] && sortVec_p_h[i1] == sortVec_ppm_pph[i3]){
                    if (std::find(temps.begin(), temps.end(), sortVec_p_h[i1]) == temps.end()){
                        temps.push_back( sortVec_p_h[i1] );
                    }
                }
            }
        }
    }*/



    if (m_triplesOn){
        counter     = 0;
        range_lower = 0;
        range_upper = 0;
        for (int hhh=0; hhh<sortVec_ppp_hhh.size(); hhh++){
            for (int ppp=0; ppp<sortVec_ppp_ppp.size(); ppp++){
                int val_hhh = sortVec_ppp_hhh[hhh];
                int val_ppp = sortVec_ppp_ppp[ppp];
                if (val_hhh == val_ppp /*&& std::find(temps.begin(), temps.end(), val_hhh) != temps.end()*/ ){      //ensures I only work on cases where hhh and ppp have equal kunique
                    for (int hhh2=0; hhh2<blockArrays_ppp_hhh.cols(); hhh2++){
                        if ( val_hhh == blockArrays_ppp_hhh(0,hhh2) ){
                            range_upper_hhh = hhh2+1;
                            counter += 1;
                        }
                    }
                    range_lower_hhh = range_upper_hhh - counter;
                    counter = 0;
                    for (int ppp2=0; ppp2<blockArrays_ppp_ppp.cols(); ppp2++){
                        if ( val_hhh == blockArrays_ppp_ppp(0,ppp2) ){
                            range_upper_ppp = ppp2+1;
                            counter += 1;
                        }
                    }
                    range_lower_ppp = range_upper_ppp - counter;
                    counter = 0;
                    boundsHolder_hhhppp_hhh.conservativeResize(Eigen::NoChange, boundsHolder_hhhppp_hhh.cols()+1);
                    boundsHolder_hhhppp_ppp.conservativeResize(Eigen::NoChange, boundsHolder_hhhppp_ppp.cols()+1);
                    boundsHolder_hhhppp_hhh.col(boundsHolder_hhhppp_hhh.cols()-1) << range_lower_hhh, range_upper_hhh;
                    //cout << "test" << endl;
                    boundsHolder_hhhppp_ppp.col(boundsHolder_hhhppp_ppp.cols()-1) << range_lower_ppp, range_upper_ppp;
                    Vhhhppp_i.push_back(val_hhh);
                }
            }
        }
        numOfKu3 = boundsHolder_hhhppp_hhh.cols();
    }


    if (world_rank==0){cout << "made indexHolders" << endl;}

    //make Vhhpp map
    /*for (int h=0; h<boundsHolder_hhpp_hh.cols(); h++){
        range_lower_hh = boundsHolder_hhpp_hh(0,h);
        range_upper_hh = boundsHolder_hhpp_hh(1,h);
        range_lower_pp = boundsHolder_hhpp_pp(0,h);
        range_upper_pp = boundsHolder_hhpp_pp(1,h);
        makeMatMap_hhpp(blockArrays_pp_hh, blockArrays_pp_pp,range_lower_hh, range_upper_hh, range_lower_pp, range_upper_pp);
    }*/
    for (int hh=0; hh<sortVec_pp_hh.size(); hh++){
        for (int pp=0; pp<sortVec_pp_pp.size(); pp++){
            if (sortVec_pp_hh[hh] == sortVec_pp_pp[pp] ){
                range_lower_hh = indexHolder_pp_hh(0,hh);
                range_upper_hh = indexHolder_pp_hh(1,hh);
                range_lower_pp = indexHolder_pp_pp(0,pp);
                range_upper_pp = indexHolder_pp_pp(1,pp);
                makeMatMap_hhpp(blockArrays_pp_hh, blockArrays_pp_pp,range_lower_hh, range_upper_hh, range_lower_pp, range_upper_pp);
            }
        }
    }
    if (world_rank==0){cout << "made Vhhpp" << endl;}

    if (m_triplesOn){
        for (int pp=0; pp<sortVec_pp_pp.size(); pp++){
            for (int ph=0; ph<sortVec_pp_ph.size(); ph++){
                if (sortVec_pp_pp[pp] == sortVec_pp_ph[ph] ){
                    range_lower_pp = indexHolder_pp_pp(0,pp);
                    range_upper_pp = indexHolder_pp_pp(1,pp);
                    range_lower_ph = indexHolder_pp_ph(0,ph);
                    range_upper_ph = indexHolder_pp_ph(1,ph);
                    makeMatMap_ppph(blockArrays_pp_pp, blockArrays_pp_ph,range_lower_pp, range_upper_pp, range_lower_ph, range_upper_ph);
                }
            }
        }
        if (world_rank==0){cout << "made Vppph" << endl;}

        for (int hh=0; hh<sortVec_pp_hh.size(); hh++){
            for (int hp=0; hp<sortVec_pp_hp.size(); hp++){
                if (sortVec_pp_hh[hh] == sortVec_pp_hp[hp] ){
                    range_lower_hh = indexHolder_pp_hh(0,hh);
                    range_upper_hh = indexHolder_pp_hh(1,hh);
                    range_lower_hp = indexHolder_pp_hp(0,hp);
                    range_upper_hp = indexHolder_pp_hp(1,hp);
                    makeMatMap_hhhp(blockArrays_pp_hh, blockArrays_pp_hp,range_lower_hh, range_upper_hh, range_lower_hp, range_upper_hp);
                }
            }
        }
        if (world_rank==0){cout << "made Vhhhp" << endl;}
    }

    //Below we make Vpppp, Vhhhh, and Vhphp, which are distinct from Vhhpp --> don't need any remapping of elements
    for (int h=0; h<indexHolder_pp_hh.cols(); h++){
        range_lower_hh = indexHolder_pp_hh(0,h);
        range_upper_hh = indexHolder_pp_hh(1,h);
        Vhhhh.push_back( makeSquareBlock(blockArrays_pp_hh, range_lower_hh, range_upper_hh) );
    }

    //int prev = 0;
    //int length = boundsHolder_hhpp_hh.cols();

    for (int hp=0; hp<indexHolder_pm_hp.cols(); hp++){
        range_lower_hp = indexHolder_pm_hp(0,hp);
        range_upper_hp = indexHolder_pm_hp(1,hp);
        Vhphp.push_back( makeSquareBlock_s(blockArrays_pm_hp, range_lower_hp, range_upper_hp) );
    }

    for (int p=0; p<indexHolder_pp_pp.cols(); p++){//I could boundsHolder here, but I think T2c has more channels that in hh
        range_lower_pp = indexHolder_pp_pp(0,p);
        range_upper_pp = indexHolder_pp_pp(1,p);
        //Vpppp.push_back( makeRektBlock(blockArrays_pp,blockArrays_pp, range_lower_pp, range_upper_pp,range_lower_pp, range_upper_pp) );
        Vpppp.push_back( makeSquareBlock(blockArrays_pp_pp, range_lower_pp, range_upper_pp) );

        //this printout is handy for large systems
        /*if (prev < p/(double)length){
            cout << 100*p/(double)length << "%  " << range_upper_pp-range_lower_pp << endl;
            prev = p/(double)length;
        }*/
        //Vpppp_i.push_back( sortVec_pp_pp[p] );
    }


    if (world_rank==0){cout << "made Vpppp" << endl;}
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

void MakeIntMat::makeMatMap_ppph(Eigen::MatrixXi& array1, Eigen::MatrixXi& array2, int range_lower1, int range_upper1, int range_lower2, int range_upper2){
    for (int i = range_lower1; i<range_upper1; i++){
        for (int j = range_lower2; j<range_upper2; j++){
            double val = m_system->assym((array1)(1,i), (array1)(2,i), (array2)(1,j), (array2)(2,j));
            if (val != 0){
                Vppph_elements[Identity_ppph((array1)(1,i), (array1)(2,i), (array2)(1,j), (array2)(2,j))] = val;
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
    return h1 + h2*m_Nh + h3*m_Nh*m_Nh + p1*m_Nh*m_Nh*m_Nh;
}

int MakeIntMat::Identity_hhpp(int h1, int h2, int p1, int p2){
    return h1 + h2*m_Nh + p1*m_Nh*m_Nh + p2*m_Nh*m_Nh*(m_Ns);
}

int MakeIntMat::Identity_ppph(int p1, int p2, int p3, int h1){
    int out = h1 + p1*m_Nh + p2*m_Nh*(m_Ns) + p3*m_Nh*(m_Ns)*(m_Ns);
    if(out<0){std::cout << "id_ppph " << out << std::endl;}
    //std::cout << out << std::endl;
    return out;
    //return h1 + p1*m_Nh + p2*m_Nh*(m_Ns-m_Nh) + p3*m_Nh*(m_Ns-m_Nh)*(m_Ns-m_Nh);
}

//there are 3 ways to calculate "out" here, but when running the program, I found no difference between them
unsigned long int MakeIntMat::Identity_hhhppp(int h1, int h2, int h3, int p1, int p2, int p3){
    /*unsigned long int out  = (unsigned long int)h1
                                + (unsigned long int)h2*m_Nh
                                + (unsigned long int)h3*m_Nh*m_Nh
                                + (unsigned long int)p1*m_Nh*m_Nh*m_Nh
                                + (unsigned long int)p2*m_Nh*m_Nh*m_Nh*(m_Ns)
                                + (unsigned long int)p3*m_Nh*m_Nh*m_Nh*(m_Ns)*(m_Ns);*/

    unsigned long int out  = (unsigned long int)h1
                                + (unsigned long int)h2*m_Nh
                                + (unsigned long int)h3*m_Nh2
                                + (unsigned long int)p1*m_Nh3
                                + (unsigned long int)p2*m_Nh3Ns
                                + (unsigned long int)p3*m_Nh3Ns2;

    /*unsigned long int out  = (unsigned long int)h1
                                + m_Nh*((unsigned long int)h2
                                + m_Nh*((unsigned long int)h3
                                + m_Nh*((unsigned long int)p1
                                + m_Ns*((unsigned long int)p2
                                + m_Ns*((unsigned long int)p3)))));*/

    return out;
}

unsigned long int MakeIntMat::Identity_hhhhhp(int h1, int h2, int h3, int h4, int h5, int p1){
    unsigned long int out = (unsigned long int)h1
                            + (unsigned long int)h2*m_Nh
                            + (unsigned long int)h3*m_Nh*m_Nh
                            + (unsigned long int)h4*m_Nh*m_Nh*m_Nh
                            + (unsigned long int)h5*m_Nh*m_Nh*m_Nh*m_Nh
                            + (unsigned long int)p1*m_Nh*m_Nh*m_Nh*m_Nh*m_Nh;
    //if(out<0){std::cout << "id_hhhhhp " << out << std::endl;}
    return out;
    //return h1 + h2*m_Nh + h3*m_Nh*m_Nh + h4*m_Nh*m_Nh*m_Nh + h5*m_Nh*m_Nh*m_Nh*m_Nh + p1*m_Nh*m_Nh*m_Nh*m_Nh*m_Nh;
}