#include "diagrams.h"

#include <eigen3/Eigen/Dense>

Diagrams::Diagrams()
{
}

void Diagrams::setSystem(class System* system){
    m_system = system;
}

/*The makeIxJBlock(ku,i1,i2,i3,i4) functions work as follows:
 * ku:      this is the channel with which we need to perform the diagram product
 * i1, i2:  these are the indices of the rows
 * i3, i4:  these are the indices of the columns
 */

void Diagrams::setIntClass(class MakeIntMat* intClass){
    m_intClass = intClass;
}

void Diagrams::setAmpClass(class MakeAmpMat* ampClass){
    m_ampClass = ampClass;
}

void Diagrams::La(){
    for (int n=0; n<m_intClass->Vhhpp_i.size(); n++){
        int ku = m_intClass->Vhhpp_i[n];
        Eigen::MatrixXd mat1 = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T_elements);    // t_ij^cd
        Eigen::MatrixXd mat2 = m_intClass->Vpppp[n];                            // v_ab^cd
        Eigen::MatrixXd product = 0.5*mat1*mat2;                  // (t_ij^cd)(t_cd^ab)

        m_ampClass->make2x2Block_inverse(product,ku,0,0,1,1, m_ampClass->T_elements_new, true);
    }
}

void Diagrams::Lb(){
    for (int n=0; n<m_intClass->Vhhpp_i.size(); n++){
        int ku = m_intClass->Vhhpp_i[n];
        Eigen::MatrixXd mat1 = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T_elements);    // t_kl^ab
        Eigen::MatrixXd mat2 = m_intClass->Vhhhh[n];                            // v_ij^kl
        Eigen::MatrixXd product = 0.5*mat2*mat1;                  // (v_ij^kl)(t_kl^ab)

        m_ampClass->make2x2Block_inverse(product,ku,0,0,1,1, m_ampClass->T_elements_new, true);
    }
}

void Diagrams::Lc(){

    //I made the script below to try and fix alignment
    //but after a meeting i've realised i may have gone about this the wrong way
    //which means i musn't use std::map, rather a vector

    /*int range_lower_hh = m_intClass->boundsHolder_hhpp_hh(0,index);
    int range_upper_hh = m_intClass->boundsHolder_hhpp_hh(1,index);
    int range_lower_pp = m_intClass->boundsHolder_hhpp_pp(0,index);
    int range_upper_pp = m_intClass->boundsHolder_hhpp_pp(1,index);

    int dim_hh = range_upper_hh - range_lower_hh;
    int dim_pp = range_upper_pp - range_lower_pp;

    Eigen::MatrixXd insertMat;
    insertMat.conservativeResize(dim_hh, dim_pp);

    int ku1; //a-j
    int ku2; //i-b
    int ku3; //c-k

    for (int h = range_lower_hh; h<range_upper_hh; h++){
        int i = m_intClass->blockArrays_hh(0,h);
        int j = m_intClass->blockArrays_hh(1,h);
        for (int g = range_lower_pp; g<range_upper_pp; g++){
            int a = m_intClass->blockArrays_pp(0,h);
            int b = m_intClass->blockArrays_pp(1,h);
            ku1 = m_system->kUnique2(a,j,1,-1);
            ku2 = m_system->kUnique2(i,b,1,-1);

            if (ku1 == ku2){ //if these are the same, we can start searching for c-k terms
                auto it = std::find(m_intClass->sortVec_ph.begin(), m_intClass->sortVec_ph.end(), ku1);
                if (it != m_intClass->sortVec_ph.end()){
                    int range_lower = m_intClass->indexHolder_ph(0,it);
                    int range_upper = m_intClass->indexHolder_ph(1,it);
                    for (int f = range_lower; f<range_upper; f++){
                        int c = m_intClass->indexHolder_ph(0,f);
                        int k = m_intClass->indexHolder_ph(1,f);
                        ku3 = m_system->kUnique2(c,k,1,-1);
                        if (ku1 == ku3){

                        }
                    }
                }
            }
        }
    }*/

    for (int i1=0; i1<m_intClass->sortVec_hp.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ph.size(); i2++){
            if ( m_intClass->sortVec_hp[i1] == m_intClass->sortVec_ph[i2] ){   // && it2 != sortVec2.end()){
                int ku = m_intClass->sortVec_hp[i1];
                Eigen::MatrixXd mat1 = m_ampClass->make2x2Block(ku,0,1,1,0, m_ampClass->T_elements);
                Eigen::MatrixXd mat2 = m_intClass->Vhphp[i1];
                //std::cout << mat1 << std::endl;
                Eigen::MatrixXd product= -mat1*mat2;
                m_ampClass->make2x2Block_inverse(product, ku, 0,1,1,0, m_ampClass->T_elements_new, true);
            }
        }
    }
}

void Diagrams::Qa(){
    for (int n=0; n<m_intClass->numOfKu; n++){
        int ku = m_intClass->Vhhpp_i[n];
        Eigen::MatrixXd mat1 = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T_elements);    // t_ij^ab
        Eigen::MatrixXd mat2 = m_intClass->make2x2Block(ku,0,0,1,1);    // v_ij^ab
        Eigen::MatrixXd mat3 = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T_elements);    // t_ij^ab
        Eigen::MatrixXd product = 0.25*mat1*mat2.transpose()*mat3;                  // (t_ij^cd)(v_cd^kl)(t_kl^ab)

        m_ampClass->make2x2Block_inverse(product,ku,0,0,1,1, m_ampClass->T_elements_new, true);
    }
    //return 0.25*product;
}

void Diagrams::Qb(){
    for (int i1=0; i1<m_intClass->sortVec_hp.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ph.size(); i2++){
            if ( m_intClass->sortVec_hp[i1] == m_intClass->sortVec_ph[i2] ){   // && it2 != sortVec2.end()){
                int ku = m_intClass->sortVec_hp[i1];
                Eigen::MatrixXd mat1 = m_ampClass->make2x2Block(ku,0,1,1,0, m_ampClass->T_elements);
                Eigen::MatrixXd mat2 = m_intClass->make2x2Block(ku,0,1,1,0);
                Eigen::MatrixXd mat3 = mat1;
                //std::cout << mat1 << std::endl;
                Eigen::MatrixXd product= 0.5*mat1*mat2*mat3.transpose();
                m_ampClass->make2x2Block_inverse(product, ku, 0,1,1,0, m_ampClass->T_elements_new, true);
            }
        }
    }
}

void Diagrams::Qc(){
    for (int i1=0; i1<m_intClass->sortVec_p.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_hhp.size(); i2++){
            if ( m_intClass->sortVec_p[i1] == m_intClass->sortVec_hhp[i2] ){
                int ku = m_intClass->sortVec_p[i1];
                Eigen::MatrixXd mat1 = m_ampClass->make3x1Block(ku,0,0,1,1, m_ampClass->T_elements);    // t_ijb^d
                Eigen::MatrixXd mat2 = m_intClass->make3x1Block(ku,0,0,1,1);                            // v_klc^d
                Eigen::MatrixXd mat3 = mat1;    // t_klc^a
                Eigen::MatrixXd product = -0.5*mat1*mat2.transpose()*mat3;                  // (t_ijb^d)(v_d^klc)(t_klc^a)
                m_ampClass->make3x1Block_inverse(product, ku, 0,0,1,1, m_ampClass->T_elements_new, true);
            }
        }
    }
}

void Diagrams::Qd(){
    for (int i1=0; i1<m_intClass->sortVec_h.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pph.size(); i2++){
            if ( m_intClass->sortVec_h[i1] == m_intClass->sortVec_pph[i2] ){
                int ku = m_intClass->sortVec_h[i1];
                Eigen::MatrixXd mat1 = m_ampClass->make3x1Block(ku,1,1,0,0, m_ampClass->T_elements);    // t_ijb^d
                Eigen::MatrixXd mat2 = m_intClass->make3x1Block(ku,1,1,0,0);                            // v_klc^d
                Eigen::MatrixXd mat3 = mat1;    // t_klc^a
                Eigen::MatrixXd product = -0.5*mat1*mat2.transpose()*mat3;                  // (t_ijb^d)(v_d^klc)(t_klc^a)
                //std::cout << product << std::endl;
                m_ampClass->make3x1Block_inverse(product, ku, 1,1,0,0, m_ampClass->T_elements_new, true);
            }
        }
    }
}

Eigen::MatrixXd Pab(Eigen::MatrixXd inMat, int ku, int i1, int i2, int i3, int i4){
    int rowDim = inMat.rows();
    int colDim = inMat.cols();

    Eigen::MatrixXd outMat;
    outMat.conservativeResize(rowDim, colDim);

    bool Mat3x1 = ( rowDim == 3 && colDim==1 );
    bool Mat2x2 = ( rowDim == 2 && colDim==2 );

    if ( Mat3x1 ){
        bool cond_hhp = (i1 == 0 && i2 == 0 && i3==1);
        bool cond_pph = (i1 == 1 && i2 == 1 && i3==0);
        bool cond_h   = (i4 == 0);
        bool cond_p   = (i4 == 1);

    }
    else{
        bool cond_hh1 = (i1 == 0 && i2 == 0);
        bool cond_hp1 = (i1 == 0 && i2 == 1);
        bool cond_ph1 = (i1 == 1 && i2 == 0);
        bool cond_pp1 = (i1 == 1 && i2 == 1);

        bool cond_hh2 = (i3 == 0 && i4 == 0);
        bool cond_hp2 = (i3 == 0 && i4 == 1);
        bool cond_ph2 = (i3 == 1 && i4 == 0);
        bool cond_pp2 = (i3 == 1 && i4 == 1);
    }
}
