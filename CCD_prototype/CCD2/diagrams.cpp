#include "diagrams.h"

#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Core"
#include <omp.h>
#include <mpi.h>
#include <chrono>
#include <time.h>
#include <unistd.h>

typedef std::chrono::high_resolution_clock Clock;   //needed for timing

Diagrams::Diagrams()
{
}

void Diagrams::setSystem(class System* system){
    m_system = system;
}

/*The make2x2Block(ku,i1,i2,i3,i4) functions work as follows:
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


// ##################################################
// ##                                              ##
// ## DOUBLES DIAGRAMS                             ##
// ##                                              ##
// ##################################################

void Diagrams::La(){
    for (int i1=0; i1<m_intClass->sortVec_pp_hh.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pp_pp.size(); i2++){
            if ( m_intClass->sortVec_pp_hh[i1]==m_intClass->sortVec_pp_pp[i2] ){

                int ku = m_intClass->sortVec_pp_hh[i1];

                /*int m = 0; int index;
            while (m<m_intClass->sortVec_pp_pp.size()){
                if (m_intClass->sortVec_pp_pp[m] == ku){
                    index = m;
                    m = m_intClass->sortVec_pp_pp.size();
                }
                m++;
            }*/

                //auto it = std::find(sortVec_pp_pp.begin(), sortVec_pp_pp.end(), ku);

                Eigen::MatrixXd mat1 = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T2_elements);            // t_ij^cd
                Eigen::MatrixXd mat2 = m_intClass->Vpppp[i2];                                                    // v_ab^cd
                Eigen::MatrixXd product = 0.5*(mat1*mat2);                                                        // (t_ij^cd)(t_cd^ab)

                m_ampClass->make2x2Block_inverse(product,ku,0,0,1,1, m_ampClass->T2_elements_new, true);
            }
        }
    }
    m_ampClass->addElementsT2(0,0);
    m_ampClass->T2_temp.clear();
}

void Diagrams::Lb(){
    for (int i1=0; i1<m_intClass->sortVec_pp_hh.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pp_pp.size(); i2++){
            if ( m_intClass->sortVec_pp_hh[i1]==m_intClass->sortVec_pp_pp[i2] ){

                int ku = m_intClass->sortVec_pp_hh[i1];

                /*int m = 0; int index;
        while (m<m_intClass->sortVec_pp_hh.size()){
            if (m_intClass->sortVec_pp_hh[m] == ku){
                index = m;
                m = m_intClass->sortVec_pp_hh.size();
            }
            m++;
        }*/

                Eigen::MatrixXd mat1 = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T2_elements);            // t_kl^ab
                Eigen::MatrixXd mat2 = m_intClass->Vhhhh[i1];                                                    // v_ij^kl
                Eigen::MatrixXd product = 0.5*(mat2*mat1);                                                        // (v_ij^kl)(t_kl^ab)

                m_ampClass->make2x2Block_inverse(product,ku,0,0,1,1, m_ampClass->T2_elements_new, true);
            }
        }
    }
    m_ampClass->addElementsT2(0,0);
    m_ampClass->T2_temp.clear();
}

void Diagrams::Lc(){
    for (int i1=0; i1<m_intClass->sortVec_pm_hp.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pm_ph.size(); i2++){
            if ( m_intClass->sortVec_pm_hp[i1] == m_intClass->sortVec_pm_ph[i2] ){
                int ku = m_intClass->sortVec_pm_ph[i2];

                Eigen::MatrixXd mat1 = m_intClass->Vhphp[i1];
                Eigen::MatrixXd mat2 = m_ampClass->make2x2Block(ku,0,1,1,0, m_ampClass->T2_elements);
                Eigen::MatrixXd product= -1*(mat1*mat2);

                m_ampClass->make2x2Block_inverse(product, ku, 0,1,1,0, m_ampClass->T2_elements_new, true);
            }
        }
    }
    m_ampClass->addElementsT2(1,1);
    m_ampClass->T2_temp.clear();
}

void Diagrams::Qa(){
    for (int n=0; n<m_intClass->numOfKu; n++){
        int ku = m_intClass->Vhhpp_i[n];

        Eigen::MatrixXd mat1 = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T2_elements);            // t_ij^ab
        Eigen::MatrixXd mat2 = m_intClass->make2x2Block(ku,0,0,1,1);                                    // v_ij^ab
        Eigen::MatrixXd mat3 = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T2_elements);            // t_ij^ab
        Eigen::MatrixXd product = 0.25*(mat1*mat2.transpose()*mat3);                                      // (t_ij^cd)(v_cd^kl)(t_kl^ab)

        m_ampClass->make2x2Block_inverse(product,ku,0,0,1,1, m_ampClass->T2_elements_new, true);
    }
    m_ampClass->addElementsT2(0,0);
    m_ampClass->T2_temp.clear();
}

void Diagrams::Qb(){
    for (int i1=0; i1<m_intClass->sortVec_pm_hp.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pm_ph.size(); i2++){
            if ( m_intClass->sortVec_pm_hp[i1] == m_intClass->sortVec_pm_ph[i2] ){
                int ku = m_intClass->sortVec_pm_hp[i1];

                Eigen::MatrixXd mat1 = m_ampClass->make2x2Block(ku,0,1,1,0, m_ampClass->T2_elements);
                Eigen::MatrixXd mat2 = m_intClass->make2x2Block(ku,0,1,1,0);
                Eigen::MatrixXd mat3 = mat1;
                Eigen::MatrixXd product= 0.5*(mat1*mat2.transpose()*mat3);

                m_ampClass->make2x2Block_inverse(product, ku, 0,1,1,0, m_ampClass->T2_elements_new, true);
            }
        }
    }
    m_ampClass->addElementsT2(1,1);
    m_ampClass->T2_temp.clear();
}

void Diagrams::Qc(){
    for (int i1=0; i1<m_intClass->sortVec_p_p.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_hhp.size(); i2++){
            if ( m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_hhp[i2] ){
                int ku = m_intClass->sortVec_p_p[i1];

                Eigen::MatrixXd mat1 = m_ampClass->make3x1Block(ku,0,0,1,1, m_ampClass->T2_elements);    // t_ijb^d
                Eigen::MatrixXd mat2 = m_intClass->make3x1Block(ku,0,0,1,1);                            // v_klc^d
                Eigen::MatrixXd mat3 = mat1;                                                            // t_klc^a
                Eigen::MatrixXd product = -0.5*(mat1*mat2.transpose()*mat3);                              // (t_ijb^d)(v_d^klc)(t_klc^a)

                m_ampClass->make3x1Block_inverse(product, ku, 0,0,1,1, m_ampClass->T2_elements_new, true);
            }
        }
    }
    m_ampClass->addElementsT2(0,1);
    m_ampClass->T2_temp.clear();
}

void Diagrams::Qd(){
    for (int i1=0; i1<m_intClass->sortVec_p_h.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_pph.size(); i2++){
            if ( m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_pph[i2] ){
                int ku = m_intClass->sortVec_p_h[i1];
                Eigen::MatrixXd mat1 = m_ampClass->make3x1Block(ku,1,1,0,0, m_ampClass->T2_elements);    // t_ijb^d
                Eigen::MatrixXd mat2 = m_intClass->make3x1Block(ku,1,1,0,0);                            // v_klc^d
                Eigen::MatrixXd mat3 = mat1;                                                            // t_klc^a
                Eigen::MatrixXd product = -0.5*(mat1*mat2.transpose()*mat3);                              // (t_ijb^d)(v_d^klc)(t_klc^a)

                m_ampClass->make3x1Block_inverse(product, ku, 1,1,0,0, m_ampClass->T2_elements_new, true);
            }
        }
    }
    m_ampClass->addElementsT2(1,0);
    m_ampClass->T2_temp.clear();
}

// ##################################################
// ##                                              ##
// ## INTERMEDIATES                                ##
// ##                                              ##
// ##################################################

void Diagrams::I1_term1(){
    for (int i1=0; i1<m_intClass->sortVec_pp_hh.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pp_pp.size(); i2++){
            if ( m_intClass->sortVec_pp_hh[i1] == m_intClass->sortVec_pp_pp[i2] ){

                Eigen::MatrixXd mat1 = m_ampClass->I1_makemat_1(i1,i2);
                Eigen::MatrixXd mat2 = m_intClass->I1_makemat(i1,i2);

                Eigen::MatrixXd I1   = 0.5*(mat1*mat2.transpose());

                Eigen::MatrixXd mat3 = m_intClass->Vhhhh[i1];

                I1.noalias() += mat3;

                Eigen::MatrixXd mat4 = m_ampClass->I1_makemat_2(i1,i2);

                Eigen::MatrixXd product = 0.5*(I1*mat4);

                m_ampClass->I1_inverse(product, i1, i2);
            }
        }
    }
    m_ampClass->addElementsT2(0,0);
    m_ampClass->T2_temp.clear();
}

void Diagrams::I2_term1(){
    for (int i1=0; i1<m_intClass->sortVec_pm_hp.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pm_ph.size(); i2++){
            if ( m_intClass->sortVec_pm_hp[i1] == m_intClass->sortVec_pm_ph[i2] ){

                Eigen::MatrixXd mat1 = m_ampClass->I2_makemat_1(i1,i2);
                Eigen::MatrixXd mat2 = m_intClass->I2_makemat(i1,i2);

                Eigen::MatrixXd I2   = 0.5*(mat1*mat2.transpose());

                Eigen::MatrixXd mat3 = m_intClass->Vhphp[i1];

                I2.noalias() += mat3; //the minus is intentional

                Eigen::MatrixXd mat4 = m_ampClass->I2_makemat_2(i1,i2);

                Eigen::MatrixXd product = I2*mat4;

                m_ampClass->I2_inverse(product, i1, i2);
            }
        }
    }
    m_ampClass->addElementsT2(1,1);
    m_ampClass->T2_temp.clear();
}

void Diagrams::I3_term1(){
    for (int i1=0; i1<m_intClass->sortVec_p_h.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_pph.size(); i2++){
            if ( m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_pph[i2] ){
                Eigen::MatrixXd mat1 = m_ampClass->I3_makemat_1(i2,i1);
                Eigen::MatrixXd mat2 = m_intClass->I3_makemat(i2,i1);
                //std::cout << "sup" << std::endl;
                Eigen::MatrixXd mat3 = m_ampClass->I3_makemat_2(i2,i1);

                Eigen::MatrixXd product   = -0.5*(mat3.transpose()*mat2*mat1.transpose());

                m_ampClass->I3_inverse(product.transpose(), i2, i1);
            }
        }
    }
    m_ampClass->addElementsT2(1,0);
    m_ampClass->T2_temp.clear();
}

void Diagrams::I4_term1(){
    for (int i1=0; i1<m_intClass->sortVec_p_p.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_hhp.size(); i2++){
            if ( m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_hhp[i2] ){

                Eigen::MatrixXd mat1 = m_ampClass->I4_makemat_1(i2,i1);
                Eigen::MatrixXd mat2 = m_intClass->I4_makemat(i2,i1);
                Eigen::MatrixXd mat3 = m_ampClass->I4_makemat_2(i2,i1);

                Eigen::MatrixXd product   = -0.5*(mat3.transpose()*mat2*mat1.transpose());

                m_ampClass->I4_inverse(product.transpose(), i2, i1);
            }
        }
    }
    m_ampClass->addElementsT2(0,1);
    m_ampClass->T2_temp.clear();
}

void Diagrams::I1_term(){
    for (int n=0; n<m_intClass->Vhhpp_i.size(); n++){
        int ku = m_intClass->Vhhpp_i[n];

        Eigen::MatrixXd mat1 = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T2_elements);
        Eigen::MatrixXd mat2 = m_intClass->make2x2Block(ku,0,0,1,1);

        Eigen::MatrixXd I1   = 0.5*mat1*mat2.transpose();

        Eigen::MatrixXd mat3 = m_intClass->Vhhhh[n].transpose();

        I1 += mat3;

        Eigen::MatrixXd product = 0.5*I1*mat1;

        //std::cout << product << std::endl;

        m_ampClass->make2x2Block_inverse(product,ku,0,0,1,1, m_ampClass->T2_elements_new, true);
    }
    m_ampClass->addElementsT2(0,0);
    m_ampClass->T2_temp.clear();
}

void Diagrams::I2_term(){
    for (int i1=0; i1<m_intClass->sortVec_pm_hp.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pm_ph.size(); i2++){
            if ( m_intClass->sortVec_pm_hp[i1] == m_intClass->sortVec_pm_ph[i2] ){
                int ku = m_intClass->sortVec_pm_ph[i2];


                Eigen::MatrixXd mat1 = m_ampClass->make2x2Block(ku,0,1,1,0, m_ampClass->T2_elements);
                Eigen::MatrixXd mat2 = m_intClass->make2x2Block(ku,0,1,1,0);

                Eigen::MatrixXd I2   = 0.5*mat2*mat1.transpose();

                Eigen::MatrixXd mat3 = m_intClass->Vhphp[i1];

                I2 -= mat3;

                Eigen::MatrixXd product = I2.transpose()*mat1;

                m_ampClass->make2x2Block_inverse(product, ku, 0,1,1,0, m_ampClass->T2_elements_new, true);
            }
        }
    }
    m_ampClass->addElementsT2(1,1);
    m_ampClass->T2_temp.clear();
}

void Diagrams::I3_term(){
    for (int i1=0; i1<m_intClass->sortVec_p_h.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_pph.size(); i2++){
            if ( m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_pph[i2] ){
                int ku = m_intClass->sortVec_p_h[i1];

                Eigen::MatrixXd mat1 = m_ampClass->make3x1Block(ku,1,1,0,0, m_ampClass->T2_elements);
                Eigen::MatrixXd mat2 = m_intClass->make3x1Block(ku,1,1,0,0);

                Eigen::MatrixXd I3   = mat2.transpose()*mat1;

                Eigen::MatrixXd product = -0.5*mat1*I3;

                m_ampClass->make3x1Block_inverse(product, ku, 1,1,0,0, m_ampClass->T2_elements_new, true);
            }
        }
    }
    m_ampClass->addElementsT2(1,0);
    m_ampClass->T2_temp.clear();
}

void Diagrams::I4_term(){
    for (int i1=0; i1<m_intClass->sortVec_p_p.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_hhp.size(); i2++){
            if ( m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_hhp[i2] ){
                int ku = m_intClass->sortVec_p_p[i1];

                Eigen::MatrixXd mat1 = m_ampClass->make3x1Block(ku,0,0,1,1, m_ampClass->T2_elements);
                Eigen::MatrixXd mat2 = m_intClass->make3x1Block(ku,0,0,1,1);

                Eigen::MatrixXd I3   = mat2.transpose()*mat1;

                Eigen::MatrixXd product = -0.5*mat1*I3;

                m_ampClass->make3x1Block_inverse(product, ku, 0,0,1,1, m_ampClass->T2_elements_new, true);
            }
        }
    }
    m_ampClass->addElementsT2(0,1);
    m_ampClass->T2_temp.clear();
}

// ##################################################
// ##                                              ##
// ## T3 contributions to T2                       ##
// ##                                              ##
// ##################################################

void Diagrams::D10b(){
    for (int i1=0; i1<m_intClass->sortVec_p_p.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_pph.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_ppm_hhp.size(); i3++){
            if ( m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_pph[i2] && m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_hhp[i3]){

                //int ku = m_intClass->sortVec_p_p[i1];
                //Eigen::MatrixXd mat1 = m_ampClass->make3x3Block(ku,0,0,1,1,1,0, m_ampClass->T3_elements);    // t_ijb^d
                //Eigen::MatrixXd mat1 = m_ampClass->make3x3Block_I(ku,0,0,1,1,1,0, m_ampClass->T3_elements_A);    // t_ijb^d
                //Eigen::MatrixXd mat2 = m_intClass->make3x1Block(ku,1,1,0,1);                            // v_klc^d
                //m_ampClass->make3x1Block_inverse_D10b(product, ku, 0,0,1,1, m_ampClass->T2_elements_new, true);

                Eigen::MatrixXd mat1 = m_ampClass->D10b_makemat(i3, i2);
                Eigen::MatrixXd mat2 = m_intClass->D10b_makemat(i2, i1);
                Eigen::MatrixXd product = -0.5*(mat1*mat2);                              // (t_ijb^d)(v_d^klc)(t_klc^a)

                m_ampClass->D10b_inverse(product, i3, i1);
            }
            }
        }
    }
    m_ampClass->addElementsT2(0,1);
    m_ampClass->T2_temp.clear();
}

void Diagrams::D10c(){
    for (int i1=0; i1<m_intClass->sortVec_p_h.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_hhp.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_ppm_pph.size(); i3++){
            if ( m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_hhp[i2] && m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_pph[i3]){

                //int ku = m_intClass->sortVec_p_h[i1];
                //Eigen::MatrixXd mat1 = m_ampClass->make3x3Block(ku,0,0,1,1,1,0, m_ampClass->T3_elements);    // t_ijb^d
                //Eigen::MatrixXd mat1 = m_ampClass->make3x3Block_I_D10c(ku,0,0,1,1,1,0, m_ampClass->T3_elements_A);    // t_ijb^d
                //Eigen::MatrixXd mat2 = m_intClass->make3x1Block(ku,0,0,1,0);                            // v_klc^d
                //m_ampClass->make3x1Block_inverse(product, ku, 1,1,0,0, m_ampClass->T2_elements_new, true);

                Eigen::MatrixXd mat1 = m_ampClass->D10c_makemat(i2, i3);
                Eigen::MatrixXd mat2 = m_intClass->D10c_makemat(i2, i1);
                Eigen::MatrixXd product = -0.5*(mat1.transpose()*mat2);                              // (t_ijb^d)(v_d^klc)(t_klc^a)

                m_ampClass->D10c_inverse(product, i3, i1);
            }
            }
        }
    }
    m_ampClass->addElementsT2(1,0);
    m_ampClass->T2_temp.clear();
}



// ##################################################
// ##                                              ##
// ## TRIPLES DIAGRAMS                             ##
// ##                                              ##
// ##################################################

void Diagrams::makeT3(){
    int index = 0; unsigned long int id;
    int i; int j; int k;
    int a; int b; int c;
    int n = 4;//omp_get_max_threads();
    std::cout << n << std::endl;
    m_ampClass->T3_elements_IV.resize(n);
    for (int channel = 0; channel<m_intClass->numOfKu3; channel++){
        //ku = m_intClass->Vhhhppp_i[channel];

        int range_lower1 = m_intClass->boundsHolder_hhhppp_hhh(0,channel);
        int range_upper1 = m_intClass->boundsHolder_hhhppp_hhh(1,channel);
        int range_lower2 = m_intClass->boundsHolder_hhhppp_ppp(0,channel);
        int range_upper2 = m_intClass->boundsHolder_hhhppp_ppp(1,channel);

        for (int hhh=range_lower1; hhh<range_upper1; hhh++){
            for (int ppp=range_lower2; ppp<range_upper2; ppp++){

                i = (m_intClass->blockArrays_ppp_hhh)(1,hhh);
                j = (m_intClass->blockArrays_ppp_hhh)(2,hhh);
                k = (m_intClass->blockArrays_ppp_hhh)(3,hhh);
                a = (m_intClass->blockArrays_ppp_ppp)(1,ppp);
                b = (m_intClass->blockArrays_ppp_ppp)(2,ppp);
                c = (m_intClass->blockArrays_ppp_ppp)(3,ppp);
                id = m_intClass->Identity_hhhppp(i,j,k,a,b,c);

                //if (m_ampClass->T3_elements_I.find(id) != m_ampClass->T3_elements_I.end()){std::cout << id << std::endl;}
                m_ampClass->T3_elements_I[id] = index;

                for (int j=0; j<n; j++){
                    m_ampClass->T3_elements_IV[j][id] = index;
                }

                index ++;
            }
        }
    }

/*    m_ampClass->T3_elements_IV.resize(omp_get_num_threads());

    for (int j=0; j<omp_get_num_threads(); j++){
        m_ampClass->T3_elements_IV[j] = m_ampClass->T3_elements_I;
    }*/

    m_ampClass->T3_elements_A.resize(m_ampClass->T3_elements_I.size(), 0);
    m_ampClass->T3_elements_A_new.resize(m_ampClass->T3_elements_I.size(), 0);
    m_ampClass->T3_elements_A_temp.resize(m_ampClass->T3_elements_I.size(), 0);

    //std::cout << m_ampClass->T3_elements_A.size() << std::endl;

    m_ampClass->T3_makeDirectMat();
}

Eigen::MatrixXi Diagrams::distributeChannels(Eigen::MatrixXi channels, int elements){

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Status status;

    //broadcast the necessesary blockArrays
    /*std::cout << m_intClass->blockArrays_p_h.cols()*m_intClass->blockArrays_p_h.rows() << std::endl;
    std::cout << m_intClass->blockArrays_ppm_pph.cols()*m_intClass->blockArrays_ppm_pph.rows()<< std::endl;
    MPI_Bcast(m_intClass->blockArrays_p_h.data(), m_intClass->blockArrays_p_h.cols()*m_intClass->blockArrays_p_h.rows(), MPI_INT, 0, MPI_COMM_WORLD);
    std::cout << "sup1" << std::endl;
    MPI_Bcast(m_intClass->blockArrays_ppm_hhp.data(), m_intClass->blockArrays_ppm_hhp.cols()*m_intClass->blockArrays_ppm_hhp.rows(), MPI_INT, 0, MPI_COMM_WORLD);
    std::cout << "sup2" << std::endl;
    MPI_Bcast(m_intClass->blockArrays_ppm_pph.data(), m_intClass->blockArrays_ppm_pph.cols()*m_intClass->blockArrays_ppm_pph.rows(), MPI_INT, 0, MPI_COMM_WORLD);
    */

    Eigen::MatrixXi matches_root;
    Eigen::MatrixXi matches_recv;

    if (world_rank == 0){ matches_root = channels;}

    int delegated_channels;
    int delegated_columns;
    int remain;
    int displs[world_size];
    int sendCount[world_size];

    if (world_rank == 0){
        remain              = (matches_root.cols() % world_size)*matches_root.rows();
        delegated_columns   = matches_root.cols()/world_size;
        delegated_channels  = delegated_columns*matches_root.rows();
    }

    MPI_Bcast(&delegated_channels, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&delegated_columns,  1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&remain,  1, MPI_INT, 0, MPI_COMM_WORLD);

    for (int rank=0; rank<world_size; rank++){
        displs[rank] = delegated_channels*rank;
        sendCount[rank] = delegated_channels;
        if (rank == world_size-1 && remain!=0){ sendCount[rank] += remain; }
    }

    if (world_rank!=world_size-1){
        matches_recv.conservativeResize(elements, delegated_columns);
    }
    else if (world_rank==world_size-1){
        matches_recv.conservativeResize(elements, delegated_columns + remain/elements);
    }


    MPI_Scatterv(matches_root.data(), sendCount, displs, MPI_INT, matches_recv.data(), sendCount[world_rank], MPI_INT, 0, MPI_COMM_WORLD);
    //MPI_Scatter(matches_root.data(), delegated_channels, MPI_INT, matches_recv.data(), delegated_channels, MPI_INT, 0, MPI_COMM_WORLD);

    return matches_recv;
}

void Diagrams::T1a(){

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Status status;


    Eigen::MatrixXi matches_root;
    Eigen::MatrixXi matches_recv;

    if (world_rank == 0){
        for (int i1=0; i1<m_intClass->sortVec_p_p.size(); i1++){
            for (int i2=0; i2<m_intClass->sortVec_ppm_pph.size(); i2++){
                for (int i3=0; i3<m_intClass->sortVec_ppm_hhp.size(); i3++){
                    if ( m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_pph[i2] && m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_hhp[i3]){
                        matches_root.conservativeResize(3, matches_root.cols()+1);
                        matches_root.col(matches_root.cols()-1) << i1,i2,i3;
                    }
                }
            }
        }
    }

    matches_recv = distributeChannels(matches_root, 3);

    int i1; int i2; int i3;
    for (int i=0; i<matches_recv.cols(); i++){
        i1 = matches_recv(0,i); i2 = matches_recv(1,i); i3 = matches_recv(2,i);

        Eigen::MatrixXd mat1 = m_intClass->T1a_makemat(i2, i1);
        Eigen::MatrixXd mat2 = m_ampClass->T1a_makemat(i3, i1);
        Eigen::MatrixXd product = mat2*mat1.transpose();

        //std::cout << product.cols() << " " << product.rows() << std::endl;

        m_ampClass->T1a_inverse(product, i3, i2);
    }

    std::vector<double> TempVec;
    if (world_rank == 0){TempVec.resize(m_ampClass->T3_elements_A_temp.size()); }
    MPI_Reduce(m_ampClass->T3_elements_A_temp.data(), TempVec.data(), m_ampClass->T3_elements_A_temp.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (world_rank == 0){ m_ampClass->T3_elements_A_temp = TempVec ;}

    /*for (int i1=0; i1<m_intClass->sortVec_p_p.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_pph.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_ppm_hhp.size(); i3++){
                if ( m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_pph[i2] && m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_hhp[i3]){

                    Eigen::MatrixXd mat1 = m_intClass->T1a_makemat(i2, i1);
                    Eigen::MatrixXd mat2 = m_ampClass->T1a_makemat(i3, i1);
                    Eigen::MatrixXd product = mat2*mat1.transpose();

                    m_ampClass->T1a_inverse(product, i3, i2);
                }
            }
        }
    }*/

    if (world_rank == 0){
        m_ampClass->addElementsT3_T1a();
    }
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}


void Diagrams::T1b(){

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Status status;


    Eigen::MatrixXi matches_root;
    Eigen::MatrixXi matches_recv;

    if (world_rank == 0){
        for (int i1=0; i1<m_intClass->sortVec_p_h.size(); i1++){
            for (int i2=0; i2<m_intClass->sortVec_ppm_hhp.size(); i2++){
                for (int i3=0; i3<m_intClass->sortVec_ppm_pph.size(); i3++){
                    if ( m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_hhp[i2] && m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_pph[i3]){
                        matches_root.conservativeResize(3, matches_root.cols()+1);
                        matches_root.col(matches_root.cols()-1) << i1,i2,i3;
                    }
                }
            }
        }
    }

    matches_recv = distributeChannels(matches_root, 3);

    int i1; int i2; int i3;
    for (int i=0; i<matches_recv.cols(); i++){
        i1 = matches_recv(0,i); i2 = matches_recv(1,i); i3 = matches_recv(2,i);

        Eigen::MatrixXd mat1 = m_intClass->T1b_makemat(i2, i1);
        Eigen::MatrixXd mat2 = m_ampClass->T1b_makemat(i3, i1);
        Eigen::MatrixXd product = -1*(mat1*mat2.transpose());

        m_ampClass->T1b_inverse(product, i2, i3);
    }

    std::vector<double> TempVec;
    if (world_rank == 0){TempVec.resize(m_ampClass->T3_elements_A_temp.size()); }
    MPI_Reduce(m_ampClass->T3_elements_A_temp.data(), TempVec.data(), m_ampClass->T3_elements_A_temp.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (world_rank == 0){ m_ampClass->T3_elements_A_temp = TempVec ;}

    /*for (int i1=0; i1<m_intClass->sortVec_p_h.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_hhp.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_ppm_pph.size(); i3++){
                if ( m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_hhp[i2] && m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_pph[i3]){

                    Eigen::MatrixXd mat1 = m_intClass->T1b_makemat(i2, i1);
                    Eigen::MatrixXd mat2 = m_ampClass->T1b_makemat(i3, i1);
                    Eigen::MatrixXd product = -1*(mat1*mat2.transpose());

                    m_ampClass->T1b_inverse(product, i2, i3);
                }
            }
        }
    }*/

    if (world_rank == 0){
        m_ampClass->addElementsT3_T1b();
    }
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}

void Diagrams::T2c(){

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Status status;


    Eigen::MatrixXi matches_root;
    Eigen::MatrixXi matches_recv;

    if (world_rank == 0){
        for (int i1=0; i1<m_intClass->sortVec_pp_pp.size(); i1++){
            for (int i2=0; i2<m_intClass->sortVec_pppm_hhhp.size(); i2++){
                if ( m_intClass->sortVec_pp_pp[i1] == m_intClass->sortVec_pppm_hhhp[i2]){
                    matches_root.conservativeResize(2, matches_root.cols()+1);
                    matches_root.col(matches_root.cols()-1) << i1,i2;
                }
            }
        }
    }

    matches_recv = distributeChannels(matches_root, 2);

    int i1; int i2;
    for (int i=0; i<matches_recv.cols(); i++){
        i1 = matches_recv(0,i); i2 = matches_recv(1,i);

        Eigen::MatrixXd mat1 = m_intClass->Vpppp[i1];
        Eigen::MatrixXd mat2 = m_ampClass->T2c_makemat(i2, i1);
        Eigen::MatrixXd product = 0.5*mat2*mat1; //mathematically I need to transpose mat1, but it's symmetric

        m_ampClass->T2c_inverse(product, i2, i1);
    }

    std::vector<double> TempVec;
    if (world_rank == 0){TempVec.resize(m_ampClass->T3_elements_A_temp.size()); }
    MPI_Reduce(m_ampClass->T3_elements_A_temp.data(), TempVec.data(), m_ampClass->T3_elements_A_temp.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (world_rank == 0){ m_ampClass->T3_elements_A_temp = TempVec ;}

    if (world_rank == 0){
        m_ampClass->addElementsT3_T2c();
    }

    /*for (int i1=0; i1<m_intClass->sortVec_pp_pp.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pppm_hhhp.size(); i2++){
            if ( m_intClass->sortVec_pp_pp[i1] == m_intClass->sortVec_pppm_hhhp[i2]){
                Eigen::MatrixXd mat1 = m_intClass->Vpppp[i1];
                Eigen::MatrixXd mat2 = m_ampClass->T2c_makemat(i2, i1);
                Eigen::MatrixXd product = 0.5*mat2*mat1; //mathematically I need to transpose mat1, but it's symmetric

                m_ampClass->T2c_inverse(product, i2, i1);
            }
        }
    }

    m_ampClass->addElementsT3_T2c();*/
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}

void Diagrams::T2d(){

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Status status;


    Eigen::MatrixXi matches_root;
    Eigen::MatrixXi matches_recv;

    if (world_rank == 0){
        for (int i1=0; i1<m_intClass->sortVec_pp_hh.size(); i1++){
                for (int i2=0; i2<m_intClass->sortVec_pppm_ppph.size(); i2++){
                    if ( m_intClass->sortVec_pp_hh[i1] == m_intClass->sortVec_pppm_ppph[i2]){
                    matches_root.conservativeResize(2, matches_root.cols()+1);
                    matches_root.col(matches_root.cols()-1) << i1,i2;
                }
            }
        }
    }

    matches_recv = distributeChannels(matches_root, 2);

    int i1; int i2;
    for (int i=0; i<matches_recv.cols(); i++){
        i1 = matches_recv(0,i); i2 = matches_recv(1,i);

        Eigen::MatrixXd mat1 = m_intClass->Vhhhh[i1];
        Eigen::MatrixXd mat2 = m_ampClass->T2d_makemat(i2, i1);
        Eigen::MatrixXd product = 0.5*mat2*mat1; //mathematically I need to transpose mat1, but it's symmetric

        m_ampClass->T2d_inverse(product, i2, i1);
    }

    std::vector<double> TempVec;
    if (world_rank == 0){TempVec.resize(m_ampClass->T3_elements_A_temp.size()); }
    MPI_Reduce(m_ampClass->T3_elements_A_temp.data(), TempVec.data(), m_ampClass->T3_elements_A_temp.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (world_rank == 0){ m_ampClass->T3_elements_A_temp = TempVec ;}

    if (world_rank == 0){
        m_ampClass->addElementsT3_T2d();
    }

    /*for (int i1=0; i1<m_intClass->sortVec_pp_hh.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pppm_ppph.size(); i2++){
            if ( m_intClass->sortVec_pp_hh[i1] == m_intClass->sortVec_pppm_ppph[i2]){
                Eigen::MatrixXd mat1 = m_intClass->Vhhhh[i1];
                Eigen::MatrixXd mat2 = m_ampClass->T2d_makemat(i2, i1);
                Eigen::MatrixXd product = 0.5*mat2*mat1; //mathematically I need to transpose mat1, but it's symmetric

                m_ampClass->T2d_inverse(product, i2, i1);
            }
        }
    }

    m_ampClass->addElementsT3_T2d();*/
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}

void Diagrams::T2e(){

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Status status;


    Eigen::MatrixXi matches_root;
    Eigen::MatrixXi matches_recv;

    if (world_rank == 0){
        for (int i1=0; i1<m_intClass->sortVec_pm_hp.size(); i1++){
                for (int i2=0; i2<m_intClass->sortVec_ppmm_hhpp.size(); i2++){  //THIS IS WRONG (17/03/17)
                    if ( m_intClass->sortVec_pm_hp[i1] == m_intClass->sortVec_ppmm_hhpp[i2]){
                    matches_root.conservativeResize(2, matches_root.cols()+1);
                    matches_root.col(matches_root.cols()-1) << i1,i2;
                }
            }
        }
    }

    matches_recv = distributeChannels(matches_root, 2);

    int i1; int i2;
    for (int i=0; i<matches_recv.cols(); i++){
        i1 = matches_recv(0,i); i2 = matches_recv(1,i);

        Eigen::MatrixXd mat1 = m_intClass->Vhphp[i1]; //I think hphp was made with sign index +- on rows and columns
        Eigen::MatrixXd mat2 = m_ampClass->T2e_makemat(i2, i1);
        Eigen::MatrixXd product = mat2*mat1; //mathematically I need to transpose mat1, but it's symmetric

        m_ampClass->T2e_inverse(product, i2, i1);
    }

    std::vector<double> TempVec;
    if (world_rank == 0){TempVec.resize(m_ampClass->T3_elements_A_temp.size()); }
    MPI_Reduce(m_ampClass->T3_elements_A_temp.data(), TempVec.data(), m_ampClass->T3_elements_A_temp.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (world_rank == 0){ m_ampClass->T3_elements_A_temp = TempVec ;}

    if (world_rank == 0){
        m_ampClass->addElementsT3_T2e();
    }

    /*for (int i1=0; i1<m_intClass->sortVec_pm_hp.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppmm_hhpp.size(); i2++){  //THIS IS WRONG (17/03/17)
            if ( m_intClass->sortVec_pm_hp[i1] == m_intClass->sortVec_ppmm_hhpp[i2]){
                Eigen::MatrixXd mat1 = m_intClass->Vhphp[i1]; //I think hphp was made with sign index +- on rows and columns
                Eigen::MatrixXd mat2 = m_ampClass->T2e_makemat(i2, i1);
                Eigen::MatrixXd product = mat2*mat1; //mathematically I need to transpose mat1, but it's symmetric

                m_ampClass->T2e_inverse(product, i2, i1);
            }
        }
    }

    //m_ampClass->addElementsT3_T2e();*/
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}

void Diagrams::T3b(){

    for (int i1=0; i1<m_intClass->sortVec_pm_hp.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pm_ph.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_pm_pp.size(); i3++){
                if ( m_intClass->sortVec_pm_hp[i1] == m_intClass->sortVec_pm_ph[i2] && m_intClass->sortVec_pm_hp[i1] == m_intClass->sortVec_pm_pp[i3]){

                    Eigen::MatrixXd mat1 = m_ampClass->T3b_makemat_1(i1, i2);
                    Eigen::MatrixXd mat2 = m_intClass->T3b_makemat(i3, i1);
                    Eigen::MatrixXd product = mat2*mat1;

                    m_ampClass->T3b_Inverse_temp(product, i3, i2);
                }
            }
        }
    }

    //now use remapped function

    for (int i1=0; i1<m_intClass->sortVec_ppm_pph.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_p_p.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_ppm_hhp.size(); i3++){
                if ( m_intClass->sortVec_ppm_pph[i1] == m_intClass->sortVec_p_p[i2] && m_intClass->sortVec_ppm_pph[i1] == m_intClass->sortVec_ppm_hhp[i3]){

                    Eigen::MatrixXd mat1 = m_ampClass->T3b_makemat_2(i1, i2);
                    Eigen::MatrixXd mat2 = m_ampClass->T3b_makemat_3(i3, i2);
                    Eigen::MatrixXd product = -mat2*mat1.transpose();

                    m_ampClass->T3b_inverse(product, i3, i1);
                }
            }
        }
    }

    m_ampClass->addElementsT3_T3b();
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
    m_ampClass->T3D_remap.clear();

}

void Diagrams::T3c(){

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Status status;

    if (world_rank==1 || world_size<2){
        //std::cout << m_intClass->sortVec_pm_hp.size() << std::endl;
        for (int i1=0; i1<m_intClass->sortVec_pm_hp.size(); i1++){
            for (int i2=0; i2<m_intClass->sortVec_pm_ph.size(); i2++){
                for (int i3=0; i3<m_intClass->sortVec_pm_hh.size(); i3++){
                    if ( m_intClass->sortVec_pm_hp[i1] == m_intClass->sortVec_pm_ph[i2] && m_intClass->sortVec_pm_hp[i1] == m_intClass->sortVec_pm_hh[i3]){

                        Eigen::MatrixXd mat1 = m_ampClass->T3c_makemat_1(i1, i2);
                        Eigen::MatrixXd mat2 = m_intClass->T3c_makemat(i3, i2);
                        Eigen::MatrixXd product = mat1*mat2.transpose();

                        //std::cout << product << std::endl;
                        m_ampClass->T3c_Inverse_temp(product.transpose(), i3, i1);
                    }
                }
            }
        }

        //now use remapped function

        for (int i1=0; i1<m_intClass->sortVec_ppm_pph.size(); i1++){
            for (int i2=0; i2<m_intClass->sortVec_p_h.size(); i2++){
                for (int i3=0; i3<m_intClass->sortVec_ppm_hhp.size(); i3++){
                    if ( m_intClass->sortVec_ppm_pph[i1] == m_intClass->sortVec_p_h[i2] && m_intClass->sortVec_ppm_pph[i1] == m_intClass->sortVec_ppm_hhp[i3]){

                        Eigen::MatrixXd mat1 = m_ampClass->T3c_makemat_2(i3, i2);
                        Eigen::MatrixXd mat2 = m_ampClass->T3c_makemat_3(i1, i2);
                        Eigen::MatrixXd product = mat1*mat2.transpose();

                        m_ampClass->T3c_inverse(product, i3, i1);
                    }
                }
            }
        }
    }

    if (world_size>=2){
        if (world_rank == 1){MPI_Send(m_ampClass->T3_elements_A_temp.data(), m_ampClass->T3_elements_A_temp.size(), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);}
        if (world_rank == 0){MPI_Recv(m_ampClass->T3_elements_A_temp.data(), m_ampClass->T3_elements_A_temp.size(), MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);}

        if (world_rank == 0){
            m_ampClass->addElementsT3_T3c();
        }
    }
    else{
        m_ampClass->addElementsT3_T3c();
    }

    /*m_ampClass->addElementsT3_T3c();*/
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
    m_ampClass->T3D_remap.clear();

}

void Diagrams::T3d(){

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Status status;

    if (world_rank==2 || world_size < 3){
        for (int i1=0; i1<m_intClass->sortVec_pp_hh.size(); i1++){
            for (int i2=0; i2<m_intClass->sortVec_pp_ph.size(); i2++){
                for (int i3=0; i3<m_intClass->sortVec_pp_pp.size(); i3++){
                    if ( m_intClass->sortVec_pp_hh[i1] == m_intClass->sortVec_pp_ph[i2] && m_intClass->sortVec_pp_hh[i1] == m_intClass->sortVec_pp_pp[i3]){

                        Eigen::MatrixXd mat1 = m_ampClass->T3d_makemat_1(i1, i3);
                        Eigen::MatrixXd mat2 = m_intClass->T3d_makemat(i3, i2);
                        Eigen::MatrixXd product = mat1*mat2;

                        m_ampClass->T3d_Inverse_temp(product, i1, i2);
                    }
                }
            }
        }

        //now use remapped function
        //std::cout << m_intClass->blockArrays_ppp_ppp.size() << std::endl;

        for (int i1=0; i1<m_intClass->sortVec_ppm_pph.size(); i1++){
            for (int i2=0; i2<m_intClass->sortVec_p_h.size(); i2++){
                for (int i3=0; i3<m_intClass->sortVec_ppm_hhp.size(); i3++){
                    if ( m_intClass->sortVec_ppm_pph[i1] == m_intClass->sortVec_p_h[i2] && m_intClass->sortVec_ppm_pph[i1] == m_intClass->sortVec_ppm_hhp[i3]){

                        Eigen::MatrixXd mat1 = m_ampClass->T3d_makemat_2(i3, i2);
                        Eigen::MatrixXd mat2 = m_ampClass->T3d_makemat_3(i1, i2);
                        Eigen::MatrixXd product = 0.5*mat1*mat2.transpose();

                        m_ampClass->T3d_inverse(product, i3, i1);
                    }
                }
            }
        }
    }

    if (world_size >= 3 ){
        if (world_rank == 2){MPI_Send(m_ampClass->T3_elements_A_temp.data(), m_ampClass->T3_elements_A_temp.size(), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);}
        if (world_rank == 0){MPI_Recv(m_ampClass->T3_elements_A_temp.data(), m_ampClass->T3_elements_A_temp.size(), MPI_DOUBLE, 2, 1, MPI_COMM_WORLD, &status);}

        if (world_rank == 0){
            m_ampClass->addElementsT3_T3d();
        }
    }
    else{
        m_ampClass->addElementsT3_T3d();
    }

    /*m_ampClass->addElementsT3_T3d();*/
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
    m_ampClass->T3D_remap.clear();

}

void Diagrams::T3e(){

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Status status;

    if (world_rank==3 || world_size<4){
        for (int i1=0; i1<m_intClass->sortVec_ppm_hhp.size(); i1++){
            for (int i2=0; i2<m_intClass->sortVec_p_p.size(); i2++){
                for (int i3=0; i3<m_intClass->sortVec_ppm_hhh.size(); i3++){
                    if ( m_intClass->sortVec_ppm_hhp[i1] == m_intClass->sortVec_p_p[i2] && m_intClass->sortVec_ppm_hhp[i1] == m_intClass->sortVec_ppm_hhh[i3]){

                        Eigen::MatrixXd mat1 = m_ampClass->T3e_makemat_1(i1, i2);
                        Eigen::MatrixXd mat2 = m_intClass->T3e_makemat(i3, i2);
                        Eigen::MatrixXd product = mat1*mat2.transpose();

                        //std::cout << product << std::endl;
                        m_ampClass->T3e_Inverse_temp(product, i1, i3);
                    }
                }
            }
        }

        //now use remapped function

        for (int i1=0; i1<m_intClass->sortVec_pp_pp.size(); i1++){
            for (int i2=0; i2<m_intClass->sortVec_pp_hh.size(); i2++){
                for (int i3=0; i3<m_intClass->sortVec_pppm_hhhp.size(); i3++){
                    if ( m_intClass->sortVec_pp_pp[i1] == m_intClass->sortVec_pp_hh[i2] && m_intClass->sortVec_pp_pp[i1] == m_intClass->sortVec_pppm_hhhp[i3]){

                        Eigen::MatrixXd mat1 = m_ampClass->T3e_makemat_2(i3, i2);
                        Eigen::MatrixXd mat2 = m_ampClass->T3e_makemat_3(i2, i1);
                        Eigen::MatrixXd product = -0.5*mat1*mat2;

                        m_ampClass->T3e_inverse(product, i3, i1);
                    }
                }
            }
        }
    }

    if (world_size >= 4){
        if (world_rank == 3){MPI_Send(m_ampClass->T3_elements_A_temp.data(), m_ampClass->T3_elements_A_temp.size(), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);}
        if (world_rank == 0){MPI_Recv(m_ampClass->T3_elements_A_temp.data(), m_ampClass->T3_elements_A_temp.size(), MPI_DOUBLE, 3, 1, MPI_COMM_WORLD, &status);}

        if (world_rank == 0){
            m_ampClass->addElementsT3_T3e();
        }
    }
    else{
        m_ampClass->addElementsT3_T3e();
    }

    /*m_ampClass->addElementsT3_T3e();*/
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
    m_ampClass->T3D_remap.clear();

}

void Diagrams::T5a(){

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Status status;


    Eigen::MatrixXi matches_root;
    Eigen::MatrixXi matches_recv;

    if (world_rank == 0){
        for (int i1=0; i1<m_intClass->sortVec_pm_hp.size(); i1++){
                for (int i2=0; i2<m_intClass->sortVec_pm_ph.size(); i2++){
                    for (int i3=0; i3<m_intClass->sortVec_ppmm_hhpp.size(); i3++){
                        if ( m_intClass->sortVec_pm_hp[i1] == m_intClass->sortVec_ppmm_hhpp[i3] && m_intClass->sortVec_pm_hp[i1] == m_intClass->sortVec_pm_ph[i2]){
                        matches_root.conservativeResize(3, matches_root.cols()+1);
                        matches_root.col(matches_root.cols()-1) << i1,i2,i3;
                    }
                }
            }
        }
    }

    matches_recv = distributeChannels(matches_root, 3);

    int i1; int i2; int i3;
    for (int i=0; i<matches_recv.cols(); i++){
        i1 = matches_recv(0,i); i2 = matches_recv(1,i); i3 = matches_recv(2,i);

        Eigen::MatrixXd mat1 = m_ampClass->T5a_makemat_1(i1, i2);
        Eigen::MatrixXd mat2 = m_intClass->T5a_makemat(i1, i2);
        Eigen::MatrixXd mat3 = m_ampClass->T5a_makemat_2(i3, i2);
        Eigen::MatrixXd product = mat3*mat2.transpose()*mat1;

        m_ampClass->T5a_inverse(product, i3, i2);
    }

    std::vector<double> TempVec;
    if (world_rank == 0){TempVec.resize(m_ampClass->T3_elements_A_temp.size()); }
    MPI_Reduce(m_ampClass->T3_elements_A_temp.data(), TempVec.data(), m_ampClass->T3_elements_A_temp.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (world_rank == 0){ m_ampClass->T3_elements_A_temp = TempVec ;}

    if (world_rank == 0){
        m_ampClass->addElementsT3_T5a();
    }

    /*for (int i1=0; i1<m_intClass->sortVec_pm_hp.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pm_ph.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_ppmm_hhpp.size(); i3++){
                if ( m_intClass->sortVec_pm_hp[i1] == m_intClass->sortVec_ppmm_hhpp[i3] && m_intClass->sortVec_pm_hp[i1] == m_intClass->sortVec_pm_ph[i2]){

                    Eigen::MatrixXd mat1 = m_ampClass->T5a_makemat_1(i1, i2);
                    Eigen::MatrixXd mat2 = m_intClass->T5a_makemat(i1, i2);
                    Eigen::MatrixXd mat3 = m_ampClass->T5a_makemat_2(i3, i2);
                    Eigen::MatrixXd product = mat3*mat2.transpose()*mat1;

                    m_ampClass->T5a_inverse(product, i3, i2);
                }
            }
        }
    }

    m_ampClass->addElementsT3_T5a();*/
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}

void Diagrams::makeT5bIndexMat(){
    for (int i1=0; i1<m_intClass->sortVec_p_h.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_pph.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_pppmm_ppphh.size(); i3++){
                if ( m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_pppmm_ppphh[i3] && m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_pph[i2]){
                    Eigen::MatrixXi mat = m_ampClass->T5b_makemat_2_I( i3,i1 );
                    m_ampClass->T3_T5b_indices.push_back( mat );
                }
            }
        }
    }
}

void Diagrams::makeT5cIndexMat(){
    for (int i1=0; i1<m_intClass->sortVec_p_p.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_hhp.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_pppmm_hhhpp.size(); i3++){
                if ( m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_pppmm_hhhpp[i3] && m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_hhp[i2]){
                    Eigen::MatrixXi mat = m_ampClass->T5c_makemat_2_I( i3,i1 );
                    m_ampClass->T3_T5c_indices.push_back( mat );
                }
            }
        }
    }
}

Eigen::MatrixXd Diagrams::contructMatT5b(int index){
    int rows = m_ampClass->T3_T5b_indices[index].rows();
    int cols = m_ampClass->T3_T5b_indices[index].cols();
    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(rows, cols);
    #pragma omp parallel for /*num_threads(2)*/
    for (int row=0; row<rows; row++){
        for (int col=0; col<cols; col++){
            returnMat(row,col) = m_ampClass->T3_elements_A[m_ampClass->T3_T5b_indices[index](row, col) ];
        }
    }
    return returnMat;
}

Eigen::MatrixXd Diagrams::contructMatT5c(int index){
    int rows = m_ampClass->T3_T5c_indices[index].rows();
    int cols = m_ampClass->T3_T5c_indices[index].cols();
    Eigen::MatrixXd returnMat;
    returnMat.conservativeResize(rows, cols);
    #pragma omp parallel for /*num_threads(2)*/
    for (int row=0; row<rows; row++){
        for (int col=0; col<cols; col++){
            returnMat(row,col) = m_ampClass->T3_elements_A[m_ampClass->T3_T5c_indices[index](row, col) ];
        }
    }
    return returnMat;
}

void Diagrams::destroy5Map(){
    m_intClass->blockArrays_pppmm_hhhpp.conservativeResize(0,0);
    m_intClass->indexHolder_pppmm_hhhpp.conservativeResize(0,0);
    m_intClass->blockArrays_pppmm_ppphh.conservativeResize(0,0);
    m_intClass->indexHolder_pppmm_ppphh.conservativeResize(0,0);
}

void Diagrams::T5b(){

    /*int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Status status;


    Eigen::MatrixXi matches_root;
    Eigen::MatrixXi matches_recv;

    if (world_rank == 0){
        for (int i1=0; i1<m_intClass->sortVec_p_h.size(); i1++){
                for (int i2=0; i2<m_intClass->sortVec_ppm_pph.size(); i2++){
                    for (int i3=0; i3<m_intClass->sortVec_pppmm_ppphh.size(); i3++){
                        if ( m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_pppmm_ppphh[i3] && m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_pph[i2]){
                        matches_root.conservativeResize(3, matches_root.cols()+1);
                        matches_root.col(matches_root.cols()-1) << i1,i2,i3;
                    }
                }
            }
        }
    }

    matches_recv = distributeChannels(matches_root);

    int i1; int i2; int i3;
    for (int i=0; i<matches_recv.cols(); i++){
        i1 = matches_recv(0,i); i2 = matches_recv(1,i); i3 = matches_recv(2,i);

        Eigen::MatrixXd mat1 = m_ampClass->T5b_makemat_1(i2, i1);
        Eigen::MatrixXd mat2 = m_intClass->T5b_makemat(i2, i1);
        Eigen::MatrixXd mat3 = m_ampClass->T5b_makemat_2(i3, i1);
        Eigen::MatrixXd product = -0.5*mat3*mat2.transpose()*mat1;

        //std::cout << mat3.cols() << " " << mat3.rows() << " " << mat3.cols()*mat3.rows() << std::endl;
        m_ampClass->T5b_inverse(product, i3, i1);
    }

    std::vector<double> TempVec;
    if (world_rank == 0){TempVec.resize(m_ampClass->T3_elements_A_temp.size()); }
    MPI_Reduce(m_ampClass->T3_elements_A_temp.data(), TempVec.data(), m_ampClass->T3_elements_A_temp.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (world_rank == 0){ m_ampClass->T3_elements_A_temp = TempVec ;}

    if (world_rank == 0){
        m_ampClass->addElementsT3_T5b();
    }*/

    int index = 0;

    Eigen::MatrixXd mat1, mat2, mat3, product, tempMat;
    for (int i1=0; i1<m_intClass->sortVec_p_h.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_pph.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_pppmm_ppphh.size(); i3++){
                if ( m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_pppmm_ppphh[i3] && m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_pph[i2]){

                    mat1 = m_ampClass->T5b_makemat_1(i2, i1);
                    mat2 = m_intClass->T5b_makemat(i2, i1);
                    //Eigen::MatrixXd mat3 = m_ampClass->T5b_makemat_2(i3, i1);
                    tempMat = mat2.transpose()*mat1;
                    mat3 = contructMatT5b(index);
                    product = -0.5*mat3*tempMat;


                    //m_ampClass->T5b_inverse(product, i3, i1);
                    m_ampClass->T5b_inverse_I(product, index);
                    index++;
                }
            }
        }
    }

    m_ampClass->addElementsT3_T5b();
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}

void Diagrams::T5c(){

    /*int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Status status;


    Eigen::MatrixXi matches_root;
    Eigen::MatrixXi matches_recv;

    if (world_rank == 0){
        for (int i1=0; i1<m_intClass->sortVec_p_p.size(); i1++){
                for (int i2=0; i2<m_intClass->sortVec_ppm_hhp.size(); i2++){
                    for (int i3=0; i3<m_intClass->sortVec_pppmm_hhhpp.size(); i3++){
                        if ( m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_pppmm_hhhpp[i3] && m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_hhp[i2]){
                        matches_root.conservativeResize(3, matches_root.cols()+1);
                        matches_root.col(matches_root.cols()-1) << i1,i2,i3;
                    }
                }
            }
        }
    }

    matches_recv = distributeChannels(matches_root);

    int i1; int i2; int i3;
    for (int i=0; i<matches_recv.cols(); i++){
        i1 = matches_recv(0,i); i2 = matches_recv(1,i); i3 = matches_recv(2,i);

        Eigen::MatrixXd mat1 = m_ampClass->T5c_makemat_1(i2, i1);
        Eigen::MatrixXd mat2 = m_intClass->T5c_makemat(i2, i1);
        Eigen::MatrixXd mat3 = m_ampClass->T5c_makemat_2(i3, i1);
        Eigen::MatrixXd product = -0.5*mat3*mat1.transpose()*mat2;


        m_ampClass->T5c_inverse(product, i3, i1);
    }

    std::vector<double> TempVec;
    if (world_rank == 0){TempVec.resize(m_ampClass->T3_elements_A_temp.size()); }
    MPI_Reduce(m_ampClass->T3_elements_A_temp.data(), TempVec.data(), m_ampClass->T3_elements_A_temp.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (world_rank == 0){ m_ampClass->T3_elements_A_temp = TempVec ;}

    if (world_rank == 0){
        m_ampClass->addElementsT3_T5c();
    }*/

    int index = 0;

    Eigen::MatrixXd mat1, mat2, mat3, product, tempMat;
    for (int i1=0; i1<m_intClass->sortVec_p_p.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_hhp.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_pppmm_hhhpp.size(); i3++){
                if ( m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_pppmm_hhhpp[i3] && m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_hhp[i2]){
                    mat1 = m_ampClass->T5c_makemat_1(i2, i1);
                    mat2 = m_intClass->T5c_makemat(i2, i1);
                    //Eigen::MatrixXd mat3 = m_ampClass->T5c_makemat_2(i3, i1);
                    tempMat = mat1.transpose()*mat2;
                    mat3 = contructMatT5c(index);
                    product = -0.5*mat3*tempMat;

                    //m_ampClass->T5c_inverse(product, i3, i1);
                    m_ampClass->T5c_inverse_I(product, index);
                    //std::cout << mat2.cols() << " " << mat2.rows() << " " << mat2.cols()*mat2.rows() << std::endl;
                    index++;
                }
            }
        }
    }

    m_ampClass->addElementsT3_T5c();
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}

void Diagrams::T5d(){

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Status status;


    Eigen::MatrixXi matches_root;
    Eigen::MatrixXi matches_recv;

    if (world_rank == 0){
        for (int i1=0; i1<m_intClass->sortVec_p_p.size(); i1++){
                for (int i2=0; i2<m_intClass->sortVec_ppm_hhp.size(); i2++){
                    for (int i3=0; i3<m_intClass->sortVec_ppm_pph.size(); i3++){
                        if ( m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_pph[i3] && m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_hhp[i2]){
                        matches_root.conservativeResize(3, matches_root.cols()+1);
                        matches_root.col(matches_root.cols()-1) << i1,i2,i3;
                    }
                }
            }
        }
    }

    matches_recv = distributeChannels(matches_root, 3);

    int i1; int i2; int i3;
    for (int i=0; i<matches_recv.cols(); i++){
        i1 = matches_recv(0,i); i2 = matches_recv(1,i); i3 = matches_recv(2,i);

        Eigen::MatrixXd mat1 = m_ampClass->T5d_makemat_1(i2, i1);
        Eigen::MatrixXd mat2 = m_intClass->T5d_makemat(i2, i1);
        Eigen::MatrixXd mat3 = m_ampClass->T5d_makemat_2(i2, i3);
        Eigen::MatrixXd product = -0.5*mat1*mat2.transpose()*mat3;

        m_ampClass->T5d_inverse(product, i2, i3);
    }

    std::vector<double> TempVec;
    if (world_rank == 0){TempVec.resize(m_ampClass->T3_elements_A_temp.size()); }
    MPI_Reduce(m_ampClass->T3_elements_A_temp.data(), TempVec.data(), m_ampClass->T3_elements_A_temp.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (world_rank == 0){ m_ampClass->T3_elements_A_temp = TempVec ;}

    if (world_rank == 0){
        m_ampClass->addElementsT3_T5d();
    }

    /*for (int i1=0; i1<m_intClass->sortVec_p_p.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_hhp.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_ppm_pph.size(); i3++){
                if ( m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_pph[i3] && m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_hhp[i2]){

                    Eigen::MatrixXd mat1 = m_ampClass->T5d_makemat_1(i2, i1);
                    Eigen::MatrixXd mat2 = m_intClass->T5d_makemat(i2, i1);
                    Eigen::MatrixXd mat3 = m_ampClass->T5d_makemat_2(i2, i3);
                    Eigen::MatrixXd product = -0.5*mat1*mat2.transpose()*mat3;

                    m_ampClass->T5d_inverse(product, i2, i3);
                }
            }
        }
    }

    m_ampClass->addElementsT3_T5d();*/
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}

void Diagrams::T5e(){

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Status status;


    Eigen::MatrixXi matches_root;
    Eigen::MatrixXi matches_recv;

    if (world_rank == 0){
        for (int i1=0; i1<m_intClass->sortVec_p_h.size(); i1++){
                for (int i2=0; i2<m_intClass->sortVec_ppm_hhp.size(); i2++){
                    for (int i3=0; i3<m_intClass->sortVec_ppm_pph.size(); i3++){
                        if ( m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_pph[i3] && m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_hhp[i2]){
                        matches_root.conservativeResize(3, matches_root.cols()+1);
                        matches_root.col(matches_root.cols()-1) << i1,i2,i3;
                    }
                }
            }
        }
    }

    matches_recv = distributeChannels(matches_root, 3);

    int i1; int i2; int i3;
    for (int i=0; i<matches_recv.cols(); i++){
        i1 = matches_recv(0,i); i2 = matches_recv(1,i); i3 = matches_recv(2,i);

        Eigen::MatrixXd mat1 = m_ampClass->T5e_makemat_1(i3, i1);
        Eigen::MatrixXd mat2 = m_intClass->T5e_makemat(i3, i1);
        Eigen::MatrixXd mat3 = m_ampClass->T5e_makemat_2(i2, i3);
        Eigen::MatrixXd product = -0.5*mat3*mat2*mat1.transpose();

        m_ampClass->T5e_inverse(product, i2, i3);
    }

    std::vector<double> TempVec;
    if (world_rank == 0){TempVec.resize(m_ampClass->T3_elements_A_temp.size()); }
    MPI_Reduce(m_ampClass->T3_elements_A_temp.data(), TempVec.data(), m_ampClass->T3_elements_A_temp.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (world_rank == 0){ m_ampClass->T3_elements_A_temp = TempVec ;}

    if (world_rank == 0){
        m_ampClass->addElementsT3_T5e();
    }

    /*for (int i1=0; i1<m_intClass->sortVec_p_h.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_hhp.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_ppm_pph.size(); i3++){
                if ( m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_pph[i3] && m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_hhp[i2]){

                    Eigen::MatrixXd mat1 = m_ampClass->T5e_makemat_1(i3, i1);
                    Eigen::MatrixXd mat2 = m_intClass->T5e_makemat(i3, i1);
                    Eigen::MatrixXd mat3 = m_ampClass->T5e_makemat_2(i2, i3);
                    Eigen::MatrixXd product = -0.5*mat3*mat2*mat1.transpose();

                    m_ampClass->T5e_inverse(product, i2, i3);
                }
            }
        }
    }

    m_ampClass->addElementsT3_T5e();*/
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}

void Diagrams::T5f(){

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Status status;


    Eigen::MatrixXi matches_root;
    Eigen::MatrixXi matches_recv;

    if (world_rank == 0){
        for (int i1=0; i1<m_intClass->sortVec_pp_hh.size(); i1++){
                for (int i2=0; i2<m_intClass->sortVec_pp_pp.size(); i2++){
                    for (int i3=0; i3<m_intClass->sortVec_pppm_ppph.size(); i3++){
                        if ( m_intClass->sortVec_pp_hh[i1] == m_intClass->sortVec_pp_pp[i2] && m_intClass->sortVec_pp_hh[i1] == m_intClass->sortVec_pppm_ppph[i3]){
                        matches_root.conservativeResize(3, matches_root.cols()+1);
                        matches_root.col(matches_root.cols()-1) << i1,i2,i3;
                    }
                }
            }
        }
    }

    matches_recv = distributeChannels(matches_root, 3);

    int i1; int i2; int i3;
    for (int i=0; i<matches_recv.cols(); i++){
        i1 = matches_recv(0,i); i2 = matches_recv(1,i); i3 = matches_recv(2,i);

        Eigen::MatrixXd mat1 = m_ampClass->T5f_makemat_1(i1, i2);
        Eigen::MatrixXd mat2 = m_intClass->T5f_makemat(i1, i2);
        Eigen::MatrixXd mat3 = m_ampClass->T5f_makemat_2(i3, i1);
        Eigen::MatrixXd product = 0.25*mat1*mat2.transpose()*mat3.transpose();

        m_ampClass->T5f_inverse(product.transpose(), i3, i1);
    }

    std::vector<double> TempVec;
    if (world_rank == 0){TempVec.resize(m_ampClass->T3_elements_A_temp.size()); }
    MPI_Reduce(m_ampClass->T3_elements_A_temp.data(), TempVec.data(), m_ampClass->T3_elements_A_temp.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (world_rank == 0){ m_ampClass->T3_elements_A_temp = TempVec ;}

    if (world_rank == 0){
        m_ampClass->addElementsT3_T5f();
    }

    /*for (int i1=0; i1<m_intClass->sortVec_pp_hh.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pp_pp.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_pppm_ppph.size(); i3++){
                if ( m_intClass->sortVec_pp_hh[i1] == m_intClass->sortVec_pp_pp[i2] && m_intClass->sortVec_pp_hh[i1] == m_intClass->sortVec_pppm_ppph[i3]){

                    Eigen::MatrixXd mat1 = m_ampClass->T5f_makemat_1(i1, i2);
                    Eigen::MatrixXd mat2 = m_intClass->T5f_makemat(i1, i2);
                    Eigen::MatrixXd mat3 = m_ampClass->T5f_makemat_2(i3, i1);
                    Eigen::MatrixXd product = 0.25*mat1*mat2.transpose()*mat3.transpose();

                    m_ampClass->T5f_inverse(product.transpose(), i3, i1);
                }
            }
        }
    }

    m_ampClass->addElementsT3_T5f();*/
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}

void Diagrams::T5g(){

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Status status;


    Eigen::MatrixXi matches_root;
    Eigen::MatrixXi matches_recv;

    if (world_rank == 0){
        for (int i1=0; i1<m_intClass->sortVec_pp_hh.size(); i1++){
                for (int i2=0; i2<m_intClass->sortVec_pp_pp.size(); i2++){
                    for (int i3=0; i3<m_intClass->sortVec_pppm_hhhp.size(); i3++){
                        if ( m_intClass->sortVec_pp_hh[i1] == m_intClass->sortVec_pp_pp[i2] && m_intClass->sortVec_pp_hh[i1] == m_intClass->sortVec_pppm_hhhp[i3]){
                        matches_root.conservativeResize(3, matches_root.cols()+1);
                        matches_root.col(matches_root.cols()-1) << i1,i2,i3;
                    }
                }
            }
        }
    }

    matches_recv = distributeChannels(matches_root, 3);

    int i1; int i2; int i3;
    for (int i=0; i<matches_recv.cols(); i++){
        i1 = matches_recv(0,i); i2 = matches_recv(1,i); i3 = matches_recv(2,i);

        Eigen::MatrixXd mat1 = m_ampClass->T5g_makemat_1(i1, i2);
        Eigen::MatrixXd mat2 = m_intClass->T5g_makemat(i1, i2);
        Eigen::MatrixXd mat3 = m_ampClass->T5g_makemat_2(i3, i2);
        Eigen::MatrixXd product = 0.25*mat3*mat2.transpose()*mat1;

        m_ampClass->T5g_inverse(product, i3, i2);
    }

    std::vector<double> TempVec;
    if (world_rank == 0){TempVec.resize(m_ampClass->T3_elements_A_temp.size()); }
    MPI_Reduce(m_ampClass->T3_elements_A_temp.data(), TempVec.data(), m_ampClass->T3_elements_A_temp.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (world_rank == 0){ m_ampClass->T3_elements_A_temp = TempVec ;}

    if (world_rank == 0){
        m_ampClass->addElementsT3_T5g();
    }

    /*for (int i1=0; i1<m_intClass->sortVec_pp_hh.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pp_pp.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_pppm_hhhp.size(); i3++){
                if ( m_intClass->sortVec_pp_hh[i1] == m_intClass->sortVec_pp_pp[i2] && m_intClass->sortVec_pp_hh[i1] == m_intClass->sortVec_pppm_hhhp[i3]){

                    Eigen::MatrixXd mat1 = m_ampClass->T5g_makemat_1(i1, i2);
                    Eigen::MatrixXd mat2 = m_intClass->T5g_makemat(i1, i2);
                    Eigen::MatrixXd mat3 = m_ampClass->T5g_makemat_2(i3, i2);
                    Eigen::MatrixXd product = 0.25*mat3*mat2.transpose()*mat1;

                    m_ampClass->T5g_inverse(product, i3, i2);
                }
            }
        }
    }

    m_ampClass->addElementsT3_T5g();*/
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}