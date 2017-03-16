#include "diagrams.h"

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <omp.h>

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
    for (int n=0; n<m_intClass->Vhhpp_i.size(); n++){
        int ku = m_intClass->Vhhpp_i[n];

        int m = 0; int index;
        while (m<m_intClass->sortVec_pp_pp.size()){
            if (m_intClass->sortVec_pp_pp[m] == ku){
                index = m;
                m = m_intClass->sortVec_pp_pp.size();
            }
            m++;
        }

        //auto it = std::find(sortVec_pp_pp.begin(), sortVec_pp_pp.end(), ku);

        Eigen::MatrixXd mat1 = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T2_elements);            // t_ij^cd
        Eigen::MatrixXd mat2 = m_intClass->Vpppp[index];                                                    // v_ab^cd
        Eigen::MatrixXd product = 0.5*mat1*mat2;                                                        // (t_ij^cd)(t_cd^ab)

        m_ampClass->make2x2Block_inverse(product,ku,0,0,1,1, m_ampClass->T2_elements_new, true);
    }
    m_ampClass->addElementsT2(0,0);
    m_ampClass->T2_temp.clear();
}

void Diagrams::Lb(){
    for (int n=0; n<m_intClass->Vhhpp_i.size(); n++){
        int ku = m_intClass->Vhhpp_i[n];

        int m = 0; int index;
        while (m<m_intClass->sortVec_pp_hh.size()){
            if (m_intClass->sortVec_pp_hh[m] == ku){
                index = m;
                m = m_intClass->sortVec_pp_hh.size();
            }
            m++;
        }

        Eigen::MatrixXd mat1 = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T2_elements);            // t_kl^ab
        Eigen::MatrixXd mat2 = m_intClass->Vhhhh[index];                                                    // v_ij^kl
        Eigen::MatrixXd product = 0.5*mat2*mat1;                                                        // (v_ij^kl)(t_kl^ab)

        m_ampClass->make2x2Block_inverse(product,ku,0,0,1,1, m_ampClass->T2_elements_new, true);
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
                Eigen::MatrixXd product= -mat1*mat2;

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
        Eigen::MatrixXd product = 0.25*mat1*mat2.transpose()*mat3;                                      // (t_ij^cd)(v_cd^kl)(t_kl^ab)

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
                Eigen::MatrixXd product= 0.5*mat1*mat2.transpose()*mat3;

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
                Eigen::MatrixXd product = -0.5*mat1*mat2.transpose()*mat3;                              // (t_ijb^d)(v_d^klc)(t_klc^a)

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
                Eigen::MatrixXd product = -0.5*mat1*mat2.transpose()*mat3;                              // (t_ijb^d)(v_d^klc)(t_klc^a)

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

void Diagrams::I1_term(){
    //#pragma omp parallel for
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
                Eigen::MatrixXd product = -0.5*mat1*mat2;                              // (t_ijb^d)(v_d^klc)(t_klc^a)

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
                Eigen::MatrixXd product = -0.5*mat1.transpose()*mat2;                              // (t_ijb^d)(v_d^klc)(t_klc^a)

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
    int index = 0; int ku; int id;
    int i; int j; int k;
    int a; int b; int c;
    for (int channel = 0; channel<m_intClass->numOfKu3; channel++){
        ku = m_intClass->Vhhhppp_i[channel];

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

                m_ampClass->T3_elements_I[id] = index;
                index ++;
            }
        }
    }

    m_ampClass->T3_elements_A.resize(m_ampClass->T3_elements_I.size(), 0);
    m_ampClass->T3_elements_A_new.resize(m_ampClass->T3_elements_I.size(), 0);
    m_ampClass->T3_elements_A_temp.resize(m_ampClass->T3_elements_I.size(), 0);

    m_ampClass->T3_makeDirectMat();
}

void Diagrams::T1a(){

    for (int i1=0; i1<m_intClass->sortVec_p_p.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_pph.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_ppm_hhp.size(); i3++){
                if ( m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_pph[i2] && m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_hhp[i3]){

                    //int ku = m_intClass->sortVec_p_p[i1];
                    //Eigen::MatrixXd mat1 = m_intClass->make3x1Block(ku,1,1,0,1);
                    //Eigen::MatrixXd mat2 = m_ampClass->make3x1Block(ku,0,0,1,1, m_ampClass->T2_elements);
                    //m_ampClass->make3x3Block_inverse_I_T1a(product, ku, 0,0,1,1,1,0, m_ampClass->T3_elements_A_temp, true);

                    Eigen::MatrixXd mat1 = m_intClass->T1a_makemat(i2, i1);
                    Eigen::MatrixXd mat2 = m_ampClass->T1a_makemat(i3, i1);
                    Eigen::MatrixXd product = mat2*mat1.transpose();

                    m_ampClass->T1a_inverse(product, i3, i2);
                }
            }
        }
    }

    //m_ampClass->addElementsT3(1,0,0,1,0,0); // 0,1,1,1,1,0

    m_ampClass->addElementsT3_T1a();
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp

}

void Diagrams::T1b(){

    for (int i1=0; i1<m_intClass->sortVec_p_h.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_hhp.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_ppm_pph.size(); i3++){
                if ( m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_hhp[i2] && m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_pph[i3]){

                    //int ku = m_intClass->sortVec_p_h[i1];
                    //Eigen::MatrixXd mat1 = m_intClass->make3x1Block(ku,0,0,1,0);
                    //Eigen::MatrixXd mat2 = m_ampClass->make3x1Block(ku,1,1,0,0, m_ampClass->T2_elements);
                    //m_ampClass->make3x3Block_inverse_I_T1b(product, ku, 0,0,1,1,1,0, m_ampClass->T3_elements_A_temp, true);

                    Eigen::MatrixXd mat1 = m_intClass->T1b_makemat(i2, i1);
                    Eigen::MatrixXd mat2 = m_ampClass->T1b_makemat(i3, i1);
                    Eigen::MatrixXd product = -mat1*mat2.transpose();

                    m_ampClass->T1b_inverse(product, i2, i3);
                }
            }
        }
    }

    m_ampClass->addElementsT3_T1b();
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp

    /*for (int i1=0; i1<m_intClass->sortVec_p_h.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_hhp.size(); i2++){
            if ( m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_hhp[i2] ){
                int ku = m_intClass->sortVec_p_h[i1];

                Eigen::MatrixXd mat1 = m_intClass->make3x1Block(ku,0,0,1,0);
                Eigen::MatrixXd mat2 = m_ampClass->make3x1Block(ku,1,1,0,0, m_ampClass->T2_elements);

                Eigen::MatrixXd product = mat1*mat2.transpose();
                //std::cout << product << std::endl;
                m_ampClass->make3x3Block_inverse(product, ku, 0,0,1,1,1,0, m_ampClass->T3_elements_new, true);
            }
        }
    }
    //m_ampClass->addElementsT3(1,1,0,0,1,1); // 1,1,0,0,1,1
    m_ampClass->addElementsT3_T1b();
    m_ampClass->T3_temp.clear();*/
}

void Diagrams::T2c(){

    for (int i1=0; i1<m_intClass->sortVec_pp_pp.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pppm_hhhp.size(); i2++){
            if ( m_intClass->sortVec_pp_pp[i1] == m_intClass->sortVec_pppm_hhhp[i2]){
                Eigen::MatrixXd mat1 = m_intClass->Vpppp[i1];
                Eigen::MatrixXd mat2 = m_ampClass->T2c_makemat(i2, i1);
                Eigen::MatrixXd product = 0.5*mat2*mat1; //mathematically I need to transpose mat1, but it's symmetric

                m_ampClass->T2c_inverse(product, i2, i1);
            }
        }
    }

    //m_ampClass->addElementsT3_T2c();
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}

void Diagrams::T2d(){

    for (int i1=0; i1<m_intClass->sortVec_pp_hh.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pppm_ppph.size(); i2++){
            if ( m_intClass->sortVec_pp_hh[i1] == m_intClass->sortVec_pppm_ppph[i2]){
                Eigen::MatrixXd mat1 = m_intClass->Vhhhh[i1];
                Eigen::MatrixXd mat2 = m_ampClass->T2d_makemat(i2, i1);
                Eigen::MatrixXd product = 0.5*mat2*mat1; //mathematically I need to transpose mat1, but it's symmetric

                m_ampClass->T2d_inverse(product, i2, i1);
            }
        }
    }

    //m_ampClass->addElementsT3_T2d();
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}

void Diagrams::T2e(){

    for (int i1=0; i1<m_intClass->sortVec_pm_hp.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppmm_hhpp.size(); i2++){
            if ( m_intClass->sortVec_pm_hp[i1] == m_intClass->sortVec_ppmm_hhpp[i2]){
                Eigen::MatrixXd mat1 = m_intClass->Vhphp[i1];
                Eigen::MatrixXd mat2 = m_ampClass->T2e_makemat(i2, i1);
                Eigen::MatrixXd product = mat2*mat1; //mathematically I need to transpose mat1, but it's symmetric

                m_ampClass->T2c_inverse(product, i2, i1);
            }
        }
    }

    //m_ampClass->addElementsT3_T2c();
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
                    Eigen::MatrixXd product = mat2*mat1.transpose();

                    m_ampClass->T3b_inverse(product, i3, i1);
                }
            }
        }
    }

    m_ampClass->addElementsT3_T3b();
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}

/*void Diagrams::T5a(){

    for (int i1=0; i1<m_intClass->sortVec_p_p.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_pph.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_ppm_hhp.size(); i3++){
                if ( m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_pph[i2] && m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_hhp[i3]){

                    Eigen::MatrixXd mat1 = m_intClass->T3b_makemat(i2, i1);
                    Eigen::MatrixXd mat2 = m_ampClass->T3b_makemat(i3, i1);
                    Eigen::MatrixXd product = mat2*mat1.transpose();

                    m_ampClass->T3b_inverse(product, i3, i2);
                }
            }
        }
    }

    //m_ampClass->addElementsT3_T1a();
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}*/

