#include "diagrams.h"

#include <eigen3/Eigen/Dense>

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

        Eigen::MatrixXd mat1 = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T2_elements);            // t_ij^cd
        Eigen::MatrixXd mat2 = m_intClass->Vpppp[n];                                                    // v_ab^cd
        Eigen::MatrixXd product = 0.5*mat1*mat2;                                                        // (t_ij^cd)(t_cd^ab)

        m_ampClass->make2x2Block_inverse(product,ku,0,0,1,1, m_ampClass->T2_elements_new, true);
    }
    m_ampClass->addElementsT2(0,0);
    m_ampClass->T2_temp.clear();
}

void Diagrams::Lb(){
    for (int n=0; n<m_intClass->Vhhpp_i.size(); n++){
        int ku = m_intClass->Vhhpp_i[n];

        Eigen::MatrixXd mat1 = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T2_elements);            // t_kl^ab
        Eigen::MatrixXd mat2 = m_intClass->Vhhhh[n];                                                    // v_ij^kl
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

void Diagrams::D10b(){
    for (int i1=0; i1<m_intClass->sortVec_p_p.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_pph.size(); i2++){
            if ( m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_pph[i2] ){
                int ku = m_intClass->sortVec_p_p[i1];
                Eigen::MatrixXd mat1 = m_ampClass->make3x3Block(ku,0,0,1,1,1,0, m_ampClass->T3_elements);    // t_ijb^d
                Eigen::MatrixXd mat2 = m_intClass->make3x1Block(ku,1,1,0,1);                            // v_klc^d
                Eigen::MatrixXd product = -0.5*mat1*mat2;                              // (t_ijb^d)(v_d^klc)(t_klc^a)

                m_ampClass->make3x1Block_inverse(product, ku, 0,0,1,1, m_ampClass->T2_elements_new, true);
            }
        }
    }
    m_ampClass->addElementsT2(0,1);
    m_ampClass->T2_temp.clear();
}

void Diagrams::D10c(){
    for (int i1=0; i1<m_intClass->sortVec_p_h.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_hhp.size(); i2++){
            if ( m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_hhp[i2] ){
                int ku = m_intClass->sortVec_p_h[i1];
                Eigen::MatrixXd mat1 = m_ampClass->make3x3Block(ku,0,0,1,1,1,0, m_ampClass->T3_elements);    // t_ijb^d
                Eigen::MatrixXd mat2 = m_intClass->make3x1Block(ku,0,0,1,0);                            // v_klc^d
                Eigen::MatrixXd product = 0.5*mat1.transpose()*mat2;                              // (t_ijb^d)(v_d^klc)(t_klc^a)
                //std::cout << mat1 << std::endl;
                m_ampClass->make3x1Block_inverse(product, ku, 1,1,0,0, m_ampClass->T2_elements_new, true);
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

void Diagrams::setFirstIt(bool argument){
    m_firstIt = argument;
}

void Diagrams::T1a(){
    /*if (m_firstIt == true){
        for (int i1=0; i1<m_intClass->sortVec_p_p.size(); i1++){
            for (int i2=0; i2<m_intClass->sortVec_ppm_pph.size(); i2++){
                if ( m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_pph[i2] ){
                    int ku = m_intClass->sortVec_p_p[i1];

                    T1a_V.push_back( m_intClass->make3x1Block(ku,1,1,0,1) );
                    //Eigen::MatrixXd mat1 = m_intClass->make3x1Block(ku,1,1,0,1);

                    T1a_T2.push_back( m_ampClass->make3x1Block(ku,0,0,1,1, m_ampClass->T2_elements) );
                    //Eigen::MatrixXd mat2 = m_ampClass->make3x1Block(ku,0,0,1,1, m_ampClass->T2_elements);

                    Eigen::MatrixXd product = T1a_T2.back()*( T1a_V.back() ).transpose();
                    //std::cout << product << std::endl;
                    m_ampClass->make3x3Block_inverse(product, ku, 0,0,1,1,1,0, m_ampClass->T3_elements_new, true);
                }
            }
        }
    }*/
    //else{
        for (int i1=0; i1<m_intClass->sortVec_p_p.size(); i1++){
            for (int i2=0; i2<m_intClass->sortVec_ppm_pph.size(); i2++){
                if ( m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_pph[i2] ){
                    int ku = m_intClass->sortVec_p_p[i1];

                    Eigen::MatrixXd mat1 = m_intClass->make3x1Block(ku,1,1,0,1);
                    Eigen::MatrixXd mat2 = m_ampClass->make3x1Block(ku,0,0,1,1, m_ampClass->T2_elements);

                    Eigen::MatrixXd product = mat2*mat1.transpose();
                    //std::cout << product << std::endl;
                    m_ampClass->make3x3Block_inverse(product, ku, 0,0,1,1,1,0, m_ampClass->T3_elements_new, true);
                }
            }
        }
    //}

    //m_ampClass->addElementsT3(1,0,0,1,0,0); // 0,1,1,1,1,0
    m_ampClass->addElementsT3_T1a();
    m_ampClass->T3_temp.clear();
}

void Diagrams::T1b(){
    for (int i1=0; i1<m_intClass->sortVec_p_h.size(); i1++){
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
    m_ampClass->T3_temp.clear();
}

// ##################################################
// ##                                              ##
// ## INTERMEDIATES                                ##
// ##                                              ##
// ##################################################

void Diagrams::I1_term(){
    for (int n=0; n<m_intClass->Vhhpp_i.size(); n++){
        int ku = m_intClass->Vhhpp_i[n];

        Eigen::MatrixXd mat1 = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T2_elements);
        Eigen::MatrixXd mat2 = m_intClass->make2x2Block(ku,0,0,1,1);

        Eigen::MatrixXd I1   = 0.5*mat1*mat2.transpose();

        Eigen::MatrixXd mat3 = m_intClass->Vhhhh[n].transpose();

        I1 += mat3;

        Eigen::MatrixXd product = 0.5*I1*mat1;

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
