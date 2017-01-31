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

void Diagrams::La(){
    for (int n=0; n<m_intClass->Vhhpp_i.size(); n++){
        int ku = m_intClass->Vhhpp_i[n];

        Eigen::MatrixXd mat1 = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T_elements);            // t_ij^cd
        Eigen::MatrixXd mat2 = m_intClass->Vpppp[n];                                                    // v_ab^cd
        Eigen::MatrixXd product = 0.5*mat1*mat2;                                                        // (t_ij^cd)(t_cd^ab)

        m_ampClass->make2x2Block_inverse(product,ku,0,0,1,1, m_ampClass->T_elements_new, true);
    }
    m_ampClass->addElements(0,0);
    m_ampClass->T_temp.clear();
}

void Diagrams::Lb(){
    for (int n=0; n<m_intClass->Vhhpp_i.size(); n++){
        int ku = m_intClass->Vhhpp_i[n];

        Eigen::MatrixXd mat1 = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T_elements);            // t_kl^ab
        Eigen::MatrixXd mat2 = m_intClass->Vhhhh[n];                                                    // v_ij^kl
        Eigen::MatrixXd product = 0.5*mat2*mat1;                                                        // (v_ij^kl)(t_kl^ab)

        m_ampClass->make2x2Block_inverse(product,ku,0,0,1,1, m_ampClass->T_elements_new, true);
    }
    m_ampClass->addElements(0,0);
    m_ampClass->T_temp.clear();
}

void Diagrams::Lc(){
    for (int i1=0; i1<m_intClass->sortVec_hp.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ph.size(); i2++){
            if ( m_intClass->sortVec_hp[i1] == m_intClass->sortVec_ph[i2] ){
                int ku = m_intClass->sortVec_ph[i2];

                Eigen::MatrixXd mat1 = m_intClass->Vhphp[i1];
                Eigen::MatrixXd mat2 = m_ampClass->make2x2Block(ku,0,1,1,0, m_ampClass->T_elements);
                Eigen::MatrixXd product= -mat1*mat2;

                m_ampClass->make2x2Block_inverse(product, ku, 0,1,1,0, m_ampClass->T_elements_new, true);
            }
        }
    }
    m_ampClass->addElements(1,1);
    m_ampClass->T_temp.clear();
}

void Diagrams::Qa(){
    for (int n=0; n<m_intClass->numOfKu; n++){
        int ku = m_intClass->Vhhpp_i[n];

        Eigen::MatrixXd mat1 = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T_elements);            // t_ij^ab
        Eigen::MatrixXd mat2 = m_intClass->make2x2Block(ku,0,0,1,1);                                    // v_ij^ab
        Eigen::MatrixXd mat3 = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T_elements);            // t_ij^ab
        Eigen::MatrixXd product = 0.25*mat1*mat2.transpose()*mat3;                                      // (t_ij^cd)(v_cd^kl)(t_kl^ab)

        m_ampClass->make2x2Block_inverse(product,ku,0,0,1,1, m_ampClass->T_elements_new, true);
    }
    m_ampClass->addElements(0,0);
    m_ampClass->T_temp.clear();
}

void Diagrams::Qb(){
    for (int i1=0; i1<m_intClass->sortVec_hp.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ph.size(); i2++){
            if ( m_intClass->sortVec_hp[i1] == m_intClass->sortVec_ph[i2] ){
                int ku = m_intClass->sortVec_hp[i1];

                Eigen::MatrixXd mat1 = m_ampClass->make2x2Block(ku,0,1,1,0, m_ampClass->T_elements);
                Eigen::MatrixXd mat2 = m_intClass->make2x2Block(ku,0,1,1,0);
                Eigen::MatrixXd mat3 = mat1;
                Eigen::MatrixXd product= 0.5*mat1*mat2.transpose()*mat3;

                m_ampClass->make2x2Block_inverse(product, ku, 0,1,1,0, m_ampClass->T_elements_new, true);
            }
        }
    }
    m_ampClass->addElements(1,1);
    m_ampClass->T_temp.clear();
}

void Diagrams::Qc(){
    for (int i1=0; i1<m_intClass->sortVec_p.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_hhp.size(); i2++){
            if ( m_intClass->sortVec_p[i1] == m_intClass->sortVec_hhp[i2] ){
                int ku = m_intClass->sortVec_p[i1];

                Eigen::MatrixXd mat1 = m_ampClass->make3x1Block(ku,0,0,1,1, m_ampClass->T_elements);    // t_ijb^d
                Eigen::MatrixXd mat2 = m_intClass->make3x1Block(ku,0,0,1,1);                            // v_klc^d
                Eigen::MatrixXd mat3 = mat1;                                                            // t_klc^a
                Eigen::MatrixXd product = -0.5*mat1*mat2.transpose()*mat3;                              // (t_ijb^d)(v_d^klc)(t_klc^a)

                m_ampClass->make3x1Block_inverse(product, ku, 0,0,1,1, m_ampClass->T_elements_new, true);
            }
        }
    }
    m_ampClass->addElements(0,1);
    m_ampClass->T_temp.clear();
}

void Diagrams::Qd(){
    for (int i1=0; i1<m_intClass->sortVec_h.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pph.size(); i2++){
            if ( m_intClass->sortVec_h[i1] == m_intClass->sortVec_pph[i2] ){
                int ku = m_intClass->sortVec_h[i1];
                Eigen::MatrixXd mat1 = m_ampClass->make3x1Block(ku,1,1,0,0, m_ampClass->T_elements);    // t_ijb^d
                Eigen::MatrixXd mat2 = m_intClass->make3x1Block(ku,1,1,0,0);                            // v_klc^d
                Eigen::MatrixXd mat3 = mat1;                                                            // t_klc^a
                Eigen::MatrixXd product = -0.5*mat1*mat2.transpose()*mat3;                              // (t_ijb^d)(v_d^klc)(t_klc^a)

                m_ampClass->make3x1Block_inverse(product, ku, 1,1,0,0, m_ampClass->T_elements_new, true);
            }
        }
    }
    m_ampClass->addElements(1,0);
    m_ampClass->T_temp.clear();
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
