#include "diagrams.h"

typedef std::chrono::high_resolution_clock Clock;

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

void Diagrams::setThreads(int numthreads){
    m_numThreads = numthreads;
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

                MatrixX mat1 = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T2_elements);            // t_ij^cd
                MatrixX mat2 = m_intClass->Vpppp[i2];                                                    // v_ab^cd
                MatrixX M1(mat1.rows(), 2*mat1.cols());
                M1 << mat1, -mat1;
                MatrixX M2(2*mat2.rows(), mat2.cols());
                M2 << mat2, -mat2;
                
                MatrixX product = 0.5*(M1*M2);

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

                MatrixX mat1 = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T2_elements);            // t_kl^ab
                MatrixX mat2 = m_intClass->Vhhhh[i1];                                                    // v_ij^kl
                MatrixX product = 0.5*(mat2*mat1);                                                        // (v_ij^kl)(t_kl^ab)

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

                MatrixX mat1 = m_intClass->Vhphp[i1];
                MatrixX mat2 = m_ampClass->make2x2Block(ku,0,1,1,0, m_ampClass->T2_elements);
                MatrixX product= -1*(mat1*mat2);

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

        MatrixX mat1 = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T2_elements);            // t_ij^ab
        MatrixX mat2 = m_intClass->make2x2Block(ku,0,0,1,1);                                    // v_ij^ab
        MatrixX mat3 = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T2_elements);            // t_ij^ab
        MatrixX product = 0.25*(mat1*mat2.transpose()*mat3);                                      // (t_ij^cd)(v_cd^kl)(t_kl^ab)

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

                MatrixX mat1 = m_ampClass->make2x2Block(ku,0,1,1,0, m_ampClass->T2_elements);
                MatrixX mat2 = m_intClass->make2x2Block(ku,0,1,1,0);
                MatrixX mat3 = mat1;
                MatrixX product= 0.5*(mat1*mat2.transpose()*mat3);

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

                MatrixX mat1 = m_ampClass->make3x1Block(ku,0,0,1,1, m_ampClass->T2_elements);    // t_ijb^d
                MatrixX mat2 = m_intClass->make3x1Block(ku,0,0,1,1);                            // v_klc^d
                MatrixX mat3 = mat1;                                                            // t_klc^a
                MatrixX product = -0.5*(mat1*mat2.transpose()*mat3);                              // (t_ijb^d)(v_d^klc)(t_klc^a)

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
                MatrixX mat1 = m_ampClass->make3x1Block(ku,1,1,0,0, m_ampClass->T2_elements);    // t_ijb^d
                MatrixX mat2 = m_intClass->make3x1Block(ku,1,1,0,0);                            // v_klc^d
                MatrixX mat3 = mat1;                                                            // t_klc^a
                MatrixX product = -0.5*(mat1*mat2.transpose()*mat3);                              // (t_ijb^d)(v_d^klc)(t_klc^a)

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

                MatrixX mat1 = m_ampClass->I1_makemat_1(i1,i2);
                MatrixX M1(mat1.rows(), 2*mat1.cols());
                M1 << mat1, -mat1;

                MatrixX mat2 = m_intClass->I1_makemat(i1,i2);
                MatrixX M2(mat2.rows(), 2*mat2.cols());
                M2 << mat2, -mat2;

                MatrixX I1   = 0.5*(M1*M2.transpose());

                MatrixX mat3 = m_intClass->Vhhhh[i1];

                I1.noalias() += mat3;

                MatrixX I1M(I1.rows(), 2*I1.cols());
                I1M << I1, -I1;

                MatrixX mat4 = m_ampClass->I1_makemat_2(i1,i2);
                MatrixX M4(2*mat4.rows(), mat4.cols());
                M4 << mat4, -mat4;

                MatrixX product = 0.5*(I1M*M4);

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

                MatrixX mat1 = m_ampClass->I2_makemat_1(i1,i2);
                MatrixX mat2 = m_intClass->I2_makemat(i1,i2);

                MatrixX I2   = 0.5*(mat1*mat2.transpose());

                MatrixX mat3 = m_intClass->Vhphp[i1];

                I2.noalias() += mat3; //the minus is intentional

                MatrixX mat4 = m_ampClass->I2_makemat_2(i1,i2);

                MatrixX product = I2*mat4;

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
                MatrixX mat1 = m_ampClass->I3_makemat_1(i2,i1);

                MatrixX mat2 = m_intClass->I3_makemat(i2,i1);
                MatrixX M2(2*mat2.rows(), mat2.cols());
                M2 << mat2, -mat2;

                MatrixX mat3 = m_ampClass->I3_makemat_2(i2,i1);
                MatrixX M3(2*mat3.rows(), mat3.cols());
                M3 << mat3, -mat3;

                MatrixX product   = -0.5*(M3.transpose()*M2*mat1.transpose());

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

                MatrixX mat1 = m_ampClass->I4_makemat_1(i2,i1);
                MatrixX M1(2*mat1.rows(), mat1.cols());
                M1 << mat1, -mat1;


                MatrixX mat2 = m_intClass->I4_makemat(i2,i1);
                MatrixX M2(2*mat2.rows(), mat2.cols());
                M2 << mat2, -mat2;

                MatrixX mat3 = m_ampClass->I4_makemat_2(i2,i1);

                MatrixX product   = -0.5*(mat3*M2.transpose()*M1);

                m_ampClass->I4_inverse(product, i2, i1);
            }
        }
    }
    m_ampClass->addElementsT2(0,1);
    m_ampClass->T2_temp.clear();
}

void Diagrams::I1_term(){
    for (int n=0; n<m_intClass->Vhhpp_i.size(); n++){
        int ku = m_intClass->Vhhpp_i[n];

        MatrixX mat1 = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T2_elements);
        MatrixX mat2 = m_intClass->make2x2Block(ku,0,0,1,1);

        MatrixX I1   = 0.5*mat1*mat2.transpose();

        MatrixX mat3 = m_intClass->Vhhhh[n].transpose();

        I1 += mat3;

        MatrixX product = 0.5*I1*mat1;

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


                MatrixX mat1 = m_ampClass->make2x2Block(ku,0,1,1,0, m_ampClass->T2_elements);
                MatrixX mat2 = m_intClass->make2x2Block(ku,0,1,1,0);

                MatrixX I2   = 0.5*mat2*mat1.transpose();

                MatrixX mat3 = m_intClass->Vhphp[i1];

                I2 -= mat3;

                MatrixX product = I2.transpose()*mat1;

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

                MatrixX mat1 = m_ampClass->make3x1Block(ku,1,1,0,0, m_ampClass->T2_elements);
                MatrixX mat2 = m_intClass->make3x1Block(ku,1,1,0,0);

                MatrixX I3   = mat2.transpose()*mat1;

                MatrixX product = -0.5*mat1*I3;

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

                MatrixX mat1 = m_ampClass->make3x1Block(ku,0,0,1,1, m_ampClass->T2_elements);
                MatrixX mat2 = m_intClass->make3x1Block(ku,0,0,1,1);

                MatrixX I3   = mat2.transpose()*mat1;

                MatrixX product = -0.5*mat1*I3;

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

    Eigen::MatrixXi matches;

    std::vector<MatrixX> mat2v;
    std::vector<MatrixX> productv;

    for (int i1=0; i1<m_intClass->sortVec_p_p.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_pph.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_ppm_hhp.size(); i3++){
                if ( m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_pph[i2] && m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_hhp[i3]){
                    matches.conservativeResize(3, matches.cols()+1);
                    matches.col(matches.cols()-1) << i1,i2,i3;

                    MatrixX mat = m_intClass->D10b_makemat(i2, i1);
                    MatrixX M(2*mat.rows(), mat.cols());
                    M << mat, -mat;

                    mat2v.push_back( M );
                }
            }
        }
    }

    productv.resize(matches.cols());

    int i1; int i2; int i3; int cols = matches.cols();
#pragma omp parallel for num_threads(m_numThreads) private(i1,i2,i3) firstprivate(cols)
    for (int i=0; i<cols; i++){
        i1 = matches(0,i); i2 = matches(1,i); i3 = matches(2,i);

        MatrixX mat1 = m_ampClass->D10b_makemat(i3, i2);
        MatrixX M1(mat1.rows(), 2*mat1.cols());
        M1 << mat1, -mat1;

        productv[i] = -0.5*(M1*mat2v[i]);
    }

    for (int i=0; i<cols; i++){
        i1 = matches(0,i); i3 = matches(2,i);
        m_ampClass->D10b_inverse(productv[i], i3, i1);
    }

    m_ampClass->addElementsT2(0,1);
    m_ampClass->T2_temp.clear();
}

void Diagrams::D10c(){

    Eigen::MatrixXi matches;

    std::vector<MatrixX> mat2v;
    std::vector<MatrixX> productv;

    for (int i1=0; i1<m_intClass->sortVec_p_h.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_hhp.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_ppm_pph.size(); i3++){
            if ( m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_hhp[i2] && m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_pph[i3]){
                    matches.conservativeResize(3, matches.cols()+1);
                    matches.col(matches.cols()-1) << i1,i2,i3;

                    MatrixX mat = m_intClass->D10c_makemat(i2, i1);
                    MatrixX M(2*mat.rows(), mat.cols());
                    M << mat, -mat;

                    mat2v.push_back( M );
                }
            }
        }
    }

    productv.resize(matches.cols());

    int i1; int i2; int i3; int cols = matches.cols();
#pragma omp parallel for num_threads(m_numThreads) private(i1,i2,i3) firstprivate(cols)
    for (int i=0; i<cols; i++){
        i1 = matches(0,i); i2 = matches(1,i); i3 = matches(2,i);

        MatrixX mat1 = m_ampClass->D10c_makemat(i2, i3);
        MatrixX M1(2*mat1.rows(), mat1.cols());
        M1 << mat1, -mat1;

        productv[i] = -0.5*(M1.transpose()*mat2v[i]);
    }

    for (int i=0; i<cols; i++){
        i1 = matches(0,i); i3 = matches(2,i);
        m_ampClass->D10c_inverse(productv[i], i3, i1);
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
    unsigned long int index = 0; unsigned long int id;
    int i; int j; int k;
    int a; int b; int c;

    for (int channel = 0; channel<m_intClass->numOfKu3; channel++){

        unsigned long int range_lower1 = m_intClass->boundsHolder_hhhppp_hhh(0,channel);
        unsigned long int range_upper1 = m_intClass->boundsHolder_hhhppp_hhh(1,channel);
        unsigned long int range_lower2 = m_intClass->boundsHolder_hhhppp_ppp(0,channel);
        unsigned long int range_upper2 = m_intClass->boundsHolder_hhhppp_ppp(1,channel);

        for (unsigned long int hhh=range_lower1; hhh<range_upper1; hhh++){
            i = (m_intClass->blockArrays_ppp_hhh)(0,hhh);
            j = (m_intClass->blockArrays_ppp_hhh)(1,hhh);
            k = (m_intClass->blockArrays_ppp_hhh)(2,hhh);
            for (unsigned long int ppp=range_lower2; ppp<range_upper2; ppp++){
                a = (m_intClass->blockArrays_ppp_ppp)(0,ppp);
                b = (m_intClass->blockArrays_ppp_ppp)(1,ppp);
                c = (m_intClass->blockArrays_ppp_ppp)(2,ppp);

                id = m_intClass->Identity_hhhppp(i,j,k,a,b,c);
                m_ampClass->T3_elements_I[id] = index;

                index ++;
            }
        }
    }

    m_ampClass->T3_elements_A.resize(m_ampClass->T3_elements_I.size(), 0);
    m_ampClass->T3_elements_A_new.resize(m_ampClass->T3_elements_I.size(), 0);

    std::cout << "Number of T3 elements: " << m_ampClass->T3_elements_A.size() << std::endl;
}

void Diagrams::T1a(){

    Eigen::MatrixXi matches;

    std::vector<MatrixX> mat1v;
    std::vector<MatrixX> mat2v;

    for (int i1=0; i1<m_intClass->sortVec_p_p.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_pph.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_ppm_hhp.size(); i3++){
                if ( m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_pph[i2] && m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_hhp[i3]){
                    matches.conservativeResize(3, matches.cols()+1);
                    matches.col(matches.cols()-1) << i1,i2,i3;

                    mat1v.push_back( m_intClass->T1a_makemat(i2, i1) );
                    mat2v.push_back( m_ampClass->T1a_makemat(i3, i1) );
                }
            }
        }
    }

    int i1; int i2; int i3; int cols = matches.cols();
#pragma omp parallel for num_threads(m_numThreads) private(i1,i2,i3) firstprivate(cols)
    for (int i=0; i<cols; i++){
        i1 = matches(0,i); i2 = matches(1,i); i3 = matches(2,i);

        MatrixX product = mat2v[i]*mat1v[i].transpose();

        m_ampClass->T1a_inverse(product, i3, i2);
    }
    
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}


void Diagrams::T1b(){

    Eigen::MatrixXi matches;

    std::vector<MatrixX> mat1v;
    std::vector<MatrixX> mat2v;

    for (int i1=0; i1<m_intClass->sortVec_p_h.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_hhp.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_ppm_pph.size(); i3++){
                if ( m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_hhp[i2] && m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_pph[i3]){
                    matches.conservativeResize(3, matches.cols()+1);
                    matches.col(matches.cols()-1) << i1,i2,i3;

                    mat1v.push_back( m_intClass->T1b_makemat(i2, i1) );
                    mat2v.push_back( m_ampClass->T1b_makemat(i3, i1) );
                }
            }
        }
    }


    int i1; int i2; int i3; int cols=matches.cols();
#pragma omp parallel for num_threads(m_numThreads) private(i1,i2,i3) firstprivate(cols)
    for (int i=0; i<cols; i++){
        i1 = matches(0,i); i2 = matches(1,i); i3 = matches(2,i);

        MatrixX product = -1*(mat1v[i]*mat2v[i].transpose());

        m_ampClass->T1b_inverse(product, i2, i3);
    }

    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}

void Diagrams::T2c(){

    Eigen::MatrixXi matches;

    for (int i1=0; i1<m_intClass->sortVec_pp_pp.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pppm_hhhp.size(); i2++){
            if ( m_intClass->sortVec_pp_pp[i1] == m_intClass->sortVec_pppm_hhhp[i2]){
                matches.conservativeResize(2, matches.cols()+1);
                matches.col(matches.cols()-1) << i1,i2;
            }
        }
    }

    int i1; int i2; int cols = matches.cols();
#pragma omp parallel for num_threads(m_numThreads) private(i1,i2) firstprivate(cols)
    for (int i=0; i<cols; i++){
        i1 = matches(0,i); i2 = matches(1,i);

        MatrixX mat1 = m_intClass->Vpppp[i1];
        MatrixX M1(2*mat1.rows(), mat1.cols());
        M1 << mat1, -mat1;

        MatrixX mat2 = m_ampClass->T2c_makemat(i2, i1);
        MatrixX M2(mat2.rows(), 2*mat2.cols());
        M2 << mat2, -mat2;

        MatrixX product = 0.5*M2*M1; //mathematically I need to transpose mat1, but it's symmetric

        m_ampClass->T2c_inverse(product, i2, i1);
    }

    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}

void Diagrams::T2d(){

    Eigen::MatrixXi matches;

    for (int i1=0; i1<m_intClass->sortVec_pp_hh.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pppm_ppph.size(); i2++){
            if ( m_intClass->sortVec_pp_hh[i1] == m_intClass->sortVec_pppm_ppph[i2]){
                matches.conservativeResize(2, matches.cols()+1);
                matches.col(matches.cols()-1) << i1,i2;
            }
        }
    }


    int i1; int i2; int cols = matches.cols();
#pragma omp parallel for num_threads(m_numThreads) private(i1,i2) firstprivate(cols)
    for (int i=0; i<cols; i++){
        i1 = matches(0,i); i2 = matches(1,i);

        MatrixX mat1 = m_intClass->Vhhhh[i1];
        MatrixX M1(2*mat1.rows(), mat1.cols());
        M1 << mat1, -mat1;

        MatrixX mat2 = m_ampClass->T2d_makemat(i2, i1);
        MatrixX M2(mat2.rows(), 2*mat2.cols());
        M2 << mat2, -mat2;

        MatrixX product = 0.5*M2*M1; //mathematically I need to transpose mat1, but it's symmetric

        m_ampClass->T2d_inverse(product, i2, i1);
    }

    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}

void Diagrams::T2e(){

    Eigen::MatrixXi matches;

    for (int i1=0; i1<m_intClass->sortVec_pm_hp.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppmm_pphh.size(); i2++){
            if ( m_intClass->sortVec_pm_hp[i1] == m_intClass->sortVec_ppmm_pphh[i2]){
                matches.conservativeResize(2, matches.cols()+1);
                matches.col(matches.cols()-1) << i1,i2;
            }
        }
    }

    int i1; int i2; int cols = matches.cols();
#pragma omp parallel for num_threads(m_numThreads) private(i1,i2) firstprivate(cols)
    for (int i=0; i<cols; i++){
#pragma omp critical
    {
        i1 = matches(0,i); i2 = matches(1,i);
        MatrixX mat1 = m_intClass->Vhphp[i1]; //hphp was made with sign index +- on rows and columns
        MatrixX mat2 = m_ampClass->T2e_makemat(i2, i1);
        MatrixX product = mat2*mat1; //mathematically I need to transpose mat1, but it's symmetric
        m_ampClass->T2e_inverse(product, i2, i1);
    }
    }
    
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}

void Diagrams::T3b(){

    for (int i1=0; i1<m_intClass->sortVec_pm_hp.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pm_ph.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_pm_pp.size(); i3++){
                if ( m_intClass->sortVec_pm_hp[i1] == m_intClass->sortVec_pm_ph[i2] && m_intClass->sortVec_pm_hp[i1] == m_intClass->sortVec_pm_pp[i3]){

                    MatrixX mat1 = m_ampClass->T3b_makemat_1(i1, i2);
                    MatrixX mat2 = m_intClass->T3b_makemat(i3, i1);

                    MatrixX product = mat2*mat1;

                    m_ampClass->T3b_Inverse_temp(product, i3, i2);
                }
            }
        }
    }

    //now use remapped function
    //we now parallellize since we are mapping back T3 amplitudes (demanding)

    Eigen::MatrixXi matches;

    std::vector<MatrixX> mat1v;
    std::vector<MatrixX> mat2v;

    for (int i1=0; i1<m_intClass->sortVec_ppm_pph.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_p_p.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_ppm_hhp.size(); i3++){
                if ( m_intClass->sortVec_ppm_pph[i1] == m_intClass->sortVec_p_p[i2] && m_intClass->sortVec_ppm_pph[i1] == m_intClass->sortVec_ppm_hhp[i3]){
                    matches.conservativeResize(3, matches.cols()+1);
                    matches.col(matches.cols()-1) << i1,i2,i3;

                    mat1v.push_back( m_ampClass->T3b_makemat_2(i1, i2) );
                    mat2v.push_back( m_ampClass->T3b_makemat_3(i3, i2) );
                }
            }
        }
    }


    int i1; int i2; int i3; int cols = matches.cols();
#pragma omp parallel for num_threads(m_numThreads) private(i1,i2,i3) firstprivate(cols)
    for (int i=0; i<cols; i++){
        i1 = matches(0,i); i2 = matches(1,i); i3 = matches(2,i);

        MatrixX product = -mat2v[i]*mat1v[i].transpose();

        m_ampClass->T3b_inverse(product, i3, i1);
    }
    
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
    m_ampClass->T3D_remap.clear();

}

void Diagrams::T3c(){

    for (int i1=0; i1<m_intClass->sortVec_pm_hp.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pm_ph.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_pm_hh.size(); i3++){
                if ( m_intClass->sortVec_pm_hp[i1] == m_intClass->sortVec_pm_ph[i2] && m_intClass->sortVec_pm_hp[i1] == m_intClass->sortVec_pm_hh[i3]){

                    MatrixX mat1 = m_ampClass->T3c_makemat_1(i1, i2);
                    MatrixX mat2 = m_intClass->T3c_makemat(i3, i2);
                    MatrixX product = mat1*mat2.transpose();

                    m_ampClass->T3c_Inverse_temp(product.transpose(), i3, i1);
                }
            }
        }
    }


    //now use remapped function
    //we now parallellize since we are mapping back T3 amplitudes (demanding)

    Eigen::MatrixXi matches;

    std::vector<MatrixX> mat1v;
    std::vector<MatrixX> mat2v;

    for (int i1=0; i1<m_intClass->sortVec_ppm_pph.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_p_h.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_ppm_hhp.size(); i3++){
                if ( m_intClass->sortVec_ppm_pph[i1] == m_intClass->sortVec_p_h[i2] && m_intClass->sortVec_ppm_pph[i1] == m_intClass->sortVec_ppm_hhp[i3]){
                    matches.conservativeResize(3, matches.cols()+1);
                    matches.col(matches.cols()-1) << i1,i2,i3;

                    mat1v.push_back( m_ampClass->T3c_makemat_2(i3, i2) );
                    mat2v.push_back( m_ampClass->T3c_makemat_3(i1, i2) );
                }
            }
        }
    }


    int i1; int i2; int i3; int cols = matches.cols();
#pragma omp parallel for num_threads(m_numThreads) private(i1,i2,i3) firstprivate(cols)
    for (int i=0; i<cols; i++){
        i1 = matches(0,i); i2 = matches(1,i); i3 = matches(2,i);

        MatrixX product = mat1v[i]*mat2v[i].transpose();

        m_ampClass->T3c_inverse(product, i3, i1);
    }
    
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
    m_ampClass->T3D_remap.clear();

}

void Diagrams::T3d(){

    for (int i1=0; i1<m_intClass->sortVec_pp_hh.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pp_ph.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_pp_pp.size(); i3++){
                if ( m_intClass->sortVec_pp_hh[i1] == m_intClass->sortVec_pp_ph[i2] && m_intClass->sortVec_pp_hh[i1] == m_intClass->sortVec_pp_pp[i3]){

                    MatrixX mat1 = m_ampClass->T3d_makemat_1(i1, i3);
                    MatrixX M1(mat1.rows(), 2*mat1.cols());
                    M1 << mat1, -mat1;

                    MatrixX mat2 = m_intClass->T3d_makemat(i3, i2);
                    MatrixX M2(2*mat2.rows(), mat2.cols());
                    M2 << mat2, -mat2;

                    MatrixX product = M1*M2;

                    m_ampClass->T3d_Inverse_temp(product, i1, i2);
                }
            }
        }
    }

    //now use remapped function
    //we now parallellize since we are mapping back T3 amplitudes (demanding)

    Eigen::MatrixXi matches;

    std::vector<MatrixX> mat1v;
    std::vector<MatrixX> mat2v;

    for (int i1=0; i1<m_intClass->sortVec_ppm_pph.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_p_h.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_ppm_hhp.size(); i3++){
                if ( m_intClass->sortVec_ppm_pph[i1] == m_intClass->sortVec_p_h[i2] && m_intClass->sortVec_ppm_pph[i1] == m_intClass->sortVec_ppm_hhp[i3]){
                    matches.conservativeResize(3, matches.cols()+1);
                    matches.col(matches.cols()-1) << i1,i2,i3;

                    mat1v.push_back( m_ampClass->T3d_makemat_2(i3, i2) );
                    mat2v.push_back( m_ampClass->T3d_makemat_3(i1, i2) );
                }
            }
        }
    }


    int i1; int i2; int i3; int cols = matches.cols();
#pragma omp parallel for num_threads(m_numThreads) private(i1,i2,i3) firstprivate(cols)
    for (int i=0; i<cols; i++){
        i1 = matches(0,i); i2 = matches(1,i); i3 = matches(2,i);

        MatrixX product = 0.5*mat1v[i]*mat2v[i].transpose();

        m_ampClass->T3d_inverse(product, i3, i1);
    }
    
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
    m_ampClass->T3D_remap.clear();

}

void Diagrams::T3e(){

    for (int i1=0; i1<m_intClass->sortVec_ppm_hhp.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_p_p.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_ppm_hhh.size(); i3++){
                if ( m_intClass->sortVec_ppm_hhp[i1] == m_intClass->sortVec_p_p[i2] && m_intClass->sortVec_ppm_hhp[i1] == m_intClass->sortVec_ppm_hhh[i3]){

                    MatrixX mat1 = m_ampClass->T3e_makemat_1(i1, i2);
                    MatrixX mat2 = m_intClass->T3e_makemat(i3, i2);
                    MatrixX product = mat1*mat2.transpose();

                    m_ampClass->T3e_Inverse_temp(product, i1, i3);
                }
            }
        }
    }

    //now use remapped function
    //we now parallellize since we are mapping back T3 amplitudes (demanding)

    Eigen::MatrixXi matches;

    std::vector<MatrixX> mat1v;
    std::vector<MatrixX> mat2v;

    for (int i1=0; i1<m_intClass->sortVec_pp_pp.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pp_hh.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_pppm_hhhp.size(); i3++){
                if ( m_intClass->sortVec_pp_pp[i1] == m_intClass->sortVec_pp_hh[i2] && m_intClass->sortVec_pp_pp[i1] == m_intClass->sortVec_pppm_hhhp[i3]){
                    matches.conservativeResize(3, matches.cols()+1);
                    matches.col(matches.cols()-1) << i1,i2,i3;

                    mat1v.push_back( m_ampClass->T3e_makemat_2(i3, i2) );
                    mat2v.push_back( m_ampClass->T3e_makemat_3(i2, i1) );
                }
            }
        }
    }


    int i1; int i2; int i3; int cols = matches.cols();
#pragma omp parallel for num_threads(m_numThreads) private(i1,i2,i3) firstprivate(cols)
    for (int i=0; i<cols; i++){
        i1 = matches(0,i); i2 = matches(1,i); i3 = matches(2,i);

        MatrixX M1(mat1v[i].rows(), 2*mat1v[i].cols());
        M1 << mat1v[i], -mat1v[i];

        MatrixX M2(2*mat2v[i].rows(), mat2v[i].cols());
        M2 << mat2v[i], -mat2v[i];

        MatrixX product = -0.5*M1*M2;

        m_ampClass->T3e_inverse(product, i3, i1);
    }
    
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
    m_ampClass->T3D_remap.clear();

}

void Diagrams::T5a(){

    Eigen::MatrixXi matches;

    std::vector<MatrixX> mat1v;
    std::vector<MatrixX> mat2v;

    for (int i1=0; i1<m_intClass->sortVec_pm_hp.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pm_ph.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_ppmm_pphh.size(); i3++){
                if ( m_intClass->sortVec_pm_hp[i1] == m_intClass->sortVec_ppmm_pphh[i3] && m_intClass->sortVec_pm_hp[i1] == m_intClass->sortVec_pm_ph[i2]){
                    matches.conservativeResize(3, matches.cols()+1);
                    matches.col(matches.cols()-1) << i1,i2,i3;

                    mat1v.push_back( m_ampClass->T5a_makemat_1(i1, i2) );
                    mat2v.push_back( m_intClass->T5a_makemat(i1, i2) );
                }
            }
        }
    }

    int i1; int i2; int i3; int cols = matches.cols();
#pragma omp parallel for num_threads(m_numThreads) private(i1,i2,i3) firstprivate(cols)
    for (int i=0; i<cols; i++){
        i1 = matches(0,i); i2 = matches(1,i); i3 = matches(2,i);

        MatrixX mat3 = m_ampClass->T5a_makemat_2(i3, i1);

        MatrixX product = mat3*mat2v[i]*mat1v[i].transpose();

        m_ampClass->T5a_inverse(product, i3, i1);
    }
    
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}

void Diagrams::makeT5bIndexMat(){
    for (int i1=0; i1<m_intClass->sortVec_p_h.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_pph.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_pppmm_ppphh.size(); i3++){
                if ( m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_pppmm_ppphh[i3] && m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_pph[i2]){
                    
                    m_ampClass->T3_T5b_indices.push_back( m_ampClass->T5b_makemat_2_I( i3,i1 ) );
                    m_ampClass->T3_T5b_indices_signs.push_back( m_ampClass->T5b_makemat_2_I_signs( i3,i1 ) );
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
                    
                    m_ampClass->T3_T5c_indices.push_back( m_ampClass->T5c_makemat_2_I( i3,i1 ) );
                    m_ampClass->T3_T5c_indices_signs.push_back( m_ampClass->T5c_makemat_2_I_signs( i3,i1 ) );
                }
            }
        }
    }
}

MatrixX Diagrams::contructMatT5b(int index){
    int rows = m_ampClass->T3_T5b_indices[index].rows();
    int cols = m_ampClass->T3_T5b_indices[index].cols();
    MatrixX returnMat;
    returnMat.conservativeResize(rows, cols);
    #pragma omp parallel for num_threads(m_numThreads) firstprivate(rows, cols)
    for (int row=0; row<rows; row++){
        for (int col=0; col<cols; col++){
            returnMat(row,col) = m_ampClass->T3_elements_A[m_ampClass->T3_T5b_indices[index](row, col) ]
                                *(double)m_ampClass->T3_T5b_indices_signs[index](row, col);
        }
    }
    return returnMat;
}

MatrixX Diagrams::contructMatT5c(int index){
    int rows = m_ampClass->T3_T5c_indices[index].rows();
    int cols = m_ampClass->T3_T5c_indices[index].cols();
    MatrixX returnMat;
    returnMat.conservativeResize(rows, cols);
    #pragma omp parallel for num_threads(m_numThreads) firstprivate(rows, cols)
    for (int row=0; row<rows; row++){
        for (int col=0; col<cols; col++){
            returnMat(row,col) = m_ampClass->T3_elements_A[m_ampClass->T3_T5c_indices[index](row, col) ]
                                *(double)m_ampClass->T3_T5c_indices_signs[index](row, col);;
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

    int index = 0;

    MatrixX mat1, mat2, mat3, product, tempMat;
    for (int i1=0; i1<m_intClass->sortVec_p_h.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_pph.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_pppmm_ppphh.size(); i3++){
                if ( m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_pppmm_ppphh[i3] && m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_pph[i2]){

                    mat1 = m_ampClass->T5b_makemat_1(i2, i1);
                    MatrixX M1(2*mat1.rows(), mat1.cols());
                    M1 << mat1, -mat1;

                    mat2 = m_intClass->T5b_makemat(i2, i1);
                    MatrixX M2(2*mat2.rows(), mat2.cols());
                    M2 << mat2, -mat2;
                    
                    tempMat = M2.transpose()*M1;

                    mat3 = contructMatT5b(index);

                    product = -0.5*mat3*tempMat;

                    m_ampClass->T5b_inverse_I(product, index);
                    index++;
                }
            }
        }
    }

    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}

void Diagrams::T5c(){

    int index = 0;

    MatrixX mat1, mat2, mat3, product, tempMat;
    for (int i1=0; i1<m_intClass->sortVec_p_p.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_hhp.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_pppmm_hhhpp.size(); i3++){
                if ( m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_pppmm_hhhpp[i3] && m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_hhp[i2]){

                    mat1 = m_ampClass->T5c_makemat_1(i2, i1);
                    MatrixX M1(2*mat1.rows(), mat1.cols());
                    M1 << mat1, -mat1;

                    mat2 = m_intClass->T5c_makemat(i2, i1);
                    MatrixX M2(2*mat2.rows(), mat2.cols());
                    M2 << mat2, -mat2;

                    tempMat = M1.transpose()*M2;

                    mat3 = contructMatT5c(index);
                    product = -0.5*mat3*tempMat;


                    m_ampClass->T5c_inverse_I(product, index);
                    index++;
                }
            }
        }
    }

    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}

void Diagrams::T5d(){

    Eigen::MatrixXi matches;

    std::vector<MatrixX> mat1v;
    std::vector<MatrixX> mat2v;

    for (int i1=0; i1<m_intClass->sortVec_p_p.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_hhp.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_ppm_pph.size(); i3++){
                if ( m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_pph[i3] && m_intClass->sortVec_p_p[i1] == m_intClass->sortVec_ppm_hhp[i2]){
                    matches.conservativeResize(3, matches.cols()+1);
                    matches.col(matches.cols()-1) << i1,i2,i3;

                    mat1v.push_back( m_ampClass->T5d_makemat_1(i2, i1) );
                    mat2v.push_back( m_intClass->T5d_makemat(i2, i1) );
                }
            }
        }
    }

    int i1; int i2; int i3; int cols = matches.cols();
#pragma omp parallel for num_threads(m_numThreads) private(i1,i2,i3) firstprivate(cols)
    for (int i=0; i<cols; i++){
        i1 = matches(0,i); i2 = matches(1,i); i3 = matches(2,i);

        MatrixX M2(2*mat2v[i].rows(), mat2v[i].cols());
        M2 << mat2v[i], -mat2v[i];

        MatrixX mat3 = m_ampClass->T5d_makemat_2(i2, i3);
        MatrixX M3(2*mat3.rows(), mat3.cols());
        M3 << mat3, -mat3;

        MatrixX product = -0.5*mat1v[i]*M2.transpose()*M3;

        m_ampClass->T5d_inverse(product, i2, i3);
    }

    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}

void Diagrams::T5e(){

    Eigen::MatrixXi matches;

    std::vector<MatrixX> mat1v;
    std::vector<MatrixX> mat2v;

    for (int i1=0; i1<m_intClass->sortVec_p_h.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_ppm_hhp.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_ppm_pph.size(); i3++){
                if ( m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_pph[i3] && m_intClass->sortVec_p_h[i1] == m_intClass->sortVec_ppm_hhp[i2]){
                    matches.conservativeResize(3, matches.cols()+1);
                    matches.col(matches.cols()-1) << i1,i2,i3;

                    mat1v.push_back( m_ampClass->T5e_makemat_1(i3, i1) );
                    mat2v.push_back( m_intClass->T5e_makemat(i3, i1) );
                }
            }
        }
    }

    int i1; int i2; int i3; int cols = matches.cols();
    #pragma omp parallel for num_threads(m_numThreads) private(i1,i2,i3) firstprivate(cols)
    for (int i=0; i<cols; i++){
        i1 = matches(0,i); i2 = matches(1,i); i3 = matches(2,i);

        MatrixX M2(2*mat2v[i].rows(), mat2v[i].cols());
        M2 << mat2v[i], -mat2v[i];

        MatrixX mat3 = m_ampClass->T5e_makemat_2(i2, i3);
        MatrixX M3(mat3.rows(), 2*mat3.cols());
        M3 << mat3, -mat3;
        
        MatrixX product = -0.5*M3*M2*mat1v[i].transpose();

        m_ampClass->T5e_inverse(product, i2, i3);
    }
    
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}

void Diagrams::T5f(){

    Eigen::MatrixXi matches;

    std::vector<MatrixX> mat1v;
    std::vector<MatrixX> mat2v;

    for (int i1=0; i1<m_intClass->sortVec_pp_hh.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pp_pp.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_pppm_ppph.size(); i3++){
                if ( m_intClass->sortVec_pp_hh[i1] == m_intClass->sortVec_pp_pp[i2] && m_intClass->sortVec_pp_hh[i1] == m_intClass->sortVec_pppm_ppph[i3]){
                    matches.conservativeResize(3, matches.cols()+1);
                    matches.col(matches.cols()-1) << i1,i2,i3;

                    mat1v.push_back( m_ampClass->T5f_makemat_1(i1, i2) );
                    mat2v.push_back( m_intClass->T5f_makemat(i1, i2) );
                }
            }
        }
    }

    int i1; int i2; int i3; int cols = matches.cols();
    #pragma omp parallel for num_threads(m_numThreads) private(i1,i2,i3) firstprivate(cols)
    for (int i=0; i<cols; i++){
        i1 = matches(0,i); i2 = matches(1,i); i3 = matches(2,i);

        MatrixX M1(mat1v[i].rows(), 2*mat1v[i].cols());
        M1 << mat1v[i], -mat1v[i];

        MatrixX M2(mat2v[i].rows(), 2*mat2v[i].cols());
        M2 << mat2v[i], -mat2v[i];

        MatrixX tempMat = M1*M2.transpose();

        MatrixX MT(tempMat.rows(), 2*tempMat.cols());
        MT << tempMat, -tempMat;

        MatrixX mat3 = m_ampClass->T5f_makemat_2(i3, i1);
        MatrixX M3(mat3.rows(), 2*mat3.cols());
        M3 << mat3, -mat3;

        MatrixX product = 0.25*M3*MT.transpose();

        m_ampClass->T5f_inverse(product, i3, i1);
    }
    
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}

void Diagrams::T5g(){

    Eigen::MatrixXi matches;

    std::vector<MatrixX> mat1v;
    std::vector<MatrixX> mat2v;

    for (int i1=0; i1<m_intClass->sortVec_pp_hh.size(); i1++){
        for (int i2=0; i2<m_intClass->sortVec_pp_pp.size(); i2++){
            for (int i3=0; i3<m_intClass->sortVec_pppm_hhhp.size(); i3++){
                if ( m_intClass->sortVec_pp_hh[i1] == m_intClass->sortVec_pp_pp[i2] && m_intClass->sortVec_pp_hh[i1] == m_intClass->sortVec_pppm_hhhp[i3]){
                    matches.conservativeResize(3, matches.cols()+1);
                    matches.col(matches.cols()-1) << i1,i2,i3;

                    mat1v.push_back( m_ampClass->T5g_makemat_1(i1, i2) );
                    mat2v.push_back( m_intClass->T5g_makemat(i1, i2) );
                }
            }
        }
    }

    int i1; int i2; int i3; int cols = matches.cols();
    #pragma omp parallel for num_threads(m_numThreads) private(i1,i2,i3) firstprivate(cols)
    for (int i=0; i<cols; i++){
        i1 = matches(0,i); i2 = matches(1,i); i3 = matches(2,i);

        MatrixX M1(2*mat1v[i].rows(), mat1v[i].cols());
        M1 << mat1v[i], -mat1v[i];

        MatrixX M2(2*mat2v[i].rows(), mat2v[i].cols());
        M2 << mat2v[i], -mat2v[i];

        MatrixX tempMat = M2.transpose()*M1;

        MatrixX MT(2*tempMat.rows(), tempMat.cols());
        MT << tempMat, -tempMat;

        MatrixX mat3 = m_ampClass->T5g_makemat_2(i3, i2);
        MatrixX M3(mat3.rows(), 2*mat3.cols());
        M3 << mat3, -mat3;

        MatrixX product = 0.25*M3*MT;

        m_ampClass->T5g_inverse(product, i3, i2);
    }
    
    std::fill(m_ampClass->T3_elements_A_temp.begin(), m_ampClass->T3_elements_A_temp.end(), 0); //reset T3 temp
}
