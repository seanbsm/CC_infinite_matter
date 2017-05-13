#include "chiral.h"
#include <iostream>

using namespace std;

extern "C" {
    // get real and imaginary components of matrix element.
    void chipot_f90_wrapper_(double *matel_real, double *matel_im,
                            int *Nparticles, double *rho,
                            int *ps, int *pt, int *ppx, int *ppy, int* ppz,
                            int *qs, int *qt, int *qpx, int *qpy, int *qpz,
                            int *rs, int *rt, int *rpx, int *rpy, int *rpz,
                            int *ss, int *st, int *spx, int *spy, int *spz);
}

CHIRAL::CHIRAL(class Master* master, double m, double L3, double L2, double L1) : System(master) /* Chiral Potential */
{
    m_Nh = master->m_Nh;
    m_Nb = master->m_Nb;
    m_dk = 2*m_Nb + 1;
    m_master = master;

    m_m  = m;
    m_L3 = L3;
    m_L2 = L2;
    m_L1 = L1;
    makeStateSpace();
}

void CHIRAL::makeStateSpace(){
    m_states.conservativeResize(m_states.rows()+5, Eigen::NoChange);    //sets the format of states
    //start for-loops
    for (int n2=0; n2<m_Nb+1; n2++){
        for (int nx=-n2; nx<n2+1; nx++){
            for (int ny=-n2; ny<n2+1; ny++){
                for (int nz=-n2; nz<n2+1; nz++){
                    if (nx*nx + ny*ny + nz*nz == n2){
                        m_states.conservativeResize(Eigen::NoChange, m_states.cols()+2);
                        m_states.col(m_states.cols()-2) << n2,nx,ny,nz, 1;
                        m_states.col(m_states.cols()-1) << n2,nx,ny,nz,-1;
                    } //end of nx^2 + ny^2 + nz^2 == n^2
                }//end of nz-loop
            }//end of ny-loop
        }//end of nx-loop
    }//end of n2-loop

    m_master->m_Ns = m_states.cols();
    below_fermi = Eigen::VectorXi::LinSpaced(m_Nh,0,m_Nh);
    above_fermi = Eigen::VectorXi::LinSpaced(m_Ns,m_Nh,m_Ns);
}

//I think using eigen here is a bit over-the-top for such a function, but whatevs~
int CHIRAL::kUnique1(int k, int s1){
    Eigen::Vector4i kk( m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k) );
    Eigen::VectorXi mom = s1*kk;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk;
    return kuni;
}

int CHIRAL::kUnique2(int k, int p, int s1, int s2){
    Eigen::Vector4i kk( m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k) );
    Eigen::Vector4i kp( m_states(1,p), m_states(2,p), m_states(3,p), m_states(4,p) );
    Eigen::VectorXi mom = s1*kk + s2*kp;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk;
    return kuni;
}

int CHIRAL::kUnique3(int k, int p, int q, int s1, int s2, int s3){
    Eigen::Vector4i kk( m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k) );
    Eigen::Vector4i kp( m_states(1,p), m_states(2,p), m_states(3,p), m_states(4,p) );
    Eigen::Vector4i kq( m_states(1,q), m_states(2,q), m_states(3,q), m_states(4,q) );
    Eigen::VectorXi mom = s1*kk + s2*kp + s3*kq;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk;
    return kuni;
}

int CHIRAL::kUnique4(int k, int p, int q, int s, int s1, int s2, int s3, int s4){
    Eigen::Vector4i kk( m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k) );
    Eigen::Vector4i kp( m_states(1,p), m_states(2,p), m_states(3,p), m_states(4,p) );
    Eigen::Vector4i kq( m_states(1,q), m_states(2,q), m_states(3,q), m_states(4,q) );
    Eigen::Vector4i ks( m_states(1,s), m_states(2,s), m_states(3,s), m_states(4,s) );
    Eigen::VectorXi mom = s1*kk + s2*kp + s3*kq + s4*ks;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk;
    return kuni;
}

int CHIRAL::kUnique5(int k, int p, int q, int s, int t, int s1, int s2, int s3, int s4, int s5){
    Eigen::Vector4i kk( m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k) );
    Eigen::Vector4i kp( m_states(1,p), m_states(2,p), m_states(3,p), m_states(4,p) );
    Eigen::Vector4i kq( m_states(1,q), m_states(2,q), m_states(3,q), m_states(4,q) );
    Eigen::Vector4i ks( m_states(1,s), m_states(2,s), m_states(3,s), m_states(4,s) );
    Eigen::Vector4i kt( m_states(1,t), m_states(2,t), m_states(3,t), m_states(4,t) );
    Eigen::VectorXi mom = s1*kk + s2*kp + s3*kq + s4*ks + s5*kt;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk;
    return kuni;
}

System::variable_type CHIRAL::f(int p){
    variable_type returnVal = h0(p);
    for (int i=0; i<m_Nh; i++){
        returnVal += assym_single(p, i);
    };
    return returnVal;
}

System::variable_type CHIRAL::h0(int p){
    variable_type energy = (variable_type)m_states(0,p);
    return energy*2.*pi*pi*m_hbarc*m_hbarc/(m_m*m_L2);
}

System::variable_type CHIRAL::assym(int p, int q, int r, int s){
    /*Eigen::Vector3i kp( m_states(1,p), m_states(2,p), m_states(3,p) );
    Eigen::Vector3i kq( m_states(1,q), m_states(2,q), m_states(3,q) );
    Eigen::Vector3i kr( m_states(1,r), m_states(2,r), m_states(3,r) );
    Eigen::Vector3i ks( m_states(1,s), m_states(2,s), m_states(3,s) );
    int sp = m_states(4,p);
    int sq = m_states(4,q);
    int sr = m_states(4,r);
    int ss = m_states(4,s);*/

    double matel_real, matel_im, rho;
    int isospin = 1;
    complex<double> result;// double temp;

    rho = m_Nh/m_L3;

    //Morten thinks +1=neutrons in the fortran code
    chipot_f90_wrapper_(&matel_real, &matel_im,
                        &m_Nh, &rho,
                        &m_states(4,p), &isospin, &m_states(1,p), &m_states(2,p), &m_states(3,p),
                        &m_states(4,q), &isospin, &m_states(1,q), &m_states(2,q), &m_states(3,q),
                        &m_states(4,r), &isospin, &m_states(1,r), &m_states(2,r), &m_states(3,r),
                        &m_states(4,s), &isospin, &m_states(1,s), &m_states(2,s), &m_states(3,s));

    //function as in wrapper
    /*chipot_f90_wrapper_(&matel_real, &matel_im,
                        &Nparticles, &rho,
                        &ps, &pt, &ppx, &ppy, &ppz,
                        &qs, &qt, &qpx, &qpy, &qpz,
                        &rs, &rt, &rpx, &rpy, &rpz,
                        &ss, &st, &spx, &spy, &spz);*/
    result.real(matel_real);
    result.imag(matel_im);

    //temp = real(result);
    return result;
}

System::variable_type CHIRAL::assym_single(int p, int q){
    return assym(p,q,p,q);
}

bool CHIRAL::vecDelta(Eigen::VectorXi v1, Eigen::VectorXi v2){
    int dim1 = v1.rows();
    int dim2 = v2.rows();
    if (dim1 != dim2){
        cout << "dimensional error" << endl;
        return 0;
    }
    else{
        bool returnInt = 1;
        for (int i =0; i<dim1; i++){
            returnInt *= ( v1(i)==v2(i) );
        }
        return returnInt;
    }
}
