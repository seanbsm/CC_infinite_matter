#include "mp.h"
#include <iostream>

using namespace std;

MP::MP(Master* master, double m, double L3, double L2, double L1) : System(master) /* Minnesota Potential */
{
    m_Nh = master->m_Nh;
    m_Nb = master->m_Nb;
    m_master = master;

    m_m  = m;
    m_L3 = L3;
    m_L2 = L2;
    m_L1 = L1;
    makeStateSpace();
}

void MP::makeStateSpace(){
    m_states.conservativeResize(m_states.rows()+6, Eigen::NoChange);    //sets the format of states
    //start for-loops
    for (int n2=0; n2<m_Nb+1; n2++){
        for (int nx=-n2; nx<n2+1; nx++){
            for (int ny=-n2; ny<n2+1; ny++){
                for (int nz=-n2; nz<n2+1; nz++){
                    if (nx*nx + ny*ny + nz*nz == n2){
                        m_states.conservativeResize(Eigen::NoChange, m_states.cols()+4);
                        m_states.col(m_states.cols()-4) << n2,nx,ny,nz, 1, 1;
                        m_states.col(m_states.cols()-3) << n2,nx,ny,nz, 1,-1;
                        m_states.col(m_states.cols()-2) << n2,nx,ny,nz,-1, 1;
                        m_states.col(m_states.cols()-1) << n2,nx,ny,nz,-1,-1;
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
int MP::kUnique2(int k, int p){
    Eigen::Matrix<int, 5, 1> kk;
    Eigen::Matrix<int, 5, 1> kp;
    kk << m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k), m_states(5,k);
    kp << m_states(1,p), m_states(2,p), m_states(3,p), m_states(4,p), m_states(5,k);
    Eigen::VectorXi mom = kk+kp;

    int val = 0;
    for (int i = 0; i<mom.rows();i++){
        if (val < mom(i)){
            val = mom(i);
        }
    }

    int dk = 2*val + 1;
    int kuni = mom(0) + mom(1)*dk + mom(2)*dk*dk + mom(3)*dk*dk*dk + mom(4)*dk*dk*dk*dk;
    return kuni;
}

double MP::f(int p){
    double returnVal = h0(p);
    for (int i=0; i<below_fermi.rows(); i++){
        returnVal += assym(p,below_fermi(i),p,below_fermi(i));
    };
    return returnVal;
}

double MP::h0(int p){
    double energy = m_states(p,0);
    return energy*2*pi*pi/(m_m*m_L2);
}

double MP::assym(int p, int q, int r, int s){

    Eigen::Vector3i kp( m_states(1,p), m_states(2,p), m_states(3,p) );
    Eigen::Vector3i kq( m_states(1,q), m_states(2,q), m_states(3,q) );
    Eigen::Vector3i kr( m_states(1,r), m_states(2,r), m_states(3,r) );
    Eigen::Vector3i ks( m_states(1,s), m_states(2,s), m_states(3,s) );
    int sp = m_states(4,p); int tp = m_states(5,p);
    int sq = m_states(4,q); int tq = m_states(5,q);
    int sr = m_states(4,r); int tr = m_states(5,r);
    int ss = m_states(4,s); int ts = m_states(5,s);

    //these test should be already performed through k_unique
    if ( vecDelta(kp+kq, kr+ks) == 0){ return 0;}   //momentum conservation
    if ( sp+sq != sr+ss){ return 0; }               //spin conservation
    if ( tp+tq != tr+ts){ return 0; }               //isospin conservation

    //I'm not certain if this indexing is correct
    double q2_dir = 0.5*((kp-kq-kr+ks).adjoint()*(kp-kq-kr+ks))(0,0);
    double q2_ex  = 0.5*((kp-kq+kr-ks).adjoint()*(kp-kq+kr-ks))(0,0);

    //cout << q2_dir << endl;

    double V_1R = V_0R_fac*exp(-q2_dir/(4*kappa_R))/m_L3;
    double V_1T = V_0T_fac*exp(-q2_dir/(4*kappa_T))/m_L3;
    double V_1S = V_0S_fac*exp(-q2_dir/(4*kappa_S))/m_L3;

    double V_2R = V_0R_fac*exp(-q2_ex /(4*kappa_R))/m_L3;
    double V_2T = V_0T_fac*exp(-q2_ex /(4*kappa_T))/m_L3;
    double V_2S = V_0S_fac*exp(-q2_ex /(4*kappa_S))/m_L3;

    //cout << a << endl;
    double returnVal = 0;
    //I'm really unsure about the spin-tests here, might be horribly wrong
    returnVal += ( V_1R + 0.5*(sp==ss)*(sq==sr)*V_1T + 0.5*(sp==ss)*(sq==sr)*V_1S )*(tp==ts)*(tq==tr);
    returnVal -= ( V_2R + 0.5*(sp==sr)*(sq==ss)*V_2T + 0.5*(sp==sr)*(sq==ss)*V_2S )*(tp==tr)*(tq==ts);

    //cout << returnVal << endl;
    return 0.5*returnVal;
}

bool MP::vecDelta(Eigen::VectorXi v1, Eigen::VectorXi v2){
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
