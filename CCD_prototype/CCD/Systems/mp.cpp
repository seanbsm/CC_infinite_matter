#include "mp.h"
#include <iostream>

using namespace std;

MP::MP(Master* master, double m, double L3, double L2, double L1) : System(master) /* Minnesota Potential */
{
    m_Nh = master->m_Nh;
    m_Nb = master->m_Nb;
    m_dk = 2*m_Nb + 1;
    m_master = master;

    m_m  = m;
    m_L3 = L3;
    m_L2 = L2;
    m_L1 = L1;
    V_0R_fac = V_0R*pow(pi/kappa_R, 1.5)/m_L3;
    V_0S_fac = -V_0S*pow(pi/kappa_S, 1.5)/m_L3; //the minus signs come from morten's code
    V_0T_fac = -V_0T*pow(pi/kappa_T, 1.5)/m_L3;
    piOverL  = pi/m_L1;
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
                        m_states.conservativeResize(Eigen::NoChange, m_states.cols()+2);
                        m_states.col(m_states.cols()-2) << n2,nx,ny,nz, 1, 1;
                        m_states.col(m_states.cols()-1) << n2,nx,ny,nz,-1, 1;
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
int MP::kUnique1(int k, int s1){
    Eigen::Matrix<int, 5, 1> kk;
    kk << m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k), m_states(5,k);
    Eigen::VectorXi mom = s1*kk;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk + mom(4)*m_dk*m_dk*m_dk*m_dk;
    return kuni;
}

int MP::kUnique2(int k, int p, int s1, int s2){
    Eigen::Matrix<int, 5, 1> kk;
    Eigen::Matrix<int, 5, 1> kp;
    kk << m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k), m_states(5,k);
    kp << m_states(1,p), m_states(2,p), m_states(3,p), m_states(4,p), m_states(5,p);
    Eigen::VectorXi mom = s1*kk+s2*kp;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk + mom(4)*m_dk*m_dk*m_dk*m_dk;
    return kuni;
}

int MP::kUnique3(int k, int p, int q, int s1, int s2, int s3){
    Eigen::Matrix<int, 5, 1> kk;
    Eigen::Matrix<int, 5, 1> kp;
    Eigen::Matrix<int, 5, 1> kq;
    kk << m_states(1,k), m_states(2,k), m_states(3,k), m_states(4,k), m_states(5,k);
    kp << m_states(1,p), m_states(2,p), m_states(3,p), m_states(4,p), m_states(5,p);
    kq << m_states(1,q), m_states(2,q), m_states(3,q), m_states(4,q), m_states(5,q);
    Eigen::VectorXi mom = s1*kk+s2*kp+s3*kq;

    int kuni = mom(0) + mom(1)*m_dk + mom(2)*m_dk*m_dk + mom(3)*m_dk*m_dk*m_dk + mom(4)*m_dk*m_dk*m_dk*m_dk;
    return kuni;
}

double MP::f(int p){
    double returnVal = h0(p);
    for (int i=0; i<m_Nh; i++){
        returnVal += assym_single(p, i);
    };
    return returnVal;
}

double MP::h0(int p){
    double energy = m_states(0,p);
    return energy*2*pi*pi*m_hbarc*m_hbarc/(m_m*m_L2);
}

/*
double MP::assym(int p, int q, int r, int s){
    Eigen::Vector3i kp( m_states(1,p), m_states(2,p), m_states(3,p) );
    Eigen::Vector3i kq( m_states(1,q), m_states(2,q), m_states(3,q) );
    Eigen::Vector3i kr( m_states(1,r), m_states(2,r), m_states(3,r) );
    Eigen::Vector3i ks( m_states(1,s), m_states(2,s), m_states(3,s) );
    int sp = m_states(4,p); int tp = m_states(5,p);
    int sq = m_states(4,q); int tq = m_states(5,q);
    int sr = m_states(4,r); int tr = m_states(5,r);
    int ss = m_states(4,s); int ts = m_states(5,s);

    //these tests should be already performed through k_unique
    if ( vecDelta(kp+kq, kr+ks) == 0){ return 0;}   //momentum conservation
    if ( sp+sq != sr+ss){ return 0; }               //spin conservation
    if ( tp+tq != tr+ts){ return 0; }               //isospin conservation

    //I'm not certain if this indexing is correct
    //I've dropped the factor of 0.5 since I think it might be a misprint
    //piOverL is needed for momentum, not just quantum numbers
    double q2_dir = piOverL*((kp-kq-kr+ks).adjoint()*(kp-kq-kr+ks))(0,0);
    double q2_ex  = piOverL*((kp-kq+kr-ks).adjoint()*(kp-kq+kr-ks))(0,0);

    double V_1R = V_0R_fac*exp(-q2_dir/(4*kappa_R));
    double V_1T = V_0T_fac*exp(-q2_dir/(4*kappa_T));
    double V_1S = V_0S_fac*exp(-q2_dir/(4*kappa_S));

    double V_2R = V_0R_fac*exp(-q2_ex /(4*kappa_R));
    double V_2T = V_0T_fac*exp(-q2_ex /(4*kappa_T));
    double V_2S = V_0S_fac*exp(-q2_ex /(4*kappa_S));

    double returnVal = 0;

    //copying morten's code, didn't help much
    bool alpha1 = (sp==ss && sq==sr);                       bool beta1 = (sp==sr && sq==ss);
    bool alpha2 = (tp==ts && tq==tr);                       bool beta2 = (tp==tr && tq==ts);
    bool alpha3 = (sp==sr && sq==ss && tp==tr && tq==ts);   bool beta3 = (sp==ss && sq==sr && tp==ts && tq==tr);
    bool alpha4 = alpha1*(tp==tr && tq==ts);                bool beta4 = beta1*(tp==ts && tq==tr);
    bool alpha5 = alpha1*alpha2;                            bool beta5 = beta1*beta2;
    bool alpha6 = alpha2*(sp==sr && sq==ss);                bool beta6 = beta2*(sp==ss && sq==sr);

    returnVal += V_1R*(alpha3-alpha5) + 0.5*V_1T*(alpha3+alpha4-alpha5-alpha6) + 0.5*V_1S*(alpha3-alpha4-alpha5+alpha6);
    returnVal -= V_2R*(beta3-beta5) + 0.5*V_2T*(beta3+beta4-beta5-beta6) + 0.5*V_2S*(beta3-beta4-beta5+beta6);

    return 0.5*returnVal;
}*/

double MP::assym(int p, int q, int r, int s){
    Eigen::Vector3i kp( m_states(1,p), m_states(2,p), m_states(3,p) );
    Eigen::Vector3i kq( m_states(1,q), m_states(2,q), m_states(3,q) );
    Eigen::Vector3i kr( m_states(1,r), m_states(2,r), m_states(3,r) );
    Eigen::Vector3i ks( m_states(1,s), m_states(2,s), m_states(3,s) );
    int sp = m_states(4,p); int tp = m_states(5,p);
    int sq = m_states(4,q); int tq = m_states(5,q);
    int sr = m_states(4,r); int tr = m_states(5,r);
    int ss = m_states(4,s); int ts = m_states(5,s);

    //these tests should be already performed through k_unique
    if ( vecDelta(kp+kq, kr+ks) == 0){ return 0;}   //momentum conservation
    if ( sp+sq != sr+ss){ return 0; }               //spin conservation
    if ( tp+tq != tr+ts){ return 0; }               //isospin conservation

    double q2_dir = 0.25*piOverL*piOverL*(kp-kq-kr+ks).squaredNorm();   // <pq||rs> (direct)
    double q2_ex = 0.25*piOverL*piOverL*(kp-kq+kr-ks).squaredNorm();    // <pq||sr> (exchange)

   // cout << q2_ex << endl;

    double VR_dir = V_0R_fac*exp(-q2_dir/(4*kappa_R));
    double VT_dir = V_0T_fac*exp(-q2_dir/(4*kappa_T));
    double VS_dir = V_0S_fac*exp(-q2_dir/(4*kappa_S));

    double VR_ex = V_0R_fac*exp(-q2_ex/(4*kappa_R));
    double VT_ex = V_0T_fac*exp(-q2_ex/(4*kappa_T));
    double VS_ex = V_0S_fac*exp(-q2_ex/(4*kappa_S));


    bool Ps_dir = (sp==ss)*(sq==sr);    //exchange spins
    bool Cs_dir = (sp==sr)*(sq==ss);    //direct spins
    bool Pt_dir = (tp==ts)*(tq==tr);    //exchange isospins
    bool Ct_dir = (tp==tr)*(tq==ts);    //direct isospins

    bool Ps_ex = (sp==sr)*(sq==ss);    //exchange spins
    bool Cs_ex = (sp==ss)*(sq==sr);    //direct spins
    bool Pt_ex = (tp==tr)*(tq==ts);    //exchange isospins
    bool Ct_ex = (tp==ts)*(tq==tr);    //direct isospins

    double returnVal = 0;

    returnVal += (VR_dir + 0.5*VT_dir + 0.5*VS_dir)*Cs_dir*Ct_dir
               + 0.5*(VT_dir - VS_dir)*Ps_dir*Ct_dir
               - (VR_dir + 0.5*VT_dir + 0.5*VS_dir)*Ps_dir*Pt_dir
               - 0.5*(VT_dir - VS_dir)*Cs_dir*Pt_dir;

    returnVal += -(VR_ex + 0.5*VT_ex + 0.5*VS_ex)*Cs_ex*Ct_ex
               - 0.5*(VT_ex - VS_ex)*Ps_ex*Ct_ex
               + (VR_ex + 0.5*VT_ex + 0.5*VS_ex)*Ps_ex*Pt_ex
               + 0.5*(VT_ex - VS_ex)*Cs_ex*Pt_ex;

    return 0.5*returnVal;
}

double MP::assym_single(int p, int q){
    return assym(p,q,p,q);
}

//can do a test by copying assym, and set int s = p, int r = q
//worked for Nh=Nb=2 and Nh=14, Nb=3
/*double MP::assym_single(int p, int q){
    Eigen::Vector3i kp( m_states(1,p), m_states(2,p), m_states(3,p) ); //r = p
    Eigen::Vector3i kq( m_states(1,q), m_states(2,q), m_states(3,q) ); //s = q
    int sp = m_states(4,p); int tp = m_states(5,p);
    int sq = m_states(4,q); int tq = m_states(5,q);

    //I'm not certain if this indexing is correct
    //I've dropped the factor of 0.5 since I think it might be a misprint
    //piOverL is needed for momentum, not just quantum numbers
    //double q2_dir = 0;
    double q2_ex  = piOverL*((kp-kq).adjoint()*(kp-kq))(0,0);

    double V_1R = V_0R_fac;
    double V_1T = V_0T_fac;
    double V_1S = V_0S_fac;

    double V_2R = V_0R_fac*exp(-q2_ex /kappa_R);
    double V_2T = V_0T_fac*exp(-q2_ex /kappa_T);
    double V_2S = V_0S_fac*exp(-q2_ex /kappa_S);

    double returnVal = 0;

    //copying morten's code, didn't help much
    bool alpha1 = (sp==sq);       bool beta1 = true;
    bool alpha2 = (tp==tq);       bool beta2 = true;
    bool alpha3 = true;           bool beta3 = (sp==sq && tp==tq);
    bool alpha4 = alpha1;         bool beta4 = beta1*(tp==tq);
    bool alpha5 = alpha1*alpha2;  bool beta5 = beta1*beta2;
    bool alpha6 = alpha2;         bool beta6 = beta2*(sp==sq);

    returnVal += V_1R*(alpha3-alpha5) + 0.5*V_1T*(alpha3+alpha4-alpha5-alpha6) + 0.5*V_1S*(alpha3-alpha4-alpha5+alpha6);
    returnVal -= V_2R*(beta3-beta5) + 0.5*V_2T*(beta3+beta4-beta5-beta6) + 0.5*V_2S*(beta3-beta4-beta5+beta6);

    return 0.5*returnVal;
}*/

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
