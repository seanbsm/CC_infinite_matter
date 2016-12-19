#include "master.h"
#include "Systems/heg.h"
#include "makeampmat.h"
#include "makeintmat.h"
#include "diagrams.h"

#include <iostream>

using namespace std;

void Master::setSize(int Nh, int Nb){
    m_Nh = Nh;
    m_Nb = Nb;
}

void Master::setSystem(class System* system){
    m_system = system;
}

double Master::Iterator(double eps, double conFac){

    MakeIntMat* Interaction = new MakeIntMat;
    MakeAmpMat* Amplituder  = new MakeAmpMat;
    Diagrams*   diagrams    = new Diagrams;
    cout << m_Ns << endl;
    Interaction->makeBlockMat(m_system, m_Nh, m_Ns);

    //would be better to simply let Amplituder inherit master?
    Amplituder->setIntClass(Interaction);
    Amplituder->setSystem(m_system);
    Amplituder->makeBlockMat();
    Amplituder->makeDenomMat();

    double ECCD     = 0;
    double ECCD_old = 0;
    for (int h = 0; h<Interaction->Vhhpp.size(); h++){
        cout << "hey" << endl;
        cout << Amplituder->Amplitudes[h].rows() << " " << Amplituder->Amplitudes[h].cols() << endl;
        cout << Amplituder->denomMat[h].rows() << " " << Amplituder->denomMat[h].cols() << endl;
        Eigen::MatrixXf temp = Amplituder->Amplitudes[h]*Amplituder->denomMat[h];
        ECCD_old += 0.25*((Interaction->Vhhpp[h].transpose())*(temp)).trace();
    }
    std::cout << "MBPT2: " << ECCD_old << std::endl;

    while (conFac > eps){
        ECCD = 0;
        for (int hh = 0; hh<Interaction->Vhhpp.size(); hh++){
            for (int pp = 0; pp<Interaction->Vpppp.size(); pp++){
                if (Interaction->Vhhpp_i[hh] == Interaction->Vpppp_i[pp]){
                    Amplituder->Amplitudes[hh] = Interaction->Vhhpp[hh] + diagrams->D3(Amplituder->Amplitudes[hh], Interaction->Vpppp[pp]);
                }
            }
            ECCD += 0.25*((Interaction->Vhhpp[hh].transpose())*(Amplituder->Amplitudes[hh])).trace();
        }
        cout << ECCD << endl;
        conFac = abs(ECCD - ECCD_old);
        ECCD_old = ECCD;

        //ECCD = 0; too good to delete
    }
    return ECCD;

}
