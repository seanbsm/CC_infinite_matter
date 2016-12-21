#include "master.h"
#include "Systems/heg.h"
#include "makeampmat.h"
#include "makeintmat.h"
#include "diagrams.h"

#include <iostream>
 #include <iomanip> //needed for std::setprecision

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

    //would be better to simply let "Amplituder" and "diagrams" inherit master?
    diagrams->setAmpClass(Amplituder);
    diagrams->setIntClass(Interaction);
    Amplituder->setIntClass(Interaction);
    Amplituder->setSystem(m_system);
    Amplituder->makeBlockMat();
    Amplituder->makeDenomMat();

    double ECCD     = 0;
    double ECCD_old = 0;
    for (int h = 0; h<Interaction->Vhhpp.size(); h++){
        //Using array<->matrix conversion costs no cpu time in Eigen, so this is fine
        Eigen::MatrixXf temp = Amplituder->Amplitudes[h].array()*Amplituder->denomMat[h].array();
        ECCD_old += ((Interaction->Vhhpp[h].transpose())*(temp)).trace();
    }
    std::cout << "MBPT2: " << std::setprecision (12) << ECCD_old << std::endl;

    int counter = 0;
    while (conFac > eps && counter < 1e3){
        ECCD = 0;
        for (int hh = 0; hh<Interaction->Vhhpp.size(); hh++){

            Amplituder->Amplitudes[hh] = ( Interaction->Vhhpp[hh]
                                           + diagrams->La(hh)
                                           + diagrams->Lb(hh) ).array()
                                         *Amplituder->denomMat[hh].array();
            /*
            for (int pp = 0; pp<Interaction->Vpppp.size(); pp++){ //all Q diagrams fall within this loop
                if (Interaction->Vhhpp_i[hh] == Interaction->Vpppp_i[pp]){
                    Amplituder->Amplitudes[hh] = ( Interaction->Vhhpp[hh] + diagrams->La(hh, pp) ).array()
                                                 *Amplituder->denomMat[hh].array();
                }
            }
            */

            /*
            for (int hp = 0; hp<Interaction->Vhphp.size(); hp++){ //this loop covers Lc
                if (Interaction->Vhhpp_i[hh] == Interaction->Vhphp_i[hp]){
                    Amplituder->Amplitudes[hh] =
                }
            }
            */

            ECCD += 0.25*((Interaction->Vhhpp[hh].transpose())*(Amplituder->Amplitudes[hh])).trace();
        }
        cout << std::setprecision (12) << ECCD << endl;
        conFac = abs(ECCD - ECCD_old);
        ECCD_old = ECCD;
        counter += 1;
        //ECCD = 0; too good to delete; you don't want to know how long i used on this
    }
    return ECCD;

}
