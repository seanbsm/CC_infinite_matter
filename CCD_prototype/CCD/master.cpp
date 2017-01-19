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
    Amplituder->setElements();
    Amplituder->makeDenomMat();

    double ECCD     = 0;
    double ECCD_old = 0;
    //for (int h = 0; h<Interaction->sortVec_hh.size(); h++){
    //std::cout << Interaction->Vhhpp_i.size() << " " << Interaction->numOfKu << std::endl;
    for (int h = 0; h<Interaction->numOfKu; h++){
        //Using array<->matrix conversion costs no cpu time in Eigen, so this is fine
        Eigen::MatrixXd Vhhpp = Interaction->make2x2Block(Interaction->Vhhpp_i[h],0,0,1,1);
        Eigen::MatrixXd temp = Vhhpp.array()*Amplituder->denomMat[h].array();
        ECCD_old += ((Vhhpp.transpose())*(temp)).trace();
    }

    //check whether or not to multiply by 0.25 for MBPT2
    std::cout << "MBPT2: " << std::setprecision (12) << 0.25*ECCD_old << std::endl;

    int counter = 0;
    //std::cout << Interaction->numOfKu << " " << Interaction->Vhhpp_i.size() << std::endl;

    /*
    //test to check inverse block function
    Eigen::MatrixXd temp1 = Amplituder->make2x2Block(Interaction->Vhhpp_i[3],0,0,1,1, Amplituder->T_elements);
    cout << temp1 << endl;
    Amplituder->make2x2Block_inverse(temp1, Interaction->Vhhpp_i[3], 0,0,1,1, Amplituder->T_elements);
    temp1 = Amplituder->make2x2Block(Interaction->Vhhpp_i[3],0,0,1,1, Amplituder->T_elements);
    cout << temp1 << endl;
    */

    while (/*conFac > eps &&*/ counter < 2e1){
        ECCD = 0;
        //Amplituder->T_elements = Amplituder->T_elements_new;    //could make an Amplituder::updateT or something
        Amplituder->T_elements_new.clear();
        for (int hh = 0; hh<Interaction->numOfKu; hh++){

            /*Notes:
             * I've been testing this program with CCD_newVmat.py, as it handles all diagrams
             * My python script also gives the exact same results as Audun's code, for all diagrams
             *
             * Vhhpp and denomMat works for all Nh and Nb
             *
             * La gives wrong results --> probably Vpppp that is wrong
             *
             * Vhhhh (i.e. diagrams Lb) works for Nh=2, not for Nh=14 (regardless of Nb)... weird
             *
             * makeSquareBlock works same as makeRektBlock, as I've tried using both for Vpppp -> same results
             *
             * thoughts: i think there's something wrong with making the blocks, but i can't figure out what
             *
             * Qa, which only uses Vhhpp, doesn't work either -> something iffy about my iteration?
             */

            /*Amplituder->Amplitudes[hh] = ( Interaction->Vhhpp[hh] + diagrams->La(hh)
                                                                  + diagrams->Lb(hh)
                                                                  + diagrams->Qa(hh)
                                         ).array()
                                        *Amplituder->denomMat[hh].array();*/
            //Amplituder->T_elements[];

            int ku = Interaction->Vhhpp_i[hh];
            //std::cout << "diagrams start" << std::endl;
            diagrams->Qb(ku);
            //std::cout << "diagrams end" << std::endl;

            Eigen::MatrixXd Vhhpp           = Interaction->make2x2Block(ku,0,0,1,1); //works
            Eigen::MatrixXd D_contributions = Amplituder->make2x2Block(ku,0,0,1,1, Amplituder->T_elements_new); //this contains the diagram contributions
            Eigen::MatrixXd temp = (Vhhpp + D_contributions).array()*Amplituder->denomMat[hh].array(); //works
            Amplituder->make2x2Block_inverse(temp, ku, 0,0,1,1, Amplituder->T_elements_new, false); //can't do the addition in inverse due to this line, get stacking of elements

            Eigen::MatrixXd Thhpp = Amplituder->make2x2Block(ku,0,0,1,1, Amplituder->T_elements_new);
            ECCD += 0.25*((Vhhpp.transpose())*(Thhpp)).trace();
            //ECCD += 0.25*((Interaction->Vhhpp[hh].transpose())*(Amplituder->Amplitudes[hh])).trace();
        }
        cout << std::setprecision (12) << ECCD << endl;
        conFac = abs(ECCD - ECCD_old);
        ECCD_old = ECCD;
        counter += 1;
        Amplituder->T_elements = Amplituder->T_elements_new;
        //ECCD = 0; too good to delete; you don't want to know how long i used on this
    }
    return ECCD;

}
