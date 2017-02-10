#include "master.h"
#include "Systems/heg.h"
#include "Systems/mp.h"
#include "makeampmat.h"
#include "makeintmat.h"
#include "diagrams.h"

#include <iostream>
#include <chrono>
#include <iomanip> //needed for std::setprecision

typedef std::chrono::high_resolution_clock Clock;   //needed for timing

using namespace std;

void Master::setSize(int Nh, int Nb){
    m_Nh = Nh;
    m_Nb = Nb;
}

void Master::setSystem(class System* system){
    m_system = system;
}

void Master::setTriples(bool argument){
    m_triplesOn = argument;
}

void Master::setIntermediates(bool argument){
    m_intermediatesOn = argument;
}

void Master::setTimer(bool argument){
    m_timerOn = argument;
}

void Master::setClasses(){
    m_ampClass = new MakeAmpMat;
    m_intClass = new MakeIntMat;
    m_diagrams = new Diagrams;
    cout << m_Ns << endl;
    m_intClass->makeBlockMat(m_system, m_Nh, m_Ns);

    m_diagrams->setAmpClass(m_ampClass);
    m_diagrams->setIntClass(m_intClass);
    m_diagrams->setSystem(m_system);
    m_ampClass->setIntClass(m_intClass);
    m_ampClass->setSystem(m_system);
    m_ampClass->setElements_T2();
    m_ampClass->makeDenomMat();
}

double Master::CC_master(double eps, double conFac){

    //would be better to simply let "m_ampClass" and "m_diagrams" inherit master?
    if (m_timerOn){

        auto t1 = Clock::now();
        setClasses();
        auto t2 = Clock::now();

        std::cout << "Time used: "
                  << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
                  << " seconds on making indexHolders" << std::endl;
    }
    else{
        setClasses();
    }

    double ECCD_old = 0;
    for (int h = 0; h<m_intClass->numOfKu; h++){
        Eigen::MatrixXd Vhhpp = m_intClass->make2x2Block(m_intClass->Vhhpp_i[h],0,0,1,1);
        //Eigen::MatrixXd Vhhpp = m_intClass->make2x2Block_alt(h);
        Eigen::MatrixXd temp = Vhhpp.array()*m_ampClass->denomMat[h].array();
        ECCD_old += ((Vhhpp.transpose())*(temp)).trace();

    }
    //check whether or not to multiply by 0.25 for MBPT2
    std::cout << "MBPT2: " << std::setprecision (12) << 0.25*ECCD_old << std::endl;



    double ECC;
    if (m_timerOn){

        auto t1 = Clock::now();
        ECC = Iterator(eps, conFac, ECCD_old);
        auto t2 = Clock::now();

        std::cout << "Time used: "
                  << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
                  << " seconds on solving CC" << std::endl;
    }
    else{
        ECC = Iterator(eps, conFac, ECCD_old);
    }

    return ECC;
}

double Master::Iterator(double eps, double conFac, double E_MBPT2){
    int counter = 0;
    double ECCD_old = E_MBPT2;
    double ECCD     = 0;

    while (conFac > eps && counter < 5e1){
        ECCD = 0;
        //could make an m_ampClass::updateT or something
        m_ampClass->T2_elements_new.clear();

        if (m_intermediatesOn){
            m_diagrams->La();
            m_diagrams->I1_term();  // Lb, Qa
            m_diagrams->I2_term();  // Lc, Qb, due to structure of blockarrays, this is no faster than calling Lc and Qb seperatly
            m_diagrams->I3_term();  // Qd
            m_diagrams->I4_term();  // Qc
        }
        else{
            //CCD diagrams
            m_diagrams->La();
            m_diagrams->Lb();
            m_diagrams->Lc();
            m_diagrams->Qa();
            m_diagrams->Qb();
            m_diagrams->Qc();
            m_diagrams->Qd();
        }

        if(m_triplesOn){
            m_diagrams->T2a();
            m_diagrams->D10b();
            m_diagrams->D10c();
        }

        for (int hh = 0; hh<m_intClass->numOfKu; hh++){
            int ku = m_intClass->Vhhpp_i[hh];

            Eigen::MatrixXd Vhhpp           = m_intClass->make2x2Block(ku,0,0,1,1);
            Eigen::MatrixXd D_contributions = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T2_elements_new);
            Eigen::MatrixXd temp = (Vhhpp + D_contributions).array()*m_ampClass->denomMat[hh].array();

            m_ampClass->make2x2Block_inverse(temp, ku, 0,0,1,1, m_ampClass->T2_elements_new, false);

            Eigen::MatrixXd Thhpp = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T2_elements_new);
            ECCD += 0.25*((Vhhpp.transpose())*(Thhpp)).trace();
        }

        cout << std::setprecision (12) << ECCD << endl;

        conFac = abs(ECCD - ECCD_old);
        ECCD_old = ECCD;
        counter += 1;

        if (0){
            double alpha = 0.2;
            std::map<int, double> T2_temp = m_ampClass->T2_elements;
            m_ampClass->T2_elements.clear();
            for(auto const& it : m_ampClass->T2_elements_new) {
                m_ampClass->T2_elements[it.first] = alpha*it.second + (1-alpha)*T2_temp[it.first];
            }
        }
        else{
            m_ampClass->T2_elements = m_ampClass->T2_elements_new;
        }

        //ECCD = 0; too good to delete; you don't want to know how long i used to find this
    }
    return ECCD;
}
