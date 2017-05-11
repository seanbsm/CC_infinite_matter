#include "master.h"
#include "Systems/heg.h"
#include "Systems/mp.h"
#include "makeampmat.h"
#include "makeintmat.h"
#include "diagrams.h"

#include <iostream>
#include <chrono>
#include <omp.h>
#include <iomanip> //needed for std::setprecision

typedef std::chrono::high_resolution_clock Clock;   //needed for timing

using namespace std;

void Master::setSize(int Nh, int Nb){
    m_Nh = Nh;
    m_Nb = Nb;
}

void Master::setThreads_forMaster(bool argument, int num){
    m_threadsOn = argument;
    if (argument){
        m_threads = num;
    }
}

void Master::setThreads(){
    m_intClass->setThreads(m_threads);
    m_ampClass->setThreads(m_threads);
    m_diagrams->setThreads(m_threads);
}

void Master::setSystem(class System* system){
    m_system = system;
}

void Master::setTriples(bool argument){
    m_triplesOn = argument;
}

void Master::setCCType(int type){
    m_CC_type = type;
}

void Master::setIntermediates(bool argument){
    m_intermediatesOn = argument;
}

void Master::setTimer(bool argument){
    m_timerOn = argument;
}

void Master::setRelaxation(bool argument, double alpha){
    m_relaxation = argument;
    m_alpha = alpha;
}

void Master::setClasses(){
    m_ampClass = new MakeAmpMat;
    m_intClass = new MakeIntMat;
    m_diagrams = new Diagrams;

    cout << m_Ns << endl;

    setThreads();
    m_intClass->setTriples(m_triplesOn);
    m_diagrams->setAmpClass(m_ampClass);
    m_diagrams->setIntClass(m_intClass);
    m_diagrams->setSystem(m_system);
    m_ampClass->setIntClass(m_intClass);
    m_ampClass->setSystem(m_system);

    if (1){
        std::cout << std::endl;
        m_intClass->makeBlockMat(m_system, m_Nh, m_Ns);
        std::cout << "Finished makeBlockMat" << std::endl;
        m_ampClass->makeFockMaps();
        std::cout << "Finished makeFockMaps" << std::endl;
        m_ampClass->makeDenomMat();
        std::cout << "Finished makeDenomMat" << std::endl;
        m_ampClass->setElements_T2();
        std::cout << "Finished setElements_T2" << std::endl;
        std::cout << std::endl;
    }
}

double Master::CC_Eref(){
    double Eref = 0;
    for (int i = 0; i<m_Nh;i++){
        Eref += m_system->f(i);
    }
    return Eref;
}

double Master::CC_master(double eps, double conFac){

    //would be better to simply let "m_ampClass" and "m_diagrams" inherit master?
    if (m_timerOn){

        auto t1 = Clock::now();
        setClasses();
        auto t2 = Clock::now();

        std::cout << "Time used: "
                  << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
                  << " seconds on system setup" << std::endl;
    }
    else{
        setClasses();
    }

    double ECCD_old = 0;

    for (int channel = 0; channel<m_intClass->numOfKu; channel++){
        Eigen::MatrixXd Vhhpp = m_intClass->make2x2Block(m_intClass->Vhhpp_i[channel],0,0,1,1);
        //Eigen::MatrixXd Vhhpp = m_intClass->make2x2Block_alt(h);
        Eigen::MatrixXd temp = Vhhpp.array()*m_ampClass->denomMat[channel].array();
        ECCD_old += ((Vhhpp.transpose())*(temp)).trace();

    }
    //m_ampClass->T2_elements = m_ampClass->T2_elements_new;


    std::cout << std::endl;
    std::cout << "MBPT2: " << std::setprecision (16) << ECCD_old << std::endl;
    std::cout << std::endl;

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

    if (m_triplesOn){
        m_diagrams->makeT3();
        m_ampClass->makeDenomMat3();
        std::cout << "finished makeDenomMat3" << std::endl;
        m_ampClass->emptyFockMaps();

        m_diagrams->makeT5bIndexMat();
        m_diagrams->makeT5cIndexMat();
        m_diagrams->destroy5Map();
    }
    else{
        counter ++;
    }

    std::cout << "Start of CC iterations: " << std::endl;
    while (conFac > eps /*&& counter < 8*/){
        ECCD = 0;
        //could make an m_ampClass::updateT or something
        m_ampClass->T2_elements_new.clear();

        //calculate CCD T2 diagrams
        if (counter != -1){
            if (m_intermediatesOn){
                m_diagrams->La();
                m_diagrams->I1_term1();  // Lb, Qa
                m_diagrams->I2_term1();  // Lc, Qb, due to structure of blockarrays, this is no faster than calling Lc and Qb seperatly
                m_diagrams->I3_term1();  // Qd
                m_diagrams->I4_term1();  // Qc
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
        }

        //calculate T2 contributions to T3 using T2_prev
        if(m_triplesOn){

            std::fill(m_ampClass->T3_elements_A_new.begin(), m_ampClass->T3_elements_A_new.end(), 0); //reset T3 new

            if (m_CC_type >= 1 ){
                m_diagrams->T1a();
                m_diagrams->T1b();
            }
            if (m_CC_type == 3 ){
                m_diagrams->T2c();
                m_diagrams->T2d();
                m_diagrams->T2e();
            }
            if (m_CC_type >= 2){
                // These diagrams require a re-alignment, done through a temporary map
                m_diagrams->T3b();
                m_diagrams->T3c();
                m_diagrams->T3d();
                m_diagrams->T3e();  //this is slower than the others because the remap is bigger
            }
            if (m_CC_type == 3){
                m_diagrams->T5a();
                m_diagrams->T5b();
                m_diagrams->T5c();
                m_diagrams->T5d();
                m_diagrams->T5e();
                m_diagrams->T5f();
                m_diagrams->T5g();
            }

            //update T3 amplitudes
#pragma omp parallel for num_threads(m_threads)
            for (int i=0; i<m_ampClass->T3_elements_A_new.size(); i++){
                m_ampClass->T3_elements_A_new[i] *= m_ampClass->denom3_elements[i];
            }

            if (m_relaxation){
                std::vector<double> T3_temp = m_ampClass->T3_elements_A;

                for(int it=0; it<m_ampClass->T3_elements_A_new.size(); it++){
                    m_ampClass->T3_elements_A[it] = m_alpha*m_ampClass->T3_elements_A_new[it] + (1-m_alpha)*T3_temp[it];
                }
            }
            else{
                m_ampClass->T3_elements_A = m_ampClass->T3_elements_A_new;
            }


            /*int zeros = 0;
            for (auto& T: m_ampClass->T3_elements_A)
                if (T == 0.){zeros ++;}
            std::cout << zeros << std::endl;*/

            //calculate T3 contributions to T2 using T3_current
            m_diagrams->D10b();
            m_diagrams->D10c();

        }

        //update T2 amplitudes
        for (int hh = 0; hh<m_intClass->numOfKu; hh++){
            int ku = m_intClass->Vhhpp_i[hh];

            Eigen::MatrixXd Vhhpp           = m_intClass->make2x2Block(ku,0,0,1,1);
            Eigen::MatrixXd D_contributions = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T2_elements_new);
            Eigen::MatrixXd temp = (Vhhpp + D_contributions).array()*m_ampClass->denomMat[hh].array();

            m_ampClass->make2x2Block_inverse(temp, ku, 0,0,1,1, m_ampClass->T2_elements_new, false);

            Eigen::MatrixXd Thhpp = m_ampClass->make2x2Block(ku,0,0,1,1, m_ampClass->T2_elements_new);
            ECCD += ((Vhhpp.transpose())*(Thhpp)).trace();
        }


        cout << std::fixed << std::setprecision (16) << ECCD << endl;

        conFac = abs(ECCD - ECCD_old);
        ECCD_old = ECCD;
        counter += 1;

        if (m_relaxation){
            spp::sparse_hash_map<unsigned long int, double> T2_temp = m_ampClass->T2_elements;
            m_ampClass->T2_elements.clear();
            for(auto const& it : m_ampClass->T2_elements_new) {
                m_ampClass->T2_elements[it.first] = m_alpha*it.second + (1-m_alpha)*T2_temp[it.first];
            }
        }
        else{
            m_ampClass->T2_elements = m_ampClass->T2_elements_new;
        }

        //ECCD = 0; too good to delete; you don't want to know how long i used to find this
    }
    return ECCD;
}
