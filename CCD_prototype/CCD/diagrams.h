#ifndef DIAGRAMS_H
#define DIAGRAMS_H

#include <eigen3/Eigen/Dense>
#include <makeintmat.h>
#include <makeampmat.h>
#include <Systems/system.h>
#include <Systems/heg.h>
#include <Systems/mp.h>
#include <Systems/pm.h>

class Diagrams
{
public:
    MakeIntMat*  m_intClass = nullptr;
    MakeAmpMat*  m_ampClass = nullptr;
    System*      m_system   = nullptr;

    Diagrams();

    //T2 diagrams
    void setIntClass(class MakeIntMat* intClass);
    void setAmpClass(class MakeAmpMat* ampClass);
    void setSystem(class System* system);
    void La(); //I made this only with the ladder diagram in mind
    void Lb();
    void Lc(); //index is a temp variable i need to check something
    void Qa();
    void Qb();
    void Qc();
    void Qd();
    void D10b();
    void D10c();

    //T3 diagrams
    void makeT3(); //performs first CCDT iteration by making T3 amplitude objects
    void T1a();
    void T1b();


    /*
     * Below are the intermediates.
     * I've called them "_terms" since they aren't actually the Ii matrices themselves,
     * but rather the terms in which they appear, eq 6.58-6.62 in Audun's thesis
    */
    void I1_term();
    void I2_term();
    void I3_term();
    void I4_term();

};

#endif // DIAGRAMS_H
