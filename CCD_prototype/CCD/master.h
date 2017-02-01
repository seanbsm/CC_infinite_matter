#ifndef MASTER_H
#define MASTER_H

#include <Systems/system.h>

class Master
{
public:
    class System*  m_system             = nullptr;
    int            m_Nh                 = 0;        //number of particles
    int            m_Nb                 = 0;
    int            m_Ns                 = 0;
    bool           m_triplesOn;
    bool           m_intermediatesOn;
    void           setSize();

    void setSystem(class System* system);
    void setSize(int Nh, int Nb);
    void setIntermediates(bool argument);
    void setTriples(bool argument);

    double Iterator(double eps, double conFac);
};

#endif // MASTER_H
