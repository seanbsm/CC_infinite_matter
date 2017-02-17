#include "testclass.h"


testClass::testClass()
{
}

void testClass::testFunc(){

    for (int i=0; i<1e4; i++){
        places[i] = (double) i;
    }
}
