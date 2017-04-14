TEMPLATE = app
CONFIG  += console
CONFIG  -= app_bundle
CONFIG  -= qt
CONFIG(release, debug|release){CONFIG += optimize_full}
QMAKE_CFLAGS = -m64
QMAKE_LFLAGS = -m64 -fopenmp
QMAKE_LFLAGS -= -O1
QMAKE_LFLAGS += -O3
QMAKE_CXXFLAGS = -m64 -fopenmp -std=c++11

SOURCES += main.cpp \
    Systems/heg.cpp \
    Systems/mp.cpp \
    Systems/pm.cpp \
    master.cpp \
    Systems/system.cpp \
    makeampmat.cpp \
    makeintmat.cpp \
    diagrams.cpp

HEADERS += \
    Systems/heg.h \
    Systems/mp.h \
    Systems/pm.h \
    master.h \
    Systems/system.h \
    makeampmat.h \
    makeintmat.h \
    diagrams.h

INCLUDEPATH += /usr/include/openmpi-x86_64
QMAKE_CXX = mpicxx #/usr/lib64/openmpi/bin/mpicxx #mpicxx
QMAKE_CXX += -O3
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX

QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = /usr/bin/mpicxx #/usr/lib64/openmpi/bin/mpicc#mpicc

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
