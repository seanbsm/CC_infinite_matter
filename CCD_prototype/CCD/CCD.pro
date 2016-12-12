TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -std=c++11

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

