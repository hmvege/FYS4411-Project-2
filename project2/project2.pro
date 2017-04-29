TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    metropolis.cpp \
    wavefunctions.cpp \
    vmc.cpp

HEADERS += \
    metropolis.h \
    wavefunctions.h \
    vmc.h
