TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    vmc.cpp \
    wavefunctions/wavefunctions.cpp \
    wavefunctions/twoelectronplain.cpp \
    wavefunctions/twoelectronjastrov.cpp \
    ratios/metropolisratio.cpp \
    ratios/importancesampler.cpp

HEADERS += \
    vmc.h \
    wavefunctions/wavefunctions.h \
    wavefunctions/twoelectronplain.h \
    wavefunctions/twoelectronjastrov.h \
    ratios/metropolisratio.h \
    ratios/importancesampler.h
