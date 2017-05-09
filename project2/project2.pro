TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    vmc.cpp \
    wavefunctions/wavefunctions.cpp \
    wavefunctions/twoelectronplain.cpp \
    wavefunctions/twoelectronjastrov.cpp \
    samplers/metropolissampler.cpp \
    samplers/importancesampler.cpp \
    samplers/uniformsampling.cpp

HEADERS += \
    vmc.h \
    wavefunctions/wavefunctions.h \
    wavefunctions/twoelectronplain.h \
    wavefunctions/twoelectronjastrov.h \
    samplers/metropolissampler.h \
    samplers/importancesampler.h \
    samplers/uniformsampling.h
