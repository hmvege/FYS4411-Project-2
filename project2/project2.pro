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
    samplers/uniformsampling.cpp \
    wavefunctions/n_electron_wf/nelectron.cpp \
    wavefunctions/n_electron_wf/hermite.cpp \
    wavefunctions/n_electron_wf/state.cpp \
    functions.cpp

HEADERS += \
    vmc.h \
    wavefunctions/wavefunctions.h \
    wavefunctions/twoelectronplain.h \
    wavefunctions/twoelectronjastrov.h \
    samplers/metropolissampler.h \
    samplers/importancesampler.h \
    samplers/uniformsampling.h \
    wavefunctions/n_electron_wf/nelectron.h \
    wavefunctions/n_electron_wf/hermite.h \
    wavefunctions/n_electron_wf/state.h \
    functions.h

LIBS += -llapack -lblas -larmadillo

# Following used to make armadillo usable on mac
LIBS += -L/usr/local/lib -larmadillo
INCLUDEPATH += /usr/local/include

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
