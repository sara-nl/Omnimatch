TEMPLATE = app
CONFIG -= qt
CONFIG += $$[fft]
#QMAKE_CFLAGS += -g
#QMAKE_LFLAGS_SHAPP += -g
#LIBS += -lprofiler
DESTDIR = ../bin
TARGET = omnimatch.bin
SOURCES = omnimatch.c
QMAKE_LINK = gcc
INCLUDEPATH += ../tom ../sarafft
DEPENDPATH += $$INCLUDEPATH
cufft {
    DEFINES += USE_GPUS
    LIBS += -lcufft -lcuda
} else:fftw2 {
    LIBS += -lsrfftw -lsfftw
} else:error("Don't know which fft implementation to use!")
LIBS += ../tom/libtom.a \
  ../sarafft/libsarafft.a \
  -lmpi

TARGETDEPS += ../tom/libtom.a \
  ../sarafft/libsarafft.a
