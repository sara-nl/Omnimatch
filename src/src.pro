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
LIBS += ../sarafft/libsarafft.a ../tom/libtom.a -lmpi
cufft {
  LIBS += ../omnicuda/libomnicuda.a -lcufft -lcuda
} else:fftw2 {
  LIBS += -lsrfftw -lsfftw
} else:error("Don't know which fft implementation to use!")
