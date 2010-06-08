TEMPLATE = app
CONFIG -= qt
#QMAKE_CFLAGS += -g
#QMAKE_LFLAGS_SHAPP += -g
#LIBS += -lprofiler
DESTDIR = ../bin
TARGET = omnimatch.bin
SOURCES = omnimatch.c
QMAKE_LINK = gcc
INCLUDEPATH += ../tom ../sarafft
LIBS += ../sarafft/libsarafft.a ../tom/libtom.a -lmpi -lsrfftw -lsfftw ../omnicuda/libomnicuda.a -lcufft -lcuda

TARGETDEPS += ../tom/libtom.a ../sarafft/libsarafft.a ../omnicuda/libomnicuda.a
