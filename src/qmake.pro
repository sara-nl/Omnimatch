TEMPLATE = app
CONFIG -= qt
#QMAKE_CFLAGS += -g
#QMAKE_LFLAGS_SHAPP += -g
#LIBS += -lprofiler
DESTDIR = ../bin
TARGET = omnimatch.bin
SOURCES = omnimatch.c



QMAKE_LINK = gcc



INCLUDEPATH += ../sarafft \
  ../tom

LIBS += ../omnicuda/libomnicuda.a \
  ../tom/libtom.a \
  ../sarafft/libsarafft.a \
  -lcufft \
  -lcuda \
  -lmpi \
  -lrfftw \
  -lfftw

TARGETDEPS += ../sarafft/libsarafft.a \
  ../omnicuda/libomnicuda.a \
  ../tom/libtom.a

