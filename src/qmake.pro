TEMPLATE = app
CONFIG -= qt
#QMAKE_CFLAGS += -g
#QMAKE_LFLAGS_SHAPP += -g
#LIBS += -lprofiler
DESTDIR = ../bin
TARGET = omnimatch.bin
SOURCES = omnimatch.c
INCLUDEPATH += ../omnicuda \
  ../tom

TARGETDEPS += ../omnicuda/libomnicuda.a \
  ../tom/libtom.a

LIBS += ../omnicuda/libomnicuda.a \
  ../tom/libtom.a \
  -lcufft \
  -lcuda \
  -lmpi \
  -lrfftw \
  -lfftw

