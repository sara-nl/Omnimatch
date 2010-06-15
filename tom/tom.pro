TEMPLATE = lib
#QMAKE_CFLAGS += -g
#QMAKE_CC = cc
TARGET = tom
CONFIG -= qt
CONFIG += static \
 staticlib
HEADERS += tom.h
SOURCES += cross.c \
	emfile.c \
	energizer.c \
	fourfilter.c \
	pastes.c \
	real_utils.c \
	shift.c \
	sort4fftw.c \
	tom_rotate3d.c
TARGETDEPS += ../sarafft/libsarafft.a

INCLUDEPATH += ../sarafft
DEPENDPATH += $$INCLUDEPATH
