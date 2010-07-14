TEMPLATE = lib
#QMAKE_CFLAGS += -g
#QMAKE_CC = cc
TARGET = tom
CONFIG -= qt
CONFIG += static
#HEADERS += nrutil.h \
#	tom.h \
#	zeit.h
INCLUDEPATH += ../include
SOURCES += cross.c \
	energizer.c \
	fourfilter.c \
	pastes.c \
	real_utils.c \
	shift.c \
	sort4fftw.c \
	tom_rotate3d.c
