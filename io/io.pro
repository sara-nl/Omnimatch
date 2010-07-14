TEMPLATE = lib
#QMAKE_CFLAGS += -g
#QMAKE_CC = cc
TARGET = io
CONFIG -= qt
CONFIG += static
HEADERS += datatypes.h\
           em.h \
           mrc.h \
           data.h
INCLUDEPATH += ../include
SOURCES += em.c \
	   mrc.c \
           data.c
