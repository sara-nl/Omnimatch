TEMPLATE = app
CONFIG -= qt
#QMAKE_CFLAGS += -g
#QMAKE_CC = cc
QMAKE_LINK = gcc
#QMAKE_LFLAGS_SHAPP += -g
INCLUDEPATH += ../tom
DEPENDPATH += ../tom
LIBS += -L../tom -ltom -lmpi
LIBS += -lsrfftw -lsfftw
#LIBS += -lrfftw -lfftw
#LIBS += -lprofiler
DESTDIR = ../bin
TARGET = omnimatch.bin
SOURCES = omnimatch.c
