TEMPLATE = lib
QMAKE_CFLAGS_STATIC_LIB -= -fPIC
QMAKE_CFLAGS -= -pipe
TARGET = omnicuda
CONFIG += static \
 staticlib
HEADERS += omnicuda.h
SOURCES += omnicuda.cu
CONFIG -= qt warn_on
QMAKE_CC = nvcc
QMAKE_LINK = gcc
