TEMPLATE = lib
CONFIG += staticlib $$[fft]
CONFIG -= qt
TARGET = sarafft
cufft {
  DEFINES += USE_GPUS
  SOURCES += saracufft.cpp
}
fftw2 {
  SOURCES += sarafftw2.c
}
HEADERS += sarafft.h
