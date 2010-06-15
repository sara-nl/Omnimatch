TEMPLATE = lib
CONFIG += staticlib $$[fft]
CONFIG -= qt
TARGET = sarafft
cufft {
    INCLUDEPATH += ../omnicuda
    DEFINES += USE_GPUS
}
SOURCES += sarafft.c
HEADERS += sarafft.h
DEPENDPATH += $$INCLUDEPATH
