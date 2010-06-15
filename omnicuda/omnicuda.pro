TEMPLATE = lib
TARGET = omnicuda
CONFIG += static \
 staticlib
HEADERS += omnicuda.h
SOURCES += omnicuda.cpp
CONFIG -= qt
INCLUDEPATH += ../sarafft
DEPENDPATH = $$INCLUDEPATH
DEFINES += USE_GPUS
