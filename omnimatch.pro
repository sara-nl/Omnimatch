TEMPLATE = subdirs
CONFIG -= gui qt core
CONFIG += ordered $$[fft]
#cufft:SUBDIRS = omnicuda
SUBDIRS += sarafft tom src
QMAKE_CLEAN += bin/*.ccf bin/*.ang bin/*.err bin/*.out bin/omnimatch*.job.[oe]*
