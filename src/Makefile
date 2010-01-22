#
#
#

CFLAGS = -g  -I. -I../include/ -I../../../include

LD = $(CC)
LDFLAGS = -g  -L. -L../../../lib -L../lib/
LIBS = -lrfftw -lfftw -lm -ltom -lmpi

OBJECTS = ./omnimatch.o  

PROGRAM = ../bin/omnimatch.exe

$(PROGRAM): $(OBJECTS)
	$(LD) $(LDFLAGS) -o $(PROGRAM) $(OBJECTS) $(LIBS)


clean:
	rm -f ./src/*.o  
	rm -f ./lib/*.o  
	rm -f ./lib/libtom.a     

libs:
	cd ./lib/ && $(MAKE)
