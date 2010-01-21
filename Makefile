#
#
#

all:
	rm -f ./src/*.o
	rm -f ./lib/*.o
	rm -f ./lib/libtom.a
	cd ./lib/ && $(MAKE)
	cd ./src/ && $(MAKE)

clean:
	rm -f ./src/*.o  
	rm -f ./lib/*.o  
	rm -f ./lib/libtom.a     

libs:
	cd ./lib/ && $(MAKE)
