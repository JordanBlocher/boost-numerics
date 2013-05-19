UTIL=../util

NC=-I/usr/local/include -I/usr/local/include/boost-numeric-bindings -I$(UTIL)
LIB=-L/usr/local/lib -L/usr/lib/atlas-base 
CFLAGS = -lm
LIBS= -lblas -latlas -llapack -lgfortran

OUT=.

all: sivash

output: data 

sivash: sivash.cpp $(UTIL)/print.h $(UTIL)/output.h $(UTIL)/deps.h 
	g++ -g sivash.cpp -o $(OUT)/sivash $(LIB) $(NC) $(CFLAGS) $(LIBS) 

data: sivash
	./sivash approx > $(OUT)/sivash.out
	./sivash matrix > matrix.out
	./sivash plot0 > plot0.dat
	./sivash plot1 > plot1.dat
	./sivash plot3d > plot3d.dat
	gnuplot plot
	gnuplot plot3d

clean:
	rm -f sivash sivash.out plot0.dat plot1.dat plot3d.dat matrix.out plot1.eps plot2.eps plot3d.eps
