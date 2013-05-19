UTIL=../util

NC=-I/usr/local/include -I/usr/local/include/boost-numeric-bindings -I$(UTIL)
LIB=-L/usr/local/lib -L/usr/local/lib/boost-numeric-bindings -L/usr/lib/atlas-base 
CFLAGS = -lm
LIBS= -lblas -latlas -llapack -lgfortran 

OUT=.

all: hyperd

output: data 

err: data2

hyperd: hyperd.cpp $(UTIL)/print.h $(UTIL)/output.h $(UTIL)/deps.h 
	g++ $(NC) -g $(CFLAGS) hyperd.cpp -o $(OUT)/hyperd $(LIBS) $(LIB)

data: hyperd
	./hyperd approx > $(OUT)/hyperd.out
	./hyperd matrix > $(OUT)/matrix.out
	./hyperd plot0 > plot0.dat
	./hyperd plot1 > plot1.dat
	./hyperd plot3d > plot3d.dat
	gnuplot plot
	gnuplot plot3d

clean:
	rm -f hyperd hyperd.out hyperd.dat matrix.out plot3d.eps plot1.eps plot0.eps 
