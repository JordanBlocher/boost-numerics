UTIL=../util

NC=-I/usr/local/include -I/usr/local/include/boost-numeric-bindings -I$(UTIL)
LIB=-L/usr/local/lib -L/usr/lib/atlas-base 
CFLAGS = -lm
LIBS= -lblas -latlas -llapack -lgfortran

OUT=.

all: hyperd

output: data 

hyperd: hyperd.cpp $(UTIL)/print.h $(UTIL)/output.h $(UTIL)/deps.h 
	g++ $(LIB) $(NC) -g $(CFLAGS) $(LIBS) hyperd.cpp -o $(OUT)/hyperd

data: hyperd
	./hyperd approx > $(OUT)/hyperd.out
	./hyperd matrix > $(OUT)/matrix.out
	./hyperd plot0 > plot0.dat
	./hyperd plot1 > plot1.dat
	./hyperd plot3d > plot3d.dat
	gnuplot plot 
	gnuplot plot3d

clean:
	rm -f hyperd hyperd.out
