UTIL=../util

NC=-I/usr/local/include -I/usr/local/include/boost-numeric-bindings -I$(UTIL)
LIB=-L/usr/local/lib -L/usr/local/lib/boost-numeric-bindings -L/usr/lib/atlas-base 
CFLAGS = -lm
LIBS= -lblas -latlas -llapack -lgfortran

OUT=.

all: splitimp

output: data 

splitimp: splitimp.cpp $(UTIL)/print.h $(UTIL)/output.h $(UTIL)/deps.h 
	g++ $(NC) -g $(CFLAGS) splitimp.cpp -o $(OUT)/splitimp $(LIBS) $(LIB)

data: splitimp
	./splitimp approx > $(OUT)/splitimp.out
	./splitimp plot0 > $(OUT)/plot0.dat
	./splitimp plot1 > $(OUT)/plot1.dat
	./splitimp plot3d > $(OUT)/plot3d.dat
	gnuplot plot
	gnuplot plot3d

clean:
	rm -f splitimp splitimp.out splitimp.dat matrix.out
