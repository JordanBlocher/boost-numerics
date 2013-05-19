NC=-I/usr/local/include -I/usr/local/include/boost-numeric-bindings
LIB=-L/usr/local/lib -L/usr/lib/atlas-base 
CFLAGS = -lm
LIBS= -lblas -latlas -llapack -lgfortran

ifeq ($(USER), 'frags') 
PLOT=gplotm
else 
PLOT=gplot
endif

all: heat

output: data1 data2 data3 plots

heat: heat.cpp ../deps.h ../print.h ../outputsmall.h
	g++ $(LIB) $(NC) -g -o heat heat.cpp $(CFLAGS) $(LIBS)

data1: heat
	./heat approx 20 1.0 1000 1 > out/data1.out
	./heat plot 20 1.0 1000 1 > data/plot1.dat
	./heat plot3d 20 1.0 1000 1 > data/plot3d1.dat

data2: heat
	./heat approx 20 0.01 1000 1 > out/data2.out
	./heat plot 20 0.01 1000 1 > data/plot2.dat
	./heat plot3d 20 0.01 1000 1 > data/plot3d2.dat
	
data3: heat
	./heat plot 40 1.0 1000 0 > data/plotd1.dat
	./heat plotd 40 1.0 1000 0 > data/plotdd1.dat
	./heat plot3d 40 1.0 1000 0 > data/plotd3d1.dat
	./heat approx 40 1.0 1000 1 > out/data3a.out
	./heat plot 40 1.0 1000 1 > data/plot3a.dat
	./heat plot3d 40 1.0 1000 1 > data/plot3d3a.dat
	./heat approx 40 0.01 1000 1 > out/data3b.out
	./heat plot 40 0.01 1000 1 > data/plot3b.dat
	./heat plot3d 40 0.01 1000 1 > data/plot3d3b.dat


plots:
	gnuplot $(PLOT)/plot
	gnuplot $(PLOT)/plot3d
	
clean:
	rm -f heat out/* data/* eps/* 
