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

output: data plots

heat: heat.cpp ../deps.h ../print.h ../output.h
	g++ $(LIB) $(NC) -g -o heat heat.cpp $(CFLAGS) $(LIBS)

data: heat
	./heat all > out/data.out
	./heat plot > data/plot.dat
	./heat ploterf > data/ploterf.dat
	./heat plotsol > data/plotsol.dat
	./heat plot3d > data/plot3d.dat

plots:
	gnuplot $(PLOT)/plot
	gnuplot $(PLOT)/ploterf
	gnuplot $(PLOT)/plotsol
	gnuplot $(PLOT)/plot3d

clean:
	rm -f heat out/* data/* eps/* 
