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
	./heat all 10 20 > out/dataa.out
	./heat plot 10 20 > data/plota.dat
	./heat plot3d 10 20 > data/plota3d.dat
	./heat all 20 100 > out/datab.out
	./heat plot 20 100 > data/plotb.dat
	./heat plot3d 20 100 > data/plotb3d.dat

plots:
	gnuplot $(PLOT)/plot
	gnuplot $(PLOT)/plot3d

clean:
	rm -f heat out/* data/* eps/* 
