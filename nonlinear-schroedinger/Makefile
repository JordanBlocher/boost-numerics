UTIL=util

NC=-I/usr/local/include -I/usr/local/include/boost-numeric-bindings -I$(UTIL)
LIB=-L/usr/local/lib -L/usr/lib/atlas-base 
CFLAGS = -lm
LIBS= -lblas -latlas -llapack -lgfortran

all: nls-explicit nls-crnk 

latex: proj2.pdf

nls-explicit: nls-explicit.cpp $(UTIL)/print.h $(UTIL)/deps.h 
	g++ $(NC) -g $(CFLAGS) nls-explicit.cpp -o nls-explicit $(LIBS) $(LIB)

nls-crnk: nls-crnk.cpp $(UTIL)/print.h $(UTIL)/deps.h 
	g++ $(NC) -g $(CFLAGS) nls-crnk.cpp -o nls-crnk $(LIBS) $(LIB)

data: 
	./nls-explicit approx > nls.out
	./nls-explicit energy > nls-explicit.out
	./nls-crnk approx >> nls.out
	./nls-crnk energy >> nls-crnk.out
	./nls-explicit plot > nls.dat
	./nls-crnk plot >> nls.dat
	gnuplot plot

clean:
	rm -f nls-explicit nls-crnk *.eps *.dat *.out proj2.pdf latex/proj2.dvi latex/proj2.aux latex/proj2.log

proj2.dvi:
	latex -output-directory=latex latex/proj2.tex

proj2.pdf: proj2.dvi
	dvipdf latex/proj2.dvi
