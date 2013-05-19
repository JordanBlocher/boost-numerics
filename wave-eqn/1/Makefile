UTIL=../util

NC=-I/usr/local/include -I/usr/local/include/boost-numeric-bindings -I$(UTIL)
LIB=-L/usr/local/lib -L/usr/local/lib/boost-numeric-bindings -L/usr/lib/atlas-base 
CFLAGS = -lm
LIBS= -lblas -latlas -llapack -lgfortran 

all: conv-diff

output: data 

convd-mat: convd-mat.cpp $(UTIL)/print.h $(UTIL)/deps.h 
	g++ $(NC) -g $(CFLAGS) convd-mat.cpp -o convd-mat $(LIBS) $(LIB)

conv-diff: conv-diff.cpp $(UTIL)/print.h $(UTIL)/deps.h 
	g++ $(NC) -g $(CFLAGS) conv-diff.cpp -o conv-diff $(LIBS) $(LIB)

data: conv-diff
	./conv-diff approx 0.1 a >> a.dat
	./conv-diff approx 0.5 a >> a.dat
	./conv-diff approx 1.0 a >> a.dat
	./conv-diff approx 2.0 a >> a.dat
	./conv-diff approx 0.2 b > b.dat
	./conv-diff approx 0.5 b >> b.dat
	./conv-diff approx 1.0 b >> b.dat
	./conv-diff approx 2.0 b >> b.dat
	gnuplot plot

data-mat: convd-mat
	./convd-mat approx 0.1 a >> a-mat.dat
	./convd-mat approx 0.5 a >> a-mat.dat
	./convd-mat approx 1.0 a >> a-mat.dat
	./convd-mat approx 2.0 a >> a-mat.dat
	./convd-mat approx 0.2 b > b-mat.dat
	./convd-mat approx 0.5 b >> b-mat.dat
	./convd-mat approx 1.0 b >> b-mat.dat
	./convd-mat approx 2.0 b >> b-mat.dat
	gnuplot plot-mat


clean:
	rm -f conv-diff a.dat b.dat a*.eps b*.eps
