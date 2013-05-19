UTIL=../util

NC=-I/usr/local/include -I/usr/local/include/boost-numeric-bindings -I$(UTIL)
LIB=-L/usr/local/lib -L/usr/local/lib/boost-numeric-bindings -L/usr/lib/atlas-base 
CFLAGS = -lm
LIBS= -lblas -latlas -llapack -lgfortran 

all: ftfs

output: data 

ftfs-mat: ftfs-mat.cpp $(UTIL)/print.h $(UTIL)/deps.h 
	g++ $(NC) -g $(CFLAGS) ftfs-mat.cpp -o ftfs-mat $(LIBS) $(LIB)

ftfs: ftfs.cpp $(UTIL)/print.h $(UTIL)/deps.h 
	g++ $(NC) -g $(CFLAGS) ftfs.cpp -o ftfs $(LIBS) $(LIB)

data: ftfs
	./ftfs approx 0.0 a > a.dat
	./ftfs approx 0.12 a >> a.dat
	./ftfs approx 0.2 a >> a.dat
	./ftfs approx 0.8 a >> a.dat
	./ftfs approx 0.0 b > b.dat
	./ftfs approx 5.0 b >> b.dat
	./ftfs approx 10.0 b >> b.dat
	./ftfs approx 20.0 b >> b.dat
	gnuplot plot

data-mat: ftfs-mat
	./ftfs-mat approx 0.0 a > a.dat
	./ftfs-mat approx 0.12 a >> a.dat
	./ftfs-mat approx 0.2 a >> a.dat
	./ftfs-mat approx 0.8 a >> a.dat
	./ftfs-mat approx 0.0 b > b.dat
	./ftfs-mat approx 5.0 b >> b.dat
	./ftfs-mat approx 10.0 b >> b.dat
	./ftfs-mat approx 20.0 b >> b.dat
	gnuplot plot


clean:
	rm -f ftfs a.dat b.dat a*.eps b*.eps
