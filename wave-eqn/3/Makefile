UTIL=../util

NC=-I/usr/local/include -I/usr/local/include/boost-numeric-bindings -I$(UTIL)
LIB=-L/usr/local/lib -L/usr/local/lib/boost-numeric-bindings -L/usr/lib/atlas-base 
CFLAGS = -lm
LIBS= -lblas -latlas -llapack -lgfortran 

all: crnk

output: data 

crnk: crnk.cpp $(UTIL)/print.h $(UTIL)/deps.h 
	g++ $(NC) -g $(CFLAGS) crnk.cpp -o crnk $(LIBS) $(LIB)

data: crnk
	./crnk approx 0.0 a > a.dat
	./crnk approx 0.12 a >> a.dat
	./crnk approx 0.2 a >> a.dat
	./crnk approx 0.8 a >> a.dat
	./crnk approx 0.0 b > b.dat
	./crnk approx 5.0 b >> b.dat
	./crnk approx 10.0 b >> b.dat
	./crnk approx 20.0 b >> b.dat
	gnuplot plot

clean:
	rm -f crnk a.dat b.dat a*.eps 
