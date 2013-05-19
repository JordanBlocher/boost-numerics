UTIL=../util

NC=-I/usr/local/include -I/usr/local/include/boost-numeric-bindings -I$(UTIL)
LIB=-L/usr/local/lib -L/usr/lib/atlas-base 
CFLAGS = -lm
LIBS= -lblas -latlas -llapack -lgfortran

OUT=.

all: shivas

output: data 

shivas: shivas.cpp $(UTIL)/print.h $(UTIL)/output.h $(UTIL)/deps.h 
	g++ $(LIB) $(NC) -g $(CFLAGS) $(LIBS) shivas.cpp -o $(OUT)/shivas

data: shivas
	./shivas approx > $(OUT)/shivas.out

clean:
	rm -f shivas shivas.out shivas.dat
