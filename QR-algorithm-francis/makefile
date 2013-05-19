NC=-I/usr/local/include -I/usr/local/include/boost-numeric-bindings
LIB=-L/usr/local/lib -L/usr/lib/atlas-base 
CFLAGS = -lm
LIBS= -lblas -latlas -llapack -llapacke -lgfortran

all: latex hessenberg qr_factorization qr_shifted

latex: qr-francis.pdf

hessenberg: hessenberg.cpp ./util/deps.h ./util/print.h
	g++ $(LIB) $(NC) -g -o hessenberg hessenberg.cpp $(CFLAGS) $(LIBS)

qr_factorization: qr_factorization.cpp ./util/deps.h ./util/print.h
	g++ $(LIB) $(NC) -g -o qr_factorization qr_factorization.cpp $(CFLAGS) $(LIBS)

qr_shifted: qr_shifted.cpp ./util/hessenberg.h ./util/deps.h ./util/print.h
	g++ $(LIB) $(NC) -g -o qr_shifted qr_shifted.cpp $(CFLAGS) $(LIBS)

clean:
	rm -f hessenberg qr_factorization qr_shifted 

qr-francis.pdf:
	pdflatex -output-directory=./latex -aux-directory=./latex ./latex/qr-francis.tex

