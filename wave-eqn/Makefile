NC=-I/usr/local/include -I/usr/local/include/boost-numeric-bindings
LIB=-L/usr/local/lib -L/usr/lib/atlas-base 
CFLAGS = -lm
LIBS= -lblas -latlas -llapack -lgfortran

all: output latex

latex: wave.pdf

clean:
	rm -f wave.pdf latex/wave.dvi latex/wave.aux latex/wave.log

wave.dvi:
	latex -output-directory=latex latex/wave.tex

wave.pdf: wave.dvi
	dvipdf latex/wave.dvi
