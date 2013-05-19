NC=-I/usr/local/include -I/usr/local/include/boost-numeric-bindings
LIB=-L/usr/local/lib -L/usr/lib/atlas-base 
CFLAGS = -lm
LIBS= -lblas -latlas -llapack -lgfortran

all: output latex

output: 2-1 3-1 5-123 

latex: hw1.pdf

2-1: 2-1/heat deps.h
	make -C 2-1/ output

3-1: 3-1/heat deps.h
	make -C 3-1/ output

5-123: 5-123/heat deps.h
	make -C 5-123/ output

5-7:
	
all-clean:
	make clean -C 2-1/
	make clean -C 3-1/
	make clean -C 5-123/
	rm -f hw1.pdf latex/hw1.dvi latex/hw1.aux latex/hw1.log

clean:
	rm -f hw1.pdf latex/hw1.dvi latex/hw1.aux latex/hw1.log

hw1.dvi:
	latex -output-directory=latex latex/hw1.tex

hw1.pdf: hw1.dvi
	dvipdf latex/hw1.dvi
