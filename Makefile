all: solver.o operator.o state.o wavefninit.o
	gcc -g solver.o operator.o state.o wavefninit.o -o soursugar.o -lgsl -lgslcblas -lm -lfftw3

operator.o: operator.c
	gcc -g -c operator.c

state.o: state.c
	gcc -g -c state.c

wavefninit.o: wavefninit.c
	gcc -g -c wavefninit.c

solver.o: solver.c
	gcc -g -c solver.c

clean:
	rm *.o
