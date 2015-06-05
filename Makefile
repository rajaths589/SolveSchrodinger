all: solver.o operator.o state.o wavefninit.o
	gcc -g ./bin/solver.o ./bin/operator.o ./bin/state.o ./bin/wavefninit.o -o ./bin/soursugar.o -lgsl -lgslcblas -lm -lfftw3	

operator.o: operator.c
	gcc -g -c operator.c -o ./bin/operator.o

state.o: state.c
	gcc -g -c state.c -o ./bin/state.o

wavefninit.o: wavefninit.c
	gcc -g -c wavefninit.c -o ./bin/wavefninit.o

solver.o: solver.c
	gcc -g -c solver.c -o ./bin/solver.o

clean:
	cd ./bin
	rm -f *.o
	cd ..
