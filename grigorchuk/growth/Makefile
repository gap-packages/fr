test: diameter.C lz4.o
	g++ -Wall -O3 -fopenmp -DTEST -o $@ $^

lz4.o: lz4-read-only/lz4.c
	gcc -Wall -O3 -c -o $@ $<

diameter6: diameter.C lz4.o
	g++ -Wall -O3 -fopenmp -o $@ $^
