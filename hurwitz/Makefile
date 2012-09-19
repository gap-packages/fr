all: layout hsolve

layout: layout.C
	g++ -g -I/usr/include/suitesparse -I/usr/local/src/libdogleg $< -o $@ -Wl,-rpath,/usr/local/src/libdogleg -L/usr/local/src/libdogleg -ldogleg -lcholmod

hsolve: hsolve.C
	g++ -g -I/usr/local/src/levmar-2.6 $< -o $@ -L/usr/local/src/levmar-2.6 -llevmar -llapack -lblas
