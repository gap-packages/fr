all: layout hsolve

layout: layout.C
	g++ -g $< -o $@ -L. -llevmar -llapack -lblas

layout-HLBFGS: layout-HLBFGS.C
	g++ -g $< -o $@ -L. -lHLBFGS

hsolve: hsolve.C
	g++ -g $< -o $@ -L. -llevmar -llapack -lblas
