all: fit.pdf

fit.pdf: plot.py salida.dat
	python plot.py

salida.dat: mcmc_solar.x
	./mcmc_solar.x 4622 10 > salida.dat

mcmc_solar.x : mcmc_solar.c
	cc mcmc_solar.c  -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas -lm -o mcmc_solar.x


