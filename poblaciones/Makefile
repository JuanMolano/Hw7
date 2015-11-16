all: final.pdf

final.pdf: salida.dat plot.py
	python plot.py

salida.dat: mcmc_loktavolterra.x new_data.txt 
	./mcmc_loktavolterra.x 95 1 > salida.dat

new_data.txt: lotka_volterra_obs.dat
	tail -n +2 lotka_volterra_obs.dat > new_data.txt

mcmc_loktavolterra.x : mcmc_lotkavolterra.c
	cc mcmc_lotkavolterra.c  -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas -lm -o mcmc_loktavolterra.x



