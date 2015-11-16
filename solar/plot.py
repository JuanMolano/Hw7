##
##  plot.py
##
##
##  Created by Juan Pablo Molano on 15/11/15.
##
##

import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("salida.dat")

plt.figure()
f,fig = plt.subplots(3,3,figsize=(10,10))

fig[0,0].scatter(data[:,0] , data[:,1], 0.01)
fig[0,0].plot(data[:,0] , data[:,7], color = 'red')
fig[0,0].tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
fig[0,0].set_title("$Fit$")
fig[0,1].plot( data[:,0] , data[:,2] , color = 'black')
fig[0,1].tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
fig[0,1].set_title("$Chi$")
fig[1,0].scatter( data[:,3] , data[:,4], 0.01)
fig[1,0].set_ylabel('a')
fig[1,0].set_xlabel('b')
fig[1,1].scatter(data[:,3] , data[:,5], 0.01)
fig[1,1].set_ylabel('a')
fig[1,1].set_xlabel('c')
fig[1,2].scatter( data[:,3] , data[:,6], 0.01)
fig[1,2].set_ylabel('a')
fig[1,2].set_xlabel('d')
fig[2,0].scatter(data[:,4] , data[:,5], 0.01)
fig[2,0].set_ylabel('b')
fig[2,0].set_xlabel('c')
fig[2,1].scatter(data[:,4] , data[:,6], 0.01)
fig[2,1].set_ylabel('b')
fig[2,1].set_xlabel('d')
fig[2,2].scatter(data[:,5] , data[:,6], 0.01)
fig[2,2].set_ylabel('c')
fig[2,2].set_xlabel('a')
fig[2,2].tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
fig[0,2].axis('off')



plt.savefig("fit.pdf",format="pdf")
plt.close()