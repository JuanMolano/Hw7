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
m = np.loadtxt("new_data.txt")

plt.figure()
f,fig = plt.subplots(3,3,figsize=(10,10))

fig[0,0].scatter(m[:,0] , m[:,1], 1)
fig[0,0].tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
fig[0,0].set_title("$Fit$")
fig[0,1].scatter( m[:,0] , m[:,2] , 1)
fig[0,1].tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
fig[0,1].set_title("$Chi$")
fig[1,0].scatter( data[:,6] , data[:,7], 1)
fig[1,0].set_ylabel('a')
fig[1,0].set_xlabel('b')
fig[1,1].scatter(data[:,6] , data[:,8], 1)
fig[1,1].set_ylabel('a')
fig[1,1].set_xlabel('c')
fig[1,2].scatter( data[:,6] , data[:,9], 1)
fig[1,2].set_ylabel('a')
fig[1,2].set_xlabel('d')
fig[2,0].scatter(data[:,7] , data[:,8], 1)
fig[2,0].set_ylabel('b')
fig[2,0].set_xlabel('c')
fig[2,1].scatter(data[:,7] , data[:,9], 1)
fig[2,1].set_ylabel('b')
fig[2,1].set_xlabel('d')
fig[2,2].scatter(data[:,8] , data[:,9], 1)
fig[2,2].set_ylabel('c')
fig[2,2].set_xlabel('a')
fig[2,2].tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
fig[0,2].axis('off')



plt.savefig("final.pdf",format="pdf")
plt.close()