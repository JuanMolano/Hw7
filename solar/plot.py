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
f,fig = plt.subplots(3,4,figsize=(12,10))

plt.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')

fig[0,0].scatter(data[:,0] , data[:,1], 0.01)
fig[0,0].plot(data[:,0] , data[:,7], color = 'red')
fig[0,0].tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
fig[0,0].set_title("$Fit$")
fig[1,0].plot( data[:,0] , data[:,2] , color = 'black')
fig[1,0].tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
fig[1,0].set_title("$Chi$")
fig[2,0].scatter( data[:,3] , data[:,2], 0.01)
fig[2,0].set_title("$a$")
fig[2,0].tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
fig[2,1].scatter(data[:,4] , data[:,2], 0.01)
fig[2,1].set_title("$b$")
fig[2,1].tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
fig[2,2].scatter( data[:,5] , data[:,2], 0.01)
fig[2,2].set_title("$c$")
fig[2,2].tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
fig[2,3].scatter(data[:,6] , data[:,2], 0.01)
fig[2,3].set_title("$d$")
fig[2,3].tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
fig[0,1].axis('off')
fig[0,2].axis('off')
fig[0,3].axis('off')
fig[1,1].axis('off')
fig[1,2].axis('off')
fig[1,3].axis('off')


plt.savefig("fit.pdf",format="pdf")
plt.close()