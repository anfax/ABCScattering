import matplotlib.pyplot as plt 
import numpy as np 
import os 
fil = 'ReactionProbabilities.dat'
if os.path.exists(fil):
    dat = np.loadtxt(fil)
    plt.plot(dat[:,0],dat[:,1])
    plt.xlabel('$E$ (eV)')
    plt.ylabel('Probability')
    plt.savefig('fig/Pro.jpg',dpi=666)
    plt.close() 