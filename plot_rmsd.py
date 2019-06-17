import matplotlib.pyplot as plt 
plt.switch_backend('agg')
import numpy as np

rmsd = np.loadtxt('rmsd.xvg',comments='#@')
plt.plot(rmsd[:,0],rmsd[:,1])
plt.savefig('rmsd.png')
