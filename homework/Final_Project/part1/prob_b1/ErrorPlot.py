import numpy as np
import matplotlib.pyplot as plt


fig, ax = plt.subplots()

data = np.loadtxt(f"JacobiData.dat")
plt.plot(data[:,0],data[:,1], label='Jacobi')

data = np.loadtxt(f"SeidelData.dat")
plt.plot(data[:,0],data[:,1], label='Seidel')

plt.xlabel("Interations")
plt.ylabel("Error")
ax.set_yscale('log')
ax.legend()
plt.title("D = 1000")
plt.show()
