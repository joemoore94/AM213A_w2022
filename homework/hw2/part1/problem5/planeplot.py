import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math

a,b,c = 0.039034, 0.363382, 0.078067
d = 0.039034 + 2*0.363382 + 3*0.078067

x = np.linspace(-5,5,10)
y = np.linspace(-5,5,10)

X,Y = np.meshgrid(x,y)
Z = (d - a*X - b*Y) / c

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(X, Y, Z, color='black', linewidth=0, alpha=0.5)

# Original points
ax.scatter([1, -3, math.pi], [2, 2, math.e], [3, 5, -math.sqrt(2)], color='r', s=20)

plt.show()
