import numpy as np
import matplotlib.pyplot as plt


fig, ax = plt.subplots()

x = [20,40,80,160,320,640,1280,2560]
x = np.array(x)
y = [0.003543,0.003166,0.002703,0.002156,0.001604,0.001173,0.000711,0.000187]
y = np.array(y)

plt.plot(x,y)

plt.xlabel("k")
plt.ylabel(r"$E_k$")

plt.show()
