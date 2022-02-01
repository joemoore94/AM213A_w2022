import matplotlib.pyplot as plt
import numpy as np

def f(x):
    return (x-2)**9

def g(x):
    return x**9 - 18*x**8 + 144*x**7 - 672*x**6 + 2016*x**5 - 4032*x**4 \
            + 5376*x**3 - 4608*x**2 + 2304*x - 512

x = np.linspace(1.920, 2.08, 160)

fig, ax = plt.subplots()

ax.plot(x, f(x))
ax.plot(x, g(x))
plt.show()

"""
Part (c)
It is clear that f(x) is not as well-conditioned as I thought. I'm not sure why though.
"""
