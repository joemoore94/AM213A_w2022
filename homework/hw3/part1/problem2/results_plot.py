import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()

r = np.loadtxt('../Cmat.dat', skiprows=1)
dd = np.loadtxt('../least_squares_data.dat')

def fun(x):
    sum = 0
    for i, val in enumerate(r):
        sum = sum + val*x**i
    return sum

x = np.linspace(0, 1, 100)

plt.plot(dd[:,0], dd[:,1],'o-', label='data')
plt.plot(x, fun(x),'-', label='fitted curve')

# sum = 0
# for i in dd:
#     sum = sum + (fun(i[0]) - i[1])**2
# sum = np.sqrt(sum)
# print("2-norm error between the fitted curve and the data:", sum)

plt.legend(loc="upper right")
plt.show()
