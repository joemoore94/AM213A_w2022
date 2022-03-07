import numpy as np

dd = np.loadtxt('../least_squares_data.dat')

A0 = np.ones(len(dd))
A1 = dd[:,0]
A = np.vstack((A0, A1))

# n-1 degree of polynomial
n = 6
for i in range(2, n):
    A = np.vstack((A, A1**i))

A = A.T
B = dd[:,1]

f = open("../Amat.dat", "w")
g = open("../Bmat.dat", "w")

f.write(f"{len(dd)} {n}\n")
g.write(f"{len(dd)} 1\n")

for i in range(len(dd)):
    for j in range(n):
        f.write(f"{A[i,j]} ")
    f.write(f"\n")
    g.write(f"{B[i]}\n")
