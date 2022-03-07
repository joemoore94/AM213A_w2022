import numpy as np

dd = np.loadtxt('../least_squares_data.dat')

A0 = np.ones(len(dd))
A1 = dd[:,0]
A = np.vstack((A0, A1))

# n-1 degree of polynomial
n = 3
for i in range(2, n):
    A = np.vstack((A, A1**i))

A = A.T
B = np.matmul(A.T,dd[:,1])
A = np.matmul(A.T,A)

f = open("../Amat.dat", "w")
g = open("../Bmat.dat", "w")

f.write(f"{n} {n}\n")
g.write(f"{n} 1\n")

for i in range(n):
    for j in range(n):
        f.write(f"{A[i,j]} ")
    f.write(f"\n")
    g.write(f"{B[i]}\n")
