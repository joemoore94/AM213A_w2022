import numpy as np
import matplotlib.pyplot as plt

# n is the number of images you wish to view at once
# n = 8
# k = int(n/2)
#
# fig, ax = plt.subplots(2, k)
# plt.gray()
# x = np.zeros((n,3355,5295))
#
# for i in range(n):
#     p = 20*2**i
#     File_data = np.loadtxt(f"Image_{p}.dat")
#     x[i] = File_data.astype(int)
#     ax[int(i/k)][i%k].imshow(x[i])
#
# plt.show()

fig, ax = plt.subplots()
plt.gray()

p = 3355

File_data = np.loadtxt(f"Image_{p}.dat")
x = File_data.astype(int)
ax.imshow(x)

ax.set_xticklabels([])
ax.set_yticklabels([])

ax.set_xticks([])
ax.set_yticks([])

plt.show()
