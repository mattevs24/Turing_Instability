import numpy as np
from matplotlib import pyplot as plt

U = np.loadtxt("solution.dat")
assert U.shape[0]==2, "The solution file should contain exactly 2 rows"

n = U.shape[1]
x = np.arange(1,n+1)/n

plt.subplot(2,1,1)
plt.plot(x, U[0,:], color="blue", label="$u$")
plt.legend()

plt.subplot(2,1,2)
plt.plot(x, U[1,:], color="red", label="$v$")
plt.xlabel("x")
plt.legend()

plt.savefig("solution.pdf")

