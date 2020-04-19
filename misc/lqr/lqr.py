import numpy as np
import matplotlib.pyplot as plt
import control
from control.matlab import *

A = np.matrix([[-0.0285, -0.0014], [-0.0371, -0.1476]])
B = np.matrix([[-0.0850, 0.0238], [0.0802, 0.4462]])

C = ctrb(A, B)
print("Controllability matrix =", C)
print("Rank =", np.linalg.matrix_rank(C))

ss = StateSpace(A,B, np.eye(2), np.zeros((2,2)))
p = pole(ss)
print("poles =", p)

# Start with a diagonal weighting
Q = np.diag([1, 1])
R = np.diag([1, 1])
K, S, E = lqr(A, B, Q, R)

print("K =", K)

T = 100
t = np.arange(T)
x0 = np.array([1, 0.5])
x = np.zeros((T, 2))
x[0] = x0

for i in range(0, T - 1):
	x[i + 1] = A @ x[i] + (B @ K) @ x[i]

plt.figure()
plt.subplot(2,1,1)
plt.plot(t, x[:,0])
plt.show()
