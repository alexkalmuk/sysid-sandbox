import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

A = np.array([[-0.0285, -0.0014], [-0.0371, -0.1476]])
B = np.array([[-0.0850, 0.0238], [0.0802, 0.4462]])
C = np.array([[0, 1], [1, 0]])
D = np.array([[0, 0], [0, 0]])

sys = signal.StateSpace(A, B, C, D)

T = 100
t = np.arange(T)
x0 = np.array([0.5, 0.5])
u = np.ones((T, 2))
x = np.ones((T, 2))

x[0] = x0

for i in range(0, T - 1):
	x[i + 1] = A @ x[i] + B @ u[i]

print(x[:,0])

plt.figure()
plt.subplot(2,1,1)
plt.plot(t, x[:,0])
#plt.subplot(2,1,2)
plt.show()
