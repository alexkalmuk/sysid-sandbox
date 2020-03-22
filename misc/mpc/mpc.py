import numpy as np
from scipy import signal

A = np.array([[-0.0285, -0.0014], [-0.0371, -0.1476]])
B = np.array([[-0.0850, 0.0238], [0.0802, 0.4462]])
C = np.array([[0, 1], [1, 0]])
D = np.array([[0, 0], [0, 0]])

sys = signal.StateSpace(A, B, C, D)
print(sys)
