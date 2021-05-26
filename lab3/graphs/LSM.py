import matplotlib.pyplot as plt 
import numpy as np


x = [1.,  1.9, 2.8, 3.7, 4.6, 5.5]
y = [2.4142, 1.0818, 0.50953, 0.11836, -0.24008, -0.66818]

def F1(x):
    return 2.575567 - 0.627578 * x

def F2(x):
    return 3.547567-1.398066 * x + 0.118537 * x * x

n = 6
x_start = 1.
x_step = .9
x1range = np.arange(x_start, x_start + n * x_step - x_step, 0.01)
x2range = np.arange(x_start, x_start + n * x_step - x_step, 0.01)

plt.scatter(x, y)
plt.plot(x1range, [F1(x) for x in x1range])
plt.plot(x2range, [F2(x) for x in x2range])
plt.grid() 
plt.show()