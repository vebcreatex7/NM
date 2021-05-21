
import matplotlib.pyplot as plt 
import numpy as np


y = [0, 0.36892, 0.85408, 1.7856, 6.3138]


table = np.loadtxt("../logs/log2.txt", dtype='f')

x_start = 0
x_step = 0.9

def S(x):
    index = int((x - x_start) / x_step)
    x = x - x_start - index * x_step
    coef = table[index]
    return coef[0] + coef[1] * x + coef[2] * x ** 2 + coef[3] * x ** 3
    

xrange = np.arange(x_start, x_start + len(y) * x_step - x_step, 0.01)

plt.plot(np.arange(x_start, x_start + len(y) * x_step, x_step), y, 'ro')
plt.plot(xrange, [S(x) for x in xrange])
plt.grid() 
plt.show()

