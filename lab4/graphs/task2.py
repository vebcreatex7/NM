import numpy as np
import matplotlib.pyplot as plt
import math



def f(x):
    return 1. / x + 1.

x_start = 1
x_step = 0.1
n = 11

table1 = np.loadtxt("/home/courage/sem6/NM/lab4/logs/log2.1.txt")



xrange = np.arange(x_start, 2, 0.001)


plt.plot(table1[0], table1[1])
plt.scatter(table1[0], table1[1])

plt.plot(xrange, [f(x) for x in xrange])

plt.grid() 
plt.show()