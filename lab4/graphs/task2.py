import numpy as np
import matplotlib.pyplot as plt
import math



def f(x):
    return 1. / x + 1.

x_0 = 1.
x_1 = 2.


table1 = np.loadtxt("/home/courage/sem6/NM/lab4/logs/log2.1.txt")
table2 = np.loadtxt("/home/courage/sem6/NM/lab4/logs/log2.2.txt")


plt.plot(table1[0], table1[1])
plt.scatter(table1[0], table1[1])
plt.plot(table2[0], table2[1])
plt.scatter(table2[0], table2[1])

xrange = np.arange(x_0, x_1, 0.001)
plt.plot(xrange, [f(x) for x in xrange])

plt.grid() 
plt.show()