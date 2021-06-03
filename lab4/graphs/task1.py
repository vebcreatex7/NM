import numpy as np
import matplotlib.pyplot as plt
import math



def f(x):
    return math.exp(x**2) + math.exp(x * math.sqrt(2)) + math.exp(-x * math.sqrt(2))

x_0 = 0.
x_1 = 1.



table1 = np.loadtxt("/home/courage/sem6/NM/lab4/logs/log1.1.txt")
table2 = np.loadtxt("/home/courage/sem6/NM/lab4/logs/log1.2.txt")
table3 = np.loadtxt("/home/courage/sem6/NM/lab4/logs/log1.3.txt")


plt.plot(table1[0], table1[1])
plt.scatter(table1[0], table1[1])
plt.plot(table2[0], table2[1])
plt.scatter(table2[0], table2[1])
plt.plot(table3[0], table3[1])
plt.scatter(table3[0], table3[1])

xrange = np.arange(x_0, x_1, 0.001)
plt.plot(xrange, [f(x) for x in xrange])

plt.grid() 
plt.show()