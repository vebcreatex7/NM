import numpy as np
import matplotlib.pyplot as plt
import math



def f(x):
    return math.exp(x**2) + math.exp(x * math.sqrt(2)) + math.exp(-x * math.sqrt(2))

x_start = 0
x_step = 0.1
n = 11

table1 = np.loadtxt("/home/courage/sem6/NM/lab4/logs/log1.1.txt")
table2 = np.loadtxt("/home/courage/sem6/NM/lab4/logs/log1.2.txt")
table3 = np.loadtxt("/home/courage/sem6/NM/lab4/logs/log1.3.txt")



xrange = np.arange(x_start, 1, 0.001)


plt.plot(table1[0], table1[1])
plt.scatter(table1[0], table1[1])
plt.plot(table2[0], table2[1])
plt.scatter(table2[0], table2[1])
plt.plot(table3[0], table3[1])
plt.scatter(table3[0], table3[1])

plt.plot(xrange, [f(x) for x in xrange])

plt.grid() 
plt.show()