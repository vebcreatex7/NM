import numpy as np
import matplotlib.pyplot as plt
import math

x_0 = 0.
x_1 = 1.

table1 = np.loadtxt("/home/courage/sem6/NM/CP/logs/log1.txt")
table2 = np.loadtxt("/home/courage/sem6/NM/CP/logs/log2.txt")

plt.plot(table1[0], table1[1], label='Shooting Method')
plt.plot(table2[0], table2[1], label='Finite Difference Method')

plt.legend()
plt.grid() 
plt.show()