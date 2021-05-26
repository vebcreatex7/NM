
import matplotlib.pyplot as plt 
import numpy as np

def f(x):
    return 0.51082 * (x - 2.) + 0.69315


n = 5
x_start = 1.
x_step = 0.5
xrange = np.arange(x_start, x_start + n * x_step - x_step, 0.01)

plt.plot((1, 1.5, 2, 2.5, 3), (0, 0.40547, 0.69315, 0.91629, 1.0986))
plt.scatter([1, 1.5, 2, 2.5, 3], [0, 0.40547, 0.69315, 0.91629, 1.0986])
plt.plot(xrange, [f(x) for x in xrange])
plt.grid() 
plt.show()