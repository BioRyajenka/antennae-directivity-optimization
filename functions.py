import numpy as np

def f1(x):
    return np.sum(x**2)

def f2(x):
    return np.abs(x).sum() + np.abs(x).prod()

def f3(x):
    sum = 0
    for i in range(x.size):
        sum = sum + x[:i].sum() ** 2
    return sum

def f4(x):
    return np.abs(x).max()

def f5(x):
    sum = 0
    for i in range(x.size - 1):
        sum = sum + (100 * (x[i + 1] - x[i] ** 2) ** 2 + (x[i] - 1) ** 2)
    return sum

def f6(x):
    return ((x + 0.5) ** 2).sum()

def f7(x):
    sum = 0
    for i in range(x.size):
        sum = sum + x[i] ** 4 * i + np.random.uniform(0, 1)
    return sum

def f8(x):
    sum = 0
    for i in range(x.size):
        sum = sum - x[i] * np.sin(abs(x[i]))
    return sum + x.size * 418.98288727243369

def f9(x):
    return (x ** 2 - 10 * np.cos(2 * np.pi * x) + 10).sum()

def f10(x):
    return -20 * np.exp(-0.2 * (x ** 2).mean()) - np.exp(np.cos(2 * np.pi * x).mean()) + 20 + np.e

def f11(x):
    ans = 0
    for i in range(x.size):
        ans = ans * np.cos(x[i] / np.sqrt(i + 1))
    return (x ** 2).sum() / 4000 - ans + 1

def u(xi, a, k):
    if xi > a:
        return k * (xi - a) ** 4
    if xi < -a:
        return k * (-xi - a) ** 4
    return 0

def f12(x):
    def y(xi):
        return 1.0 + (xi + 1.0) / 4.0
    
    sumy = 0
    for i in range(x.size - 1):
        sumy = sumy + (y(x[i]) - 1.0) ** 2 * (1.0 + 10.0 * (np.sin(np.pi * y(x[i + 1]))) ** 2)
    
    sumu = 0
    for i in range(x.size):
        sumu = sumu + u(x[i], 10.0, 100.0)
        
    return np.pi / x.size * (10.0 * (np.sin(np.pi * y(x[0]))) ** 2 + sumy + (y(x[-1]) - 1.0) ** 2) + sumu

def f13(x):
    sumy = 0
    for i in range(x.size - 1):
        sumy = sumy + (x[i] - 1.0) ** 2 * (1.0 + (np.sin(3.0 * np.pi * x[i + 1])) ** 2)
        
    sumu = 0
    for i in range(x.size):
        sumu = sumu + u(x[i], 5.0, 100.0)
        
    return 0.1 * (np.sin(3.0 * np.pi * x[0]) ** 2 + sumy 
                  + (x[-1] - 1.0) ** 2 * (1.0 + np.sin(2.0 * np.pi * x[-1]) ** 2)) + sumu

lbound = np.array([-100, -10, -100, -100, -30, -100, -1.28, -500, -5.12, -32, -600, -50, -50])
rbound = np.array([100, 10, 100, 100, 30, 100, 1.28, 500, 5.12, 32, 600, 50, 50])

#plot
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def plot_f(n):
    func = globals()[f"f{n}"]
    lb = lbound[n - 1]
    rb = rbound[n - 1]

    xs = np.linspace(lb, rb, 100)
    ys = np.linspace(lb, rb, 100)

    Xs, Ys = np.meshgrid(xs, ys)

    nu = np.vectorize(lambda a, b: func(np.array([a, b])))(Xs, Ys)

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(Xs, Ys, nu)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()