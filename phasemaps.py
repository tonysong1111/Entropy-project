"""
Tony Song, Phase maps: to be imported into entropy plot

Updated July 2, 2020
"""
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from sympy import *


def draw_n(width, n=500):
    """
    draws n phases from a specified interval (initialization)

    :param width: real, width of the interval
    :param left: real, left endpoint of the interval
    :param n: unsigned integer, number of objects (default 500 if not specified)
    :return: np array of phases
    """
    np.random.seed(2)
    left = np.random.uniform(0, 2*np.pi, 1)
    arr = np.array(np.random.uniform(left, left+width, n))
    return np.mod(arr, 2*np.pi)


def iterate(f, phases):
    """
    maps ensemble of phases, then translate by random amount

    :param phases: (n,) np array, ensemble of real-valued phases between 0 and 1
    :return: another ensemble of real-valued phases
    """
    shift = np.random.uniform(0, 2*np.pi, 1)
    phases = f(phases)
    return np.mod(phases + shift, 2*np.pi)


def closest(i, j, k):
    """
    returns minimum difference between neighbor phases

    :param i: real
    :param j: real
    :param k: real
    :return: real, the smallest difference
    """
    d1 = j-i
    d2 = k-j
    if d1 == 0 or d2 == 0:
        return 1e-30
    elif d1 >= d2:
        return d2
    else:
        return d1


def entropy(phases):
    """
    returns the differential entropy of the ensemble

    :param phases: (n,) np array of real-valued phases between 0 and 1
    :return: real, entropy
    """
    sp = np.sort(phases)
    len_sp = len(sp)
    res = 0
    for i in range(len_sp):
        if i == 0:
            res += np.log(sp[1] - sp[0] + 1e-30)
        elif i == len_sp - 1:
            res += np.log(sp[-1] - sp[-2] + 1e-30)
        else:
            res += np.log(closest(sp[i-1], sp[i], sp[i+1]))
    return res / len_sp + np.log(2*len_sp - 2) + np.euler_gamma


def evolve(f, phases, num_iter):
    entropy_list = np.zeros(num_iter+1)
    for i in range(num_iter+1):
        entropy_list[i] = entropy(phases)
        phases = iterate(f, phases)
    entropy_list = np.delete(entropy_list, 0)
    return entropy_list


def plotresults(f, phases, num_iter, n_seq):
    """
    plots evolution of entropy for an ensemble of objects

    :param f: function, phase map of one variable
    :param phases: (n,) np array of real-valued phases
    :param num_iter: unsigned int, number of iterations
    :param n_seq: unsigned int, number of time-sequences to average over
    :return: void
    """
    fig, ax1 = plt.subplots()
    xs = np.linspace(1, num_iter, num=num_iter)
    ys = []
    for i in range(n_seq):
        res = evolve(f, phases, num_iter)
        ys.append(res)
    np.array(ys).reshape((num_iter, n_seq))
    averages = np.mean(ys, axis=0)
    asym_h = np.mean(averages[int(len(averages) * 0.8):])
    print(f"The asymptotic entropy is approximately {asym_h}")
    colors = cm.spring(np.linspace(0, 1, n_seq))
    for y, col in zip(ys, colors):
        plt.scatter(xs, y, s=1, color=col, alpha=0.5)
    plt.plot(xs, averages, c='black')
    ax1.set_xlabel('iteration')
    ax1.set_ylabel('entropy')
    plt.show()
    return asym_h


def gen(x_coord):
    """
    generate the cubed, squared, and first power of some value

    :param x_coord: real number
    :return: tuple of 3 real numbers
    """
    return np.array([[x_coord**3, x_coord**2, x_coord]])


def generate_cubic(ps):
    """
    generate the three coefficients for a cubic (d = 0)

    :param p1: coordinate pair (x, y)
    :param p2: ^
    :param p3: ^
    :return: tuple of (a, b, c)
    """
    rows = ()
    col = np.zeros(3)
    for i, p in enumerate(ps):
        rows += (gen(p[0]),)
        col[i] = p[1]
    matrix = np.concatenate(rows, axis=0)
    coeffs = np.linalg.solve(matrix, col)
    return coeffs[0], coeffs[1], coeffs[2]


class Cubic:

    xs = np.linspace(0, 2*np.pi, num=10000)
    coefficients = (0, 0, 0, 0)

    def __init__(self, f):
        self.f = f

    @classmethod
    def set_coefficients(cls, a, b, c, d=0):
        cls.coefficients = (a, b, c, d)

    def set_cubic(self, x):
        return self.f(coeffs=self.coefficients, x=x)

    def lyapunov(self):
        x = Symbol('x')
        f = self.set_cubic(x)
        df = f.diff(x)
        df = lambdify(x, df)
        lyap = sum(np.log(abs(df(self.xs)))) * (1 / len(self.xs))
        return lyap

    def show(self):
        f = self.set_cubic(self.xs)
        plt.plot(self.xs, f, color='b')
        plt.show()
