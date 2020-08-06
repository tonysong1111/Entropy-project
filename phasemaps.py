"""
Tony Song, Phase maps: to be imported into entropy plot

Updated July 2, 2020
"""
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from sympy import *
import seaborn as sns
from scipy.optimize import curve_fit
from collections import Counter


def draw_n(width, n=500):
    """
    draws n phases from a specified interval (initialization)

    :param width: real, width of the interval
    :param left: real, left endpoint of the interval
    :param n: unsigned integer, number of objects (default 500 if not specified)
    :return: np array of phases
    """
    np.random.seed(3)
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
        return 1e-15
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
            res += np.log(sp[1] - sp[0] + 1e-15)
        elif i == len_sp - 1:
            res += np.log(sp[-1] - sp[-2] + 1e-15)
        else:
            res += np.log(closest(sp[i-1], sp[i], sp[i+1]))
    return res / len_sp + np.log(2*len_sp - 2) + np.euler_gamma


def evolve(f, phases, num_iter):
    entropy_list = np.zeros(num_iter+1)
    num_bodies = len(phases)
    phase_list = np.array([[0]*num_bodies])
    for i in range(num_iter+1):
        entropy_list[i] = entropy(phases)
        phases = iterate(f, phases)
        phase_list = np.append(phase_list, np.array([phases]), axis=0)
    entropy_list = np.delete(entropy_list, 0)
    phase_list = phase_list[1:-1].reshape((num_iter, num_bodies))
    return phase_list, entropy_list


def change(arr, before, after):
    return arr[after] - arr[before]


def average_entropy(arr):
    """
    compute the asymptotic entropy from the last 20% of iterations

    :param arr: num_seq by num_iter array of entropies
    :return: 1 by num_iter array of average entropies
    """
    averages = np.mean(arr, axis=0)
    asym_h = np.mean(averages[int(len(averages) * 0.8):])
    print(f"The asymptotic entropy is approximately {asym_h}")
    return averages, asym_h


def plot_scatter(arr, averages, num_iter, n_seq, ax):
    xs = np.linspace(1, num_iter, num=num_iter)
    colors = cm.winter(np.linspace(0, 1, num=n_seq))
    for i, c in enumerate(colors):
        ax.scatter(xs, arr[i], color=c, s=1, alpha=0.2)
    ax.plot(xs, averages, c='black')
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Entropy")


def plot_hist(q, ax, hist=True, kde=False, color=None, label=None):
    sns.distplot(q, hist=hist, kde=kde, kde_kws={'shade': False, 'linewidth': 1}, ax=ax,
                 color=color, label=label)


def plot_circle(phases, ax):
    xs = np.cos(phases)
    ys = np.sin(phases)

    xss = np.linspace(0, 2*np.pi, num=10000)
    ax.scatter(xs, ys, c='r', marker='o', alpha = 0.1, zorder=10, s=3)
    ax.plot(np.cos(xss), np.sin(xss), color='black', lw=1, zorder=0)
    ax.set_aspect('equal', 'box')
    return xs, ys


def phase_evolution(f, phases, num_iter):
    phase_list, _ = evolve(f, phases, num_iter)
    return phase_list


def load_results(f, phases, num_iter, n_seq):
    ys = []
    body_phases = []

    for i in range(n_seq):
        phase_list, res = evolve(f, phases, num_iter)
        ys.append(res)
        body_phases.append(phase_list)
    ys = np.array(ys).reshape((n_seq, num_iter))
    body_phases = np.array(body_phases)
    return ys, body_phases


def plot_phases(f, phases, num_iter, n_seq, ax):
    """
    loads data for n_seq of time sequences
    this function plots only one of those time sequences
    """
    # load the data
    ys, body_phase = load_results(f, phases, num_iter, n_seq)
    selection = [int(num_iter * (1 - 0.05 * i)) - 1 for i in range(10)]

    # plot one trajectory and select phase distributions
    iteration = int(np.random.uniform(1, n_seq, 1))
    xs = np.linspace(1, num_iter, num_iter)
    count_distinct_phases = [count_distinct(body_phase[iteration][i]) for i in range(num_iter)]

    s = count_distinct_phases[0]
    for i in range(len(count_distinct_phases)):
        if i == 0:
            pass
        elif count_distinct_phases[i] == s:
            count_distinct_phases[i] = -9999
        else:
            s = count_distinct_phases[i]

    ax.plot(xs, ys[iteration], '-o', c='b', alpha=0.3, linewidth=1,
            markersize=2, zorder=0)
    ax.scatter([xs[i] for i in selection], [ys[iteration, i] for i in selection],
               color='r', s=2, zorder=5)
    for i, count in enumerate(count_distinct_phases):
        if count == -9999:
            pass
        else:
            ax.annotate(count, (xs[i], ys[iteration][i]), fontsize=4)

    # ax_inset1 = ax.inset_axes([0.15, 0.25, 0.15, 0.15], zorder=10)
    # ax_inset2 = ax.inset_axes([0.30, 0.25, 0.15, 0.15], zorder=10)
    # ax_inset3 = ax.inset_axes([0.45, 0.25, 0.15, 0.15], zorder=10)
    # ax_inset4 = ax.inset_axes([0.60, 0.25, 0.15, 0.15], zorder=10)
    # ax_inset5 = ax.inset_axes([0.75, 0.25, 0.15, 0.15], zorder=10)
    # ax_inset6 = ax.inset_axes([0.15, 0.05, 0.15, 0.15], zorder=10)
    # ax_inset7 = ax.inset_axes([0.30, 0.05, 0.15, 0.15], zorder=10)
    # ax_inset8 = ax.inset_axes([0.45, 0.05, 0.15, 0.15], zorder=10)
    # ax_inset9 = ax.inset_axes([0.60, 0.05, 0.15, 0.15], zorder=10)
    # ax_inset10 = ax.inset_axes([0.75, 0.05, 0.15, 0.15], zorder=10)
    # inset_dict = {'inset1': ax_inset1, 'inset2': ax_inset2, 'inset3': ax_inset3, 'inset4': ax_inset4,
    #               'inset5': ax_inset5, 'inset6': ax_inset6, 'inset7': ax_inset7, 'inset8': ax_inset8,
    #               'inset9': ax_inset9, 'inset10': ax_inset10}
    # inset_dict_keys = ['inset1', 'inset2', 'inset3', 'inset4', 'inset5', 'inset6', 'inset7',
    #                    'inset8', 'inset9', 'inset10']
    #
    # for i, key in enumerate(inset_dict_keys):
    #     plot_circle(body_phase[iteration, selection[i]], inset_dict[key])
    #     inset_dict[key].text(0, 0, f'n = {selection[i] + 1}', fontsize=4.5,
    #                          ha='center', va='center')
    #     inset_dict[key].axis('off')


def plot_results(f, phases, num_iter, n_seq):
    """
    plots evolution of entropy for an ensemble of objects

    :param f: function, phase map of one variable
    :param phases: (n,) np array of real-valued phases
    :param num_iter: unsigned int, number of iterations
    :param n_seq: unsigned int, number of time-sequences to average over
    :return: void
    """
    # set up the axes
    fig, ax1 = plt.subplots(nrows=2, ncols=3, figsize=(9,6))

    # load entropy evolution data
    ys, body_phases = load_results(f, phases, num_iter, n_seq)
    xs = np.linspace(1, num_iter, num=num_iter)
    ys_t = ys.T

    # compute the asymptotic entropy from the last 20% iterations
    averages, asym_h = average_entropy(ys)

    # plot scatter of trajectories and average
    plot_scatter(ys, averages, num_iter, n_seq, ax1[0][0])

    # # plot absolute entropy
    # abs_entropy = abs(ys_t[-1])
    # plot_hist(abs_entropy, ax1[0][1], True, True)
    # ax1[0][1].axvline(abs(asym_h), color='r', lw=1)
    # ax1[0][1].set_xlabel('absolute entropy')
    #
    # # plot the change in entropy from 2nd to last to last iteration
    # entropy_change_last = change(ys_t, -1, -2)
    # plot_hist(entropy_change_last, ax1[0][2], True, True)
    # ax1[0][2].set_xlabel('change in entropy')
    #
    # # plot lambda vs count
    # xss = np.linspace(0, 2*np.pi, num=10000)
    # slopes = local_slope(f, xss)
    # plot_hist(slopes, ax1[1][0], True, True)
    # ax1[1][0].set_xlabel('local lyapunov exponent')
    #
    # # entropy vs. change in entropy at the final stage
    # indices = [-1, -31, -61]
    # scatter_color = ['black', 'red', 'blue']
    # for i, c in zip(indices, scatter_color):
    #     ax1[1][1].scatter(ys_t[i-1], change(ys_t, i-1, i), color=c, s=1, label=f"n={int(xs[i])}")
    # ax1[1][1].legend(loc=0, prop={'size': 5})
    # ax1[1][1].set_xlabel('entropy $(n-1)^{th}$ iteration')
    # ax1[1][1].set_ylabel('change in entropy')

    # plot just one circle from the data
    #plot_circle(body_phases[10, -10], ax1[1][2])

    # delete extra set of axes
    fig.delaxes(ax1[1][2])

    plt.tight_layout()
    plt.show()
    return asym_h


def local_slope(f, ps):
    x = Symbol('x')
    df = f(x).diff(x)
    df = lambdify(x, df)
    return [np.log(abs(df(p))) for p in ps]


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
    :return: tuple of (a, b, c)
    """
    rows = ()
    col = np.zeros(3)
    ps.append((2*np.pi, 2*np.pi))
    for i, p in enumerate(ps):
        rows += (gen(p[0]),)
        col[i] = p[1]
    matrix = np.concatenate(rows, axis=0)
    coeffs = np.linalg.solve(matrix, col)
    return coeffs[0], coeffs[1], coeffs[2]


def plot_phasemaps(cubic_set, colors):
    fig, ax_main = plt.subplots()
    ax = ax_main.inset_axes([.1, .7, .3, .3])
    len_cubic_set = len(cubic_set)

    # plot phase maps
    for i in range(len_cubic_set):
        cubic_set[i].show(ax_main, colors[i])
    ax_main.set_xlabel(r"$\phi$")
    ax_main.set_ylabel(r"$\psi(\phi)$")

    # plot slope histograms as inset to phase maps
    xss = np.linspace(0, 2*np.pi, num=10000)
    for i in range(len(cubic_set)):
        slopes = local_slope(cubic_set[i].set_cubic, xss)
        avg_lyapunov = round(cubic_set[i].lyapunov(), 3)
        label = f"{avg_lyapunov}"
        plot_hist(slopes, ax, False, True, colors[i], label)
        #ax.axvline(avg_lyapunov, color=colors[i], lw=1)
    ax.set_xlabel(r"local lyapunov exponent $\lambda$", fontsize=5)
    ax.tick_params(axis='x', labelsize=5)
    ax.tick_params(axis='y', labelsize=5)
    ax.legend(title=r'$\langle\lambda\rangle$', prop={'size': 5})
    plt.show()


def count_distinct(arr):
    """

    :param arr: array of phases
    :return: integer, number of distinct values
    """
    return len(Counter(arr).keys())


def plot_decay(f, step_size=0.0001):
    xs = np.linspace(0, 2*np.pi, num=100000)

    left = -10
    right = 3
    lambdas = np.linspace(left, right, num=int((right - left)/step_size))
    slopes = local_slope(f, xs)
    density = np.zeros(len(lambdas) - 1)
    for i in range(len(lambdas) - 1):
        count = 0
        for slope in slopes:
            if lambdas[i + 1] >= slope >= lambdas[i]:
                count += 1
        density[i] = count

    fig, ax = plt.subplots()
    ax.plot(lambdas[:-1], density, c='b', lw=1)
    plt.show()


def inv_sq_fit(x, a, b):
    return a / (x**2) + b


def inv_sq_decrease(x, e):
    return 1/(x**2) - 1/(x-e)**2

def plot_decay2(f):
    fig, ax = plt.subplots()
    xs = np.linspace(0, 2*np.pi, num=100000)
    slopes = local_slope(f, xs)
    xs, ys = sns.distplot(slopes, hist=False, kde=True,
                          kde_kws={'shade': False, 'linewidth': 2},
                          color='b').get_lines()[0].get_data()

    peak = np.argmax(ys)
    xs_d = []
    ys_d = []
    increment = 0.05
    x = xs[peak-10]
    y = ys[peak-10]
    xs_d = []
    ys_d = []
    for i in range(10):
        y -= inv_sq_decrease(x, increment)
        x -= increment
        xs_d.append(x)
        ys_d.append(y)
    ax.scatter(xs_d, ys_d, color='r', s=15, zorder=10)
    # params, params_cov = curve_fit(inv_sq_fit, x_data, y_data)
    # ax.plot(x_data, inv_sq_fit(x_data, params[0], params[1]))
    plt.show()


class Cubic:

    xs = np.linspace(0, 2*np.pi, num=10000)
    coefficients = (0, 0, 0, 0)

    def __init__(self, f):
        self.f = f

    def set_coefficients(self, a, b, c, d=0):
        self.coefficients = (a, b, c, d)

    def print_coefficients(self):
        print(self.coefficients)

    def set_cubic(self, x):
        return self.f(coeffs=self.coefficients, x=x)

    def lyapunov(self):
        x = Symbol('x')
        f = self.set_cubic(x)
        df = f.diff(x)
        df = lambdify(x, df)
        lyap = sum(np.log(abs(df(self.xs)))) * (1 / len(self.xs))
        return lyap

    def show(self, ax, color='b'):
        f = self.set_cubic(self.xs)
        ax.plot(self.xs, f, color=color)
