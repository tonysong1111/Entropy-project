"""
Tony Song, Entropy project

Updated: July 2nd, 2020
"""
import phasemaps as pm
import numpy as np
import matplotlib.pyplot as plt


def cubicmap(coeffs, x):
    """
    :param coeffs: tuple of 4 real numbers for the cubic function
    :param x: real number variable
    :return: cubic function
    """
    a, b, c, d = coeffs
    return a*x**3 + b*x**2 + c*x + d


def prompt_cubic():
    p_list = []
    line = input("Enter a list of coordinates, in the form of x y:\n")
    while line != "":
        p_list.append(tuple(map(float, line.split())))
        line = input()
    if len(p_list) != 2:
        print("you must enter exactly two coordinate pairs")
        return -1
    return p_list


def prompt_config():
    width = float(input("Enter width of initial interval: "))
    num_obj = int(input("Enter number of objects: "))
    num_iter = int(input("Enter number of iterations: "))
    num_seq = int(input("Enter number of time sequences: "))
    return num_obj, num_iter, num_seq, width


if __name__ == "__main__":
    ### define model specifications ###
    # p_list = prompt_cubic()
    # coeffs = pm.generate_cubic(p_list)
    # num_obj, num_iter, num_seq, width = prompt_config()
    p_list = [(3.141, 2.356), (4.712, 1.385)]
    coeffs = pm.generate_cubic(p_list)
    num_obj, num_iter, num_seq, width = (300, 150, 50, 2*np.pi)

    ### define cubic function ###
    #cubic = pm.Cubic(cubicmap)
    #cubic.set_coefficients(coeffs[0], coeffs[1], coeffs[2])

    ### print the lyapunov exponent ###
    #fstring = f"The average Lyapunov exponent is {cubic.lyapunov()}"
    #print(fstring)

    ### plot results ###
    random = pm.draw_n(width, num_obj)
    #pm.plotresults(cubic.set_cubic, random, num_iter, num_seq)

    # set up a set of cubic functions with lyapunov exponents in multiple regimes
    len_cubic_set = 4
    colors = ['black', 'r', 'b', 'g']
    cubic_set = [pm.Cubic(cubicmap) for i in range(len_cubic_set)]
    p_lists = []
    for i in range(4):
        x2, y2 = p_list[1]
        p_lists.append([p_list[0], (x2, y2 + 0.25*i)])
    coeffss = [pm.generate_cubic(p) for p in p_lists]
    for i, cubic in enumerate(cubic_set):
        cubic.set_coefficients(coeffss[i][0], coeffss[i][1], coeffss[i][2])

    # plot phase maps and their local lyapunov distributions as inset
    pm.plot_phasemaps(cubic_set, colors)

    for i in range(4):
        pm.plotresults(cubic_set[i].set_cubic, random, num_iter, num_seq)

    #hist_vals = [0]*100

    #for i in range(100):
    #    hist_vals[i] = pm.plotresults(cubic.set_cubic, random, num_iter, num_seq)

    #pm.plt.hist(hist_vals)
    #pm.plt.show()

