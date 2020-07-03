"""
Tony Song, Entropy project

Updated: July 2nd, 2020
"""
import phasemaps as pm


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
    if len(p_list) != 3:
        print("you must enter exactly three coordinate pairs")
        return -1
    return p_list


def prompt_config():
    width = float(input("Enter width of initial interval: "))
    num_obj = int(input("Enter number of objects: "))
    num_iter = int(input("Enter number of iterations: "))
    num_seq = int(input("Enter number of time sequences: "))
    return num_obj, num_iter, num_seq, width


if __name__ == "__main__":
    p_list = prompt_cubic()
    coeffs = pm.generate_cubic(p_list)
    num_obj, num_iter, num_seq, width = prompt_config()

    ### define cubic function ###
    cubic = pm.Cubic(cubicmap)
    cubic.set_coefficients(coeffs[0], coeffs[1], coeffs[2])

    ### print the lyapunov exponent ###
    fstring = f"The average Lyapunov exponent is {cubic.lyapunov()}"
    print(fstring)

    ### plot results ###
    random = pm.draw_n(width, num_obj)
    pm.plotresults(cubic.set_cubic, random, num_iter, num_seq)

    #hist_vals = [0]*100

    #for i in range(100):
    #    hist_vals[i] = pm.plotresults(cubic.set_cubic, random, num_iter, num_seq)

    #pm.plt.hist(hist_vals)
    #pm.plt.show()

