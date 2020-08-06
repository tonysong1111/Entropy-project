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


def set_triangle(val):
    def trianglemap(x):
        res = []
        for phase in x:
            if phase <= np.pi:
                res.append(phase)
            elif phase <= val:
                res.append(-phase + 2*np.pi)
            else:
                m3 = val / (2*np.pi - val)
                b3 = 2*np.pi - (2*np.pi*val / (2*np.pi - val))
                res.append(m3*phase + b3)
        return res
    return trianglemap


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
    ### define model specifications from user###
    # p_list = prompt_cubic()
    # coeffs = pm.generate_cubic(p_list)
    # num_obj, num_iter, num_seq, width = prompt_config()

    ### hard code values for checking
    p_list = [(3.141, 2.356), (4.712, 0.5)]
    #p_list = [(3.141, 2.356), (4.712, 1.385)]
    #coeffs = pm.generate_cubic(p_list)

    num_obj, num_iter, num_seq, width = (300, 300, 50, 0.0001)

    ## define cubic function ###
    #cubic = pm.Cubic(cubicmap)
    #cubic.set_coefficients(coeffs[0], coeffs[1], coeffs[2])

    # show cubic function
    # fig, ax = plt.subplots(ncols=2, nrows=2)
    # x1, y1 = p_list[0]
    # x2, y2 = p_list[1]
    #
    # increment = 0.55
    # sub_axs = [ax[0][0], ax[0][1], ax[1][0], ax[1][1]]
    # for i, sub_ax in enumerate(sub_axs):
    #     p_list = [(x1, y1 - i * increment), (x2, y2 + i * increment)]
    #     coeffs = pm.generate_cubic(p_list)
    #     cubic = pm.Cubic(cubicmap)
    #     cubic.set_coefficients(coeffs[0], coeffs[1], coeffs[2])
    #     cubic.show(sub_ax)
    #     sub_ax.plot([p_list[0][0], p_list[1][0]], [p_list[0][1], p_list[1][1]], 'ro')
    #     label = rf"$\langle\lambda\rangle={round(cubic.lyapunov(),3)}$"
    #     sub_ax.text(0.45, 0.8, label, fontsize=10, transform=sub_ax.transAxes)
    # plt.savefig("graphics/phasemaps_four.png", bbox_inches='tight', dpi=1000)
    # plt.show()

    ## print the lyapunov exponent ###
    #fstring = f"The average Lyapunov exponent is {cubic.lyapunov()}"
    #print(fstring)

    ### plot results ###
    # initialize some randome phases according to parameters
    random = pm.draw_n(width, num_obj)

    # make one plot (regular)
    # fig, ax = plt.subplots()
    # ys, _ = pm.load_results(cubic.set_cubic, random, num_iter, num_seq)
    # averages, avg_H = pm.average_entropy(ys    fig, ax = plt.subplots()
    # pm.plot_scatter(ys, averages, num_iter, num_seq, ax)
    # ax.text(0.5, 0.17, rf"$\langle\lambda\rangle={round(cubic.lyapunov(), 3)}$", fontsize=15,
    #         transform=ax.transAxes)
    # plt.savefig("graphics/lyap0.141_fullstart.png", bbox_inches='tight', dpi=1000)
    # plt.show()


    ### plot for same phase map but different starting points
    # fig, ax = plt.subplots()
    # xs = np.linspace(1, num_iter, num=num_iter)
    #
    # objects_full = pm.draw_n(0.0001, num_obj)
    # objects_concentrated = pm.draw_n(2*np.pi, num_obj)
    # ys1, body_phases1 = pm.load_results(cubic.set_cubic, objects_full,
    #                                     num_iter, num_seq)
    # ys2, body_phases2 = pm.load_results(cubic.set_cubic, objects_concentrated,
    #                                     num_iter, num_seq)
    # averages1, asym_H1 = pm.average_entropy(ys1)
    # averages2, asym_H2 = pm.average_entropy(ys2)
    # ax.plot(xs, averages1, c='b', label=r'$w=$1e-4')
    # ax.plot(xs, averages2, c='r', label=r'$w=2\pi$')
    # ax.text(0.65, 0.68, rf'$\langle H\rangle_\infty={round(asym_H2, 3)}$', color='r',
    #         transform=ax.transAxes)
    # ax.text(0.65, 0.475, rf'$\langle H\rangle_\infty={round(asym_H1, 3)}$', color='b',
    #         transform=ax.transAxes)
    # plt.legend()
    # #plt.savefig("graphics/different_starting.png", bbox_inches='tight', dpi=1000)
    # plt.show()




    ### plot for triangle phase map
    # phasemap1 = set_triangle(3.669)
    # phasemap2 = set_triangle(5.989)
    # fig, ax = plt.subplots(nrows=2, ncols=2)
    #
    # ys1, body_phases1 = pm.load_results(phasemap1, random, num_iter, num_seq)
    # averages1, asym_H1 = pm.average_entropy(ys1)
    # ys2, body_phases2 = pm.load_results(phasemap2, random, num_iter, num_seq)
    # averages2, asym_H2 = pm.average_entropy(ys2)
    #
    # xs = np.linspace(0, 2*np.pi, num=10000)
    #
    # ax[0][0].plot(xs, phasemap1(xs), c='b', lw=2)
    # ax[0][0].text(0.1, 0.9, r"$\psi(\phi)$ with $\langle\lambda\rangle=0.141$",
    #               transform=ax[0][0].transAxes, size=8)
    # pm.plot_scatter(ys1, averages1, num_iter, num_seq, ax[0][1])
    # ax[0][1].text(0.7, 0.1, rf"$\langle H\rangle_\infty=${round(asym_H1,3)}",
    #               transform=ax[0][1].transAxes, size=8)
    #
    # ax[1][0].plot(xs, phasemap2(xs), c='b', lw=2)
    # pm.plot_scatter(ys2, averages2, num_iter, num_seq, ax[1][1])
    # ax[1][1].text(0.7, 0.1, rf"$\langle H\rangle_\infty=${round(asym_H2, 3)}",
    #               transform=ax[1][1].transAxes, size=8)
    #
    # for sub_ax in [ax[0][0], ax[1][0]]:
    #     sub_ax.set_aspect("equal")
    #
    # fig.tight_layout()
    # #plt.savefig("graphics/trianglemap_comparison.png", bbox_inches='tight', dpi=1000)
    # plt.show()

    #lambda vs entropy
    # fig, ax = plt.subplots()
    # num = 350
    # lambdas = [0]*num
    # asym_hs = [0]*num
    # increment = 0.005
    # x1, y1 = p_list[0]
    # x2, y2 = p_list[1]
    # cubic = pm.Cubic(cubicmap)
    # for j in range(num):
    #     p_list = [(x1, y1 - j * increment), (x2, y2 + j * increment)]
    #     coeffs = pm.generate_cubic(p_list)
    #     cubic = pm.Cubic(cubicmap)
    #     cubic.set_coefficients(coeffs[0], coeffs[1], coeffs[2])
    #     lambdas[j] = cubic.lyapunov()
    #     ys, _ = pm.load_results(cubic.set_cubic, random, num_iter, num_seq)
    #     _, asym_hs[j] = pm.average_entropy(ys)
    #     print(j)
    # ax.scatter(lambdas, asym_hs, color='b', s=5)
    # ax.set_xlabel(r"$\langle\lambda\rangle$")
    # ax.set_ylabel(r"$\langle H\rangle_\infty$")
    # plt.savefig("graphics/big_lambdaVSentropy_300iter.png", bbox_inches='tight', dpi=1000)
    # plt.show()


    #pm.plot_results(cubic.set_cubic, random, num_iter, num_seq)

    # plot phase circles for individual trajectories
    # p_list = [(3.141, 2.356), (4.712, 1.385+0.25)]
    # coeffs = pm.generate_cubic(p_list)
    #
    # cubic = pm.Cubic(cubicmap)
    # cubic.set_coefficients(coeffs[0], coeffs[1], coeffs[2])
    #
    fig, ax = plt.subplots(ncols=2, nrows=2)
    ax_list = [ax[0][0], ax[0][1], ax[1][0], ax[1][1]]
    np.random.seed(1)
    x1, y1 = p_list[0]
    x2, y2 = p_list[1]

    increment = 0.3
    sub_axs = [ax[0][0], ax[0][1], ax[1][0], ax[1][1]]
    for i, sub_ax in enumerate(sub_axs):
        p_list = [(x1, y1 - i * increment), (x2, y2 + i * increment)]
        coeffs = pm.generate_cubic(p_list)
        cubic = pm.Cubic(cubicmap)
        cubic.set_coefficients(coeffs[0], coeffs[1], coeffs[2])
        pm.plot_phases(cubic.set_cubic, random, num_iter, num_seq, sub_ax)
        sub_ax.text(0.45, 0.1, rf"$\langle\lambda\rangle={round(cubic.lyapunov(), 3)}$",
                    fontsize=10, transform=sub_ax.transAxes)
    #plt.suptitle(rf"$\langle\lambda\rangle=${round(cubic.lyapunov(), 3)}")
    plt.savefig("graphics/lostpoints_manylyapunov.png", bbox_inches='tight', dpi=1000)
    plt.show()

    # set up a set of cubic functions with lyapunov exponents in multiple regimes
    # len_cubic_set = 4
    # colors = ['black', 'r', 'b', 'g']
    # cubic_set = [pm.Cubic(cubicmap) for i in range(len_cubic_set)]
    # p_lists = []
    # for i in range(4):
    #     x2, y2 = p_list[1]
    #     p_lists.append([p_list[0], (x2, y2 + 0.25*i)])
    # coeffss = [pm.generate_cubic(p) for p in p_lists]
    # for i, cubic in enumerate(cubic_set):
    #     cubic.set_coefficients(coeffss[i][0], coeffss[i][1], coeffss[i][2])

    # plot phase maps and their local lyapunov distributions as inset
    #pm.plot_phasemaps(cubic_set, colors)

    # plot results for cubic_set
    #for i in range(4):
    #    pm.plotresults(cubic_set[i].set_cubic, random, num_iter, num_seq)

    #hist_vals = [0]*100

    #for i in range(100):
    #    hist_vals[i] = pm.plotresults(cubic.set_cubic, random, num_iter, num_seq)

    #pm.plt.hist(hist_vals)
    #pm.plt.show()

    # plot phases distributed over a circle
    #pm.plot_circle(np.random.uniform(0, 1, 300))

    #pm.plot_decay2(cubic.set_cubic)

    #see what happens starting from uniform
    # ncols = 10
    # nrows = 10
    # phase_list, entropy_list = pm.evolve(phasemap1, random, ncols*nrows)
    # #print(entropy_list)
    # #print(coeffs)
    # fig, ax = plt.subplots(nrows=nrows, ncols=ncols)
    # for i in range(nrows):
    #     for j in range(ncols):
    #         pm.plot_circle(phase_list[ncols*i + j], ax[i][j])
    #         ax[i][j].axis('off')
    #         ax[i][j].text(0, 0, f'n = {ncols*i + j + 1}', fontsize=4.5,
    #                       ha='center', va='center')
    # #r_coeffs = [round(c, 3) for c in coeffs]
    # #title = rf"$\psi(\phi)=${r_coeffs[0]}$\phi^3-${abs(r_coeffs[1])}$\phi^2+${r_coeffs[2]}$\phi,\quad\lambda=${round(cubic.lyapunov(), 3)}"
    # title = r"triangle map with $\langle\lambda\rangle=0.141$ (top row map)"
    # plt.suptitle(title)
    # #plt.savefig('graphics/trianglemap_fig100.png', bbox_inches='tight', dpi=1000)
    # plt.show()