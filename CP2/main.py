import math
import numpy as np
import scipy
from numpy.random import rand
import matplotlib.pyplot as plt
from game_of_life import Game_Of_Life
import sys

def animation(n, nsweeps, initial_condition):

    g = Game_Of_Life(n, initial_condition)
    fig = plt.figure()
    im=plt.imshow(g.state, animated=True)

    for i in range (nsweeps):

        g.sweep() #perform flip
        plt.cla()
        im=plt.imshow(g.state, animated=True, vmin = 0, vmax = 1)
        plt.draw()
        plt.pause(0.0001)

def equil_histo(n, nsweeps):
    g = Game_Of_Life(n, 'r')
    equil_t = []

    for i in range(100):
        while True:
            for j in range(nsweeps):
                init_sites = g.active_states()
                g.sweep()

                current_sites = g.active_states()
                if init_sites == current_sites:
                    equil_t.append(i)
                    break

    plt.hist(equil_t)
    plt.show()

def main():

    n = int(sys.argv[1])
    nsweeps = int(sys.argv[2])
    initial_condition = sys.argv[3]

    #animation(n, nsweeps, initial_condition)
    equil_histo(n, nsweeps)

main()