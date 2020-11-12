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
    equil_t = []

    for i in range(500):
        print (i)
        g = Game_Of_Life(n, 'r')
        counter = 0
        for j in range(5000):
            init_sites = g.active_states()
            g.sweep()
            current_sites = g.active_states()
            if init_sites == current_sites:
                counter += 1
                if counter == 5:
                    equil_t.append(j-4)
                    break
            else:
                counter = 0

    np.savetxt('GOL_equilibration_t_histo.dat', np.column_stack([equil_t]))

def com(n, nsweeps):
    x = []
    y = []
    t = []
    vx = []
    vy = []
    vt = []
    velocity = []
    boundaries = []
    g = Game_Of_Life(n, 'g')
    for i in range (nsweeps):
        x1, y1 = g.centre_mass()
        if x1 > 0 and y1 > 0:
            x.append(x1)
            y.append(y1)
            t.append(i)
        if i > 320 and i < 420:
            vx.append(x1)
            vy.append(y1)
            vt.append(i)
        g.sweep()

    np.savetxt('GOL_com_coords.dat', np.column_stack([x, y, t]))
    xslope, xintercept = np.polyfit(vt, vx, 1)
    yslope, yintercept = np.polyfit(vt, vy, 1)

    vel = math.sqrt(xslope**2 + yslope**2)
    velocity.append(vel)
    np.savetxt('GOL_com_velocity.dat', velocity)

    plt.scatter(t, x, s=3, label = 'x radius')
    plt.scatter(t, y, s=3, label = 'y radius')
    plt.legend()
    plt.ylabel('X and Y radii')
    plt.xlabel('Time (Sweeps)')
    plt.title('X and Y Radii vs Time')
    plt.savefig('GOL_com_vs_time.png')
    plt.show()

def main():

    n = int(sys.argv[1])
    nsweeps = int(sys.argv[2])
    initial_condition = sys.argv[3]

    animation(n, nsweeps, initial_condition)
    #equil_histo(n, nsweeps)
    #com(n, nsweeps)

main()