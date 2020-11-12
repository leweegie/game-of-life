import math
import numpy as np
import scipy
from numpy.random import rand
import matplotlib.pyplot as plt

class Game_Of_Life(object):

    def __init__(self, n, initial_condition):
        self.n = n
        if initial_condition == 'r':
            self.state = self.random_lattice(self.n)
        if initial_condition == 'o':
            self.state = self.oscillator(self.n)
        if initial_condition == 'g':
            self.state = self.glider(self.n)

    def random_lattice(self, n):
        lattice = np.random.randint(2, size = (n, n))
        return lattice

    def oscillator(self, n):
        lattice = np.zeros([n,n])
        lattice[n/2, n/2] = 1
        lattice[n/2, n/2 - 1] = 1
        lattice[n/2, n/2 + 1] = 1
        return lattice

    def glider(self, n):
        lattice = np.zeros([n, n])
        lattice[n/2, n/2 - 1] = 1
        lattice[n/2, n/2 + 1] = 1
        lattice[n/2 + 1, n/2] = 1
        lattice[n/2 - 1, n/2 + 1] = 1
        lattice[n/2 + 1, n/2 + 1] = 1
        return lattice

    def nb(self, i, j):
        neighbours = []

        top         = self.state[i, (j-1)%self.n]
        bottom      = self.state[i, (j+1)%self.n]
        left        = self.state[(i-1)%self.n, j]
        right       = self.state[(i+1)%self.n, j]
        top_left    = self.state[(i-1)%self.n, (j-1)%self.n]
        top_right   = self.state[(i+1)%self.n, (j-1)%self.n]
        bottom_left = self.state[(i-1)%self.n, (j+1)%self.n]
        bottom_right= self.state[(i+1)%self.n, (j+1)%self.n]

        neighbours.append(top)
        neighbours.append(bottom)
        neighbours.append(left)
        neighbours.append(right)
        neighbours.append(top_left)
        neighbours.append(bottom_left)
        neighbours.append(top_right)
        neighbours.append(bottom_right)

        return neighbours

    def active_states(self):
        active = np.sum(self.state)
        return active

    def sweep(self):
        #0 = dead 1 = alive
        dummy = np.zeros([self.n, self.n])
        for i in range (self.n):
            for j in range (self.n):

                counter = 0
                neighbours = self.nb(i,j)
                
                for x in range(len(neighbours)):
                    if neighbours[x] == 1:
                        counter+=1

                if self.state[i,j] == 0:
                    if counter == 3:
                        dummy[i,j] = 1

                if self.state[i, j] == 1:
                    if counter < 2:
                        dummy[i, j] = 0
                    if counter > 3:
                        dummy[i, j] = 0
                    if counter == 2:
                        dummy[i, j] = 1
                    if counter == 3:
                        dummy[i, j] = 1

        self.state = np.copy(dummy)

    def centre_mass(self):
        locations = np.where(self.state == 1)
        if (abs(np.max(locations[0]) - np.min(locations[0]))) > 3:
            return -1, -1
        elif (abs(np.max(locations[1]) - np.min(locations[1]))) > 3:
            return -1, -1
        else:
            x_sum = np.sum(locations[0])
            y_sum = np.sum(locations[1])
            x_center = x_sum / len(locations[0])
            y_center = y_sum / len(locations[1])
            return x_center, y_center