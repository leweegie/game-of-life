import math
import numpy as np
import scipy
from numpy.random import rand
import matplotlib.pyplot as plt

class Sirs(object):

    def __init__(self, n, p1, p2, p3):
        self.n = n
        self.state = self.random_lattice(self.n)
        self.p1 = float(p1)
        self.p2 = float(p2)
        self.p3 = float(p3)

    def random_lattice(self, n):
        lattice = np.random.randint(3, size = (n, n)) - 1
        return lattice

    def immune_lattice(self, p4):
        non_immune_p = (1 - p4)/3
        self.state = np.random.choice(4, size = (self.n, self.n), p = [non_immune_p, non_immune_p, non_immune_p, p4]) - 1

    def nb(self, i, j):
        neighbours = []
        check = False

        top         = self.state[i, (j-1)%self.n]
        bottom      = self.state[i, (j+1)%self.n]
        left        = self.state[(i-1)%self.n, j]
        right       = self.state[(i+1)%self.n, j]

        neighbours.append(top)
        neighbours.append(bottom)
        neighbours.append(left)
        neighbours.append(right)

        for x in range (len(neighbours)):
            if neighbours[x] ==0:
                check = True
            
        return check

    def active_states(self):
        active = np.sum(self.state)
        return active

    def sweep(self):
        #-1 = susceptible 0 = infected 1 = recovered
        #print(self.p1, self.p2, self.p3)
        for x in range (self.n * self.n):
            i = np.random.randint(self.n)
            j = np.random.randint(self.n)
            
            neighbours = self.nb(i, j)
            counter = 0

            if self.state[i, j] == -1 and neighbours:
                random1 = float(rand(1))
                if self.p1 >= random1:
                        self.state[i, j] = 0

            elif self.state[i, j] == 0:
                random2 = float(rand(1))
                if self.p2 >= random2:
                    self.state[i, j] = 1
                
            elif self.state[i, j] == 1:
                random3 = float(rand(1))
                if self.p3 >= random3:
                    self.state[i, j] = -1

    def average_infected(self):
        counter = 0
        for i in range (self.n):
            for j in range (self.n):
                if self.state[i,j] == 0:
                    counter += 1
        return counter