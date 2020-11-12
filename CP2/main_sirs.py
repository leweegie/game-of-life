import math
import numpy as np
import scipy
from scipy import stats
from numpy.random import rand
import matplotlib.pyplot as plt
from sirs import Sirs
import sys
from astropy.stats import jackknife_resampling
from astropy.stats import jackknife_stats

def animation(n, nsweeps, p1, p2, p3):

    g = Sirs(n, p1, p2, p3)
    fig = plt.figure()
    im=plt.imshow(g.state, animated=True)

    for i in range (nsweeps):

        g.sweep() #perform flip
        plt.cla()
        im=plt.imshow(g.state, animated=True, vmin = -1, vmax = 1)
        plt.draw()
        plt.pause(0.0001)
    print(g.state)

def jackknife(data, n):
    resamples = jackknife_resampling(data)
    x_resamples = [np.var(resamples[i])/(n*n) for i in range (len(resamples))]
    v = np.var(data)/(n*n)
    return np.sqrt(np.sum([(x_resamples[i] - v) * (x_resamples[i]-v) for i in range (len(x_resamples))]))

def average_infected_sites(n):

    p1 = np.linspace(0,1,21)
    p2 = 0.5    
    p3 = np.linspace(0,1,21)
    fractions = np.zeros((21,21))
    variances = np.zeros((21,21))

    for i in range(len(p1)):
        print i
        for j in range(len(p3)):
            infected = []
            g = Sirs(n,p1[i], p2, p3[j])
            for x in range (1000):
                g.sweep()
                if x > 100:
                    no_inf = g.average_infected()
                    infected.append(no_inf)
                    
                    if no_inf == 0:
                        break
            if no_inf == 0:
                v = 0
                fraction = 0
            else:
                average = np.mean(np.asarray(infected))
                v = (np.var(np.asarray(infected))/(n*n))
                fraction = average / (n*n)
            variances[i][j] = v
            fractions[i][j] = fraction

    np.savetxt('average_infected.dat', fractions)
    np.savetxt('variance.dat', variances)

def p4_b(n):

    p1 = np.linspace(0.2,0.5,31)
    p2 = 0.5    
    p3 = 0.5
    fractions = []
    v_errors = []

    for i in range(len(p1)):
        print i
        infected = []
        g = Sirs(n, p1[i], p2, p3)
        for x in range (10000):
            g.sweep()
            if x > 100:
                inf_no = g.average_infected()
                infected.append(inf_no)
                if inf_no == 0:
                    break
        if inf_no == 0:
            v = 0
            v_e = 0
        else:
            v = np.var(np.array(infected)) / (n*n)
            v_e = jackknife(np.array(infected), n)
        fractions.append(v)
        v_errors.append(v_e)

    np.savetxt('p1_vs_v.dat', np.column_stack([p1, fractions, v_errors]))

def p5(n):

    p1 = 0.5
    p2 = 0.5
    p3 = 0.5
    p4 = np.linspace(0,1,101)
    fractions = np.zeros((101,5))
    means = []
    av_errors = []

    for z in range(5):
        for i in range(len(p4)):
            print i
            infected = []
            g = Sirs(n, p1, p2, p3)
            g.immune_lattice(p4[i])
            for x in range (1000):
                g.sweep()
                if x > 100:
                    inf_no = g.average_infected()
                    infected.append(inf_no)
                    if inf_no == 0:
                        break
            if inf_no == 0:
                av_inf = 0
                av_e = 0
            else:
                av_inf = np.mean(infected) / (n*n)
            fractions[i, z] = av_inf

    for i in range(len(p4)):
        mean = np.mean(fractions[i])
        error = scipy.stats.sem(fractions[i])
        means.append(mean)
        av_errors.append(error)

    np.savetxt('in_f_vs_av_inf.dat', np.column_stack([p4, means, av_errors]))

def plot_histogram():

    equilibration_times = np.loadtxt('GOL_equilibration_t_histo.dat')
    plt.hist(equilibration_times, 75)
    plt.ylabel('Bin Count')
    plt.xlabel('Number of Runs')
    plt.title('Equilibration times')
    plt.savefig('GOL_histogram.png')
    plt.show()

def plot_contour(filename, graphname, title):
    
    #fig = plt.figure()
    p1 = np.linspace(0,1,21)
    p3 = np.linspace(0,1,21)
    X,Y = np.meshgrid(p1,p3)
    fig, ax = plt.subplots(1,1)
    
    
    data = np.loadtxt(filename)
    data = np.transpose(data)
    
    cp = ax.contourf(X, Y, data)
    plt.colorbar(cp)
    plt.ylabel('p3')
    plt.xlabel('p1')
    plt.title(title)
    plt.savefig(graphname)
    plt.show()
        
def plot_heatmap(filename, graphname, title):

    fig, ax = plt.subplots(1,1)
    data = np.loadtxt(filename)
    data = np.transpose(data)
    im=plt.imshow(data, origin = 'lower', extent = [0,1,0,1])
    plt.colorbar(im)
    plt.ylabel('p3')
    plt.xlabel('p1')
    plt.title(title)
    plt.savefig(graphname)
    plt.show()

def plot_function_with_errors(filename, graphname, title, x, y):

    X = np.loadtxt(filename)[:,0]
    Y = np.loadtxt(filename)[:,1]
    errors = np.loadtxt(filename)[:,2]

    plt.plot(X, Y, linewidth = 0.5)
    plt.errorbar(X, Y, linewidth = 0.5, yerr = errors, ecolor = 'red')
    plt.ylabel(y)
    plt.xlabel(x)
    plt.title(title)
    plt.savefig(graphname)
    plt.show()
    
def main():

    n = int(sys.argv[1])
    nsweeps = int(sys.argv[2])
    p1 = sys.argv[3]
    p2 = sys.argv[4]
    p3 = sys.argv[5]

    animation(n, nsweeps, p1, p2, p3)
    #equil_histo(n, nsweeps)
    #average_infected_sites(n)
    #p4_b(n)
    #p5(n)
    #plot_histogram()

    #plot_contour('SIRS_average_infected_fractions.dat', 'SIRS_average_infected_fractions_contour.png', 'Average Infected Fractions')
    #plot_contour('SIRS_variance.dat', 'SIRS_variance_contour.png', 'Variance')

    #plot_heatmap('SIRS_average_infected_fractions.dat', 'SIRS_average_infected_fractions_heatmap.png', 'Average Infected Fractions')
    #plot_heatmap('SIRS_variance.dat', 'SIRS_variance_heatmap.png', 'Variance')

    #plot_function_with_errors('SIRS_variance_vs_p1.dat', 'SIRS_variance_vs_p1.png', 'Variance vs p1', 'Variance', 'p1')
    #plot_function_with_errors('SIRS_average_inf_vs_immune.dat', 'SIRS_average_inf_vs_immune.png', 'Average Infected Fraction vs Immune Fraction', 'Average Infected Fraction', 'Immune Fraction')

    #average_infected_sites_test(n)

main()