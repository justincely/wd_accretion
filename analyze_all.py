import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

import glob
import os
import numpy as np



if __name__ == "__main__":
    simulations = glob.glob('simul*')
    all_a = list(set([item.split('_')[1] for item in simulations]))
    all_mass = list(set([item.split('_')[2] for item in simulations]))

    print all_a
    print all_mass

    fig = plt.figure(figsize=(40, 40))

    finished_count = 0

    n = 0
    for a in all_a:
        for m in sorted(all_mass):
            print a, m
            simulations = glob.glob('simul_{}_{}*'.format(a, m))
            
            collisions = []
            for sim in simulations:
                parfile = os.path.join(sim, 'out_parameters.txt')
                if not os.path.exists(parfile): continue

                finished_count += 1
                print "finished: ", finished_count
                with open(parfile, 'r') as txt:
                    for line in txt.readlines():
                        if 'collided' in line:
                            collisions.append(line)
            if not len(collisions):
                continue
            n+=1
            times = np.array([float(item.strip().split()[-2]) for item in collisions])
            print a, m, len(collisions), times
            ax = fig.add_subplot(len(all_a), len(all_mass), n)
            ax.hist(times, bins=np.logspace(0, 8, 30), alpha=.5, label='a={} m={}'.format(a, m.replace('-', '.')))
            ax.grid(True)
            ax.set_xlim(1, 2E8)
            ax.set_ylim(0, 30)
            ax.set_xscale('log')
            ax.set_ylabel('Disruptions')
            
            ax.legend(loc='upper left', shadow=True, numpoints=1)
    ax.set_xlabel('Time (years)')
    fig.savefig('disruptions_all.pdf', bbox_inches='tight')


