"""Create scenarios for running the MERCURY code

Evenly sample the 2:1 resonance for the Sun - Jupiter system under different
solar system distributions.
    1. Sampling Jupter's semi-major axis at 1x, 2x, and 3x
    2. Sampling Jupiter's mass from .1x to 10x
    3. Using 10,000 - 20,000 small bodies

Fully sampling the libration width will allow later analysis to simply convolve
different asteroid distributions with this data.

Run for 250Myr

asteroid inclinations from -30 to 30 degrees

Metals are in the UV or optical, dust in the IR.  Metal polluted are DAZ, DBZ, DZ

"""

import os
import numpy as np
import shutil
import random
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

#-- solar masses
JUPITER_MASS = 9.54791938424326609E-04
#-- in fortran already
J_HILL_RADIUS = 3
J_DENSITY = 1.33

#-- AU
A_JUP = 5.20336301
E_JUP = .0489
I_JUP = 1.30530
G_JUP = 14.75385
N_JUP = 100.55615
M_JUP = 0

#---
N_SMALL = 100
E_MAX = .6

N_SIMUL = 200

#-------------------------------------------------------------------------------

def write_input_file(outname, bodies, coords, filetype='small'):

    if isinstance(bodies, str):
        bodies = [bodies]

    if isinstance(coords, str):
        coords = [coords]

    with open(outname, 'w') as outtxt:
        if filetype == 'small':
            outtxt.write(")O+_06 Small-body initial data  (WARNING: Do not delete this line!!)\n")
            outtxt.write(") Lines beginning with <)> are ignored.\n")
            outtxt.write(")---------------------------------------------------------------------\n")
            outtxt.write("style (Cartesian, Asteroidal, Cometary) = Asteroidal\n")
            outtxt.write(")---------------------------------------------------------------------\n")

        if filetype == 'big':
            outtxt.write(")O+_06 Big-body initial data  (WARNING: Do not delete this line!!)\n")
            outtxt.write(") Lines beginning with `)' are ignored.\n")
            outtxt.write(")---------------------------------------------------------------------\n")
            outtxt.write("style (Cartesian, Asteroidal, Cometary) = Asteroidal\n")
            outtxt.write("epoch (in days) = 0.0000000\n")
            outtxt.write(")---------------------------------------------------------------------\n")

        for body, coord in zip(bodies, coords):
            outtxt.write(body + '\n')
            outtxt.write(coord + '\n')

#-------------------------------------------------------------------------------

def libration_width(e, a, m_sun, m_planet):
    """

    Parameters:
    -----------
    e : float, int
        eccentricity
    a : float, int
        semi-major axis
    m_sun : float
        mass of the central body
    m_planet : flaot
        mass of the large planet

    Returns:
    --------
    max_l : float
        +/- maximum libration width

    """

    #-- for 2:1 resonance
    j2 = -1
    jacobi = abs((m_planet / m_sun) * (-.749964))

    max_l = ((16.0 / 3.0) * jacobi * e) ** (1/2.)
    max_l *= (1 + (jacobi/(27.0 * j2**2 *(e**3))) )**(1/2.)
    max_l *= a
    max_l -= (2.0/(9.0 * j2 * e)) * jacobi * a

    return max_l

#-------------------------------------------------------------------------------

def rel_period(a1, a2):
    """

    Parameters:
    -----------
    a1 : float
        semi-major axis of first body
    a2 : float
        semi-major axis of second body

    Returns:
    --------
    T : float
        relative orbital period

    """

    return np.sqrt(a1**3 / a2**3)

#-------------------------------------------------------------------------------

def make_condor_submit(root, n_simul):
    """Write out the condor submit file"""

    with open(os.path.split(root)[-1]+'mercury.submit', 'w') as outtxt:
        outtxt.write('################\n')
        outtxt.write('# Submit the simulations with condor\n')
        outtxt.write('################\n')
        outtxt.write('\n')
        outtxt.write('Executable  = call_mercury\n')
        outtxt.write('Universe    = vanilla\n')
        outtxt.write('priority    = 10\n')
        outtxt.write('getenv      = true\n')
        outtxt.write('\n')
        outtxt.write('notification = Complete\n')
        outtxt.write('notify_user  = ely@stsci.edu \n')
        outtxt.write('RUN_PREFIX   = {}\n'.format(root))
        outtxt.write('\n')
        outtxt.write('transfer_executable = true\n')
        outtxt.write('\n')
        outtxt.write('error      = $(RUN_PREFIX)$(Process)/condor.err\n')
        outtxt.write('output     = $(RUN_PREFIX)$(Process)/condor.out\n')
        outtxt.write('log        = $(RUN_PREFIX)$(Process)/condor.log\n')
        outtxt.write('arguments  = $(RUN_PREFIX)$(Process)\n')
        outtxt.write('\n')
        outtxt.write('queue {}\n'.format(n_simul))        

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    for a_factor in np.arange(1, 8, .5):
        print a_factor
        masses_sample = [.1, .5, 1, 2, 4, 6, 8, 10]
        for mass_factor in masses_sample:
            print mass_factor
    
            rootname = 'data/simul_{}_{}_'.format(a_factor, mass_factor)
            make_condor_submit(rootname, N_SIMUL)

            for N in xrange(N_SIMUL):
                #-- Big body
                N = str(N)

                simul_name = rootname + N

                if not os.path.exists(simul_name):
                    os.makedirs(simul_name)

                big_outname = 'big_{}_{}.in'.format(a_factor, str(mass_factor).replace('.', '-'))
                print "Generating big file: {}".format(big_outname)

                body = 'JUPITER r={:1.5E} d={:1.5E} m={:1.15E} '.format(J_HILL_RADIUS,
                                                                        J_DENSITY,
                                                                        JUPITER_MASS * mass_factor)
                #-- a, e, i, g, n, m
                coords = '    {:1.15E} {:1.15E} {:1.15E} {:1.15E} {:1.15E} {:1.15E} 0 0 0 '.format(A_JUP * a_factor,
                                                                                                   E_JUP,
                                                                                                   I_JUP,
                                                                                                   G_JUP,
                                                                                                   N_JUP,
                                                                                                   M_JUP)
                write_input_file(big_outname, body, coords, 'big')

                #-- Small bodies

                a_res = ((1/2.)**2 * (A_JUP * a_factor)**3) ** (1/3.)
                all_e = np.random.random_sample(N_SMALL) * E_MAX
                #-- Msun is MS mass
                lib_width = libration_width(E_MAX, a_res, .54, JUPITER_MASS * mass_factor)

                max_a = a_res + (2 * lib_width)
                min_a = a_res - (2 * lib_width)
                print a_res, lib_width, min_a, max_a
                all_a = np.array([random.uniform(min_a, max_a) for i in xrange(N_SMALL)])

                all_i = np.array([random.uniform(-30, 30) for i in xrange(N_SMALL)])
                all_g = np.array([random.uniform(0, 360) for i in xrange(N_SMALL)])
                all_n = np.array([random.uniform(0, 360) for i in xrange(N_SMALL)])
                all_m = np.array([random.uniform(0, 360) for i in xrange(N_SMALL)])

                bodies = []
                coords = []
                for num, a, e, i, g, h, m in zip(xrange(N_SMALL),
                                                all_a,
                                                all_e,
                                                all_i,
                                                all_g,
                                                all_n,
                                                all_m):
                    bodies.append("{}   ep=0".format(num))
                    coords.append("    {}    {}    {}    {}    {}    {}    0  0  0".format(a, e, i, g, h, m))

                small_outname = 'small_{}_{}.in'.format(a_factor, str(mass_factor).replace('.', '-'))
                write_input_file(small_outname, bodies, coords, 'small')

                #-- Plot the distribution
                fig = plt.figure()
                ax = fig.add_subplot(1, 1, 1)
                ax.plot(all_a, all_e, '.', alpha=.6, color='black')
                ax.axvline(x=a_res, color='r', ls='-', lw=4, label="Resonance")

                lib_e = np.linspace(0, 1, 200)
                lib_a = np.linspace(all_a.min(), all_a.max(), 200)
                ax.plot(a_res - libration_width(lib_e, lib_a, .54, JUPITER_MASS * mass_factor), lib_e, color='b', ls='-', lw=4, label='Libration_width')
                ax.plot(a_res + libration_width(lib_e, lib_a, .54, JUPITER_MASS * mass_factor), lib_e,  color='b', ls='-', lw=4)

                #ax.axvline(x=a_res - lib_width, color='r', ls='--', lw=4, label="Lib width at e=.6")
                #ax.axvline(x=a_res + lib_width, color='r', ls='--', lw=4)
                ax.set_ylim(0, .6)
                ax.set_xlim(all_a.min(), all_a.max())

                ax.set_xlabel('Semi-major Axis: a (AU)')
                ax.set_ylabel('Eccentricity: e')
                ax.set_title('Asteroid Distribution: a={}, m={}'.format(a_factor,
                                                          str(mass_factor).replace('.', '-')))
                ax.legend(shadow=True, numpoints=1, loc='upper right')
                fig.savefig('small_dist_{}_{}.pdf'.format(a_factor,
                                                          str(mass_factor).replace('.', '-')), bbox_inches='tight')
                plt.close(fig)
                #-----------------------


                with open('files.in', 'w') as out_txt:
                    out_txt.write(' {}\n'.format(big_outname))
                    out_txt.write(' {}\n'.format(small_outname))
                    out_txt.write(' in_parameters.txt\n')
                    out_txt.write(' out_positions.txt\n')
                    out_txt.write(' out_encounters.txt\n')
                    out_txt.write(' out_parameters.txt\n')
                    out_txt.write(' dump_big.txt\n')
                    out_txt.write(' dump_small.txt\n')
                    out_txt.write(' dump_parameters.txt\n')
                    out_txt.write(' dump_restart.txt\n')


                needed_files = ['mercury_code/Makefile',
                                'mercury_code/mercury6_3.for',
                                'mercury_code/mercury.inc',
                                'mercury_code/close6.for',
                                'mercury_code/close.in',
                                'mercury_code/element6.for',
                                'mercury_code/element.in',
                                'mercury_code/in_parameters.txt',
                                'mercury_code/element.in',
                                'mercury_code/message.in',
                                'mercury_code/swift.inc',
                                'files.in',
   				                small_outname,
           			            big_outname]

                for filename in needed_files:
                    shutil.copy(filename, simul_name)

                os.remove(small_outname)
                os.remove(big_outname)
                os.remove('files.in')
                os.system('rm -f small_dist*.pdf')
