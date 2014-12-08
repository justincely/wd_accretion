"""Create scenarios for running the MERCURY code

Evenly sample the 2:1 resonance for the Sun - Jupiter system under different
solar system distributions.
    1. Sampling Jupter's semi-major axis at 1x, 2x, and 3x
    2. Sampling Jupiter's mass from .1x to 10x
    3. Using 10,000 - 20,000 small bodies

Fully sampling the libration width will allow later analysis to simply convolve
different asteroid distributions with this data.

Run for 250Myr

"""

import numpy as np
import random
import matplotlib.pyplot as plt
plt.ioff()

#-- solar masses
JUPITER_MASS = 9.54791938424326609E-04
#-- in fortran already
J_HILL_RADIUS = r='3.d0'
J_DENSITY = 1.33

#-- AU
A_JUP = 5.20336301
E_JUP = .0489
I_JUP = 1.30530
G_JUP = 14.75385
N_JUP = 100.55615
M_JUP = 0

#---
N_SMALL = 10000
E_MAX = .6

#-------------------------------------------------------------------------------

def generate_small(n_bodies, mass_j, radius_j, n_files):
    """Create the list of small bodies



    """

    pass

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
            outtxt.write("style (Cartesian, Asteroidal, Cometary) = Aseroidal\n")
            outtxt.write("epoch (in days) = 0\n")
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
    max_l -= (2/(9 * j2 * e)) * jacobi * a

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

if __name__ == "__main__":
    for a_factor in [1, 2, 3]:
        for mass_factor in [.1, .5, 1, 5, 10]:
            #-- Big body
            big_outname = 'big_{}_{}.in'.format(a_factor, str(mass_factor).replace('.', '-'))
            print "Generating big file: {}".format(big_outname)

            body = 'JUPITER m={} r={} d={}'.format(JUPITER_MASS * mass_factor,
                                                   J_HILL_RADIUS,
                                                   J_DENSITY)
            #-- a, e, i, g, n, m
            coords = '    {} {} {} {} {} {}'.format(A_JUP * a_factor,
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
            lib_width = libration_width(E_MAX, a_res, 1, JUPITER_MASS * mass_factor)
            print a_res, lib_width

            max_a = a_res + 2 * lib_width
            min_a = a_res - 2 * lib_width
            all_a = np.array([random.uniform(min_a, max_a) for i in xrange(N_SMALL)])



            #-- Plot the distribution
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.plot(all_a, all_e, '.', alpha=.6, color='black')
            ax.axvline(x=a_res, color='r', ls='-', lw=4, label="Resonance")
            ax.axvline(x=a_res - lib_width, color='r', ls='--', lw=4, label="Lib width at e=.6")
            ax.axvline(x=a_res + lib_width, color='r', ls='--', lw=4)
            ax.set_ylim(0, .6)
            ax.set_xlim(all_a.min(), all_a.max())

            ax.set_xlabel('a (AU)')
            ax.set_ylabel('e')
            ax.set_title('Asteroid Distribution around resonance')
            fig.savefig('small_dist_{}_{}.pdf'.format(a_factor,
                                                      str(mass_factor).replace('.', '-')))
            plt.close(fig)
            #--
