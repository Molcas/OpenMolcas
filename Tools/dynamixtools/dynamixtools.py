#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#                                                                      *
# Copyright (C) 2018, Alessio Valentini                                *
#               2018, Luis Manuel Frutos                               *
#               2021,2022, Jonathan Richard Church                     *
#               2021,2022, Igor Schapiro                               *
#***********************************************************************

import numpy as np
import random
import argparse
import os
import sys
import math

def printDict(dictionary):
    '''
    pretty printer for dictionary
    dictionary :: Dictionary
    '''
    for x in dictionary:
        print('{} -> {}'.format(x,dictionary[x]))

def test_initial_things():
    '''
    This function returns a dictionary with test inputs instead of reading them
    dictio :: Dictionary
    These values can be used as a test. I need to set the seed and the results.
    '''
    dictio = {}
    dictio['amu_to_kg'] = 1.66053892173E-27
    dictio['joule_to_kcal'] = 1.439326E+20
    dictio['bohr_to_m'] = 5.2917724900001E-11
    dictio['degrN'] = 3
    dictio['atomN'] = 3
    dictio['T'] = 298.15
    dictio['hbar'] = 1.0545718E-34
    dictio['kb'] = 1.38064852E-23
    dictio['beta'] = 1 / (dictio['T'] * 1.38064852E-23)
    dictio['c'] = 299792458.00
    dictio['cm_to_SI']=dictio['c'] * 100
    dictio['NCMatx'] = np.array([
        [[0.00,0.00,0.07],[0.00,-0.43,-0.56],[0.00,0.43,-0.56]],
        [[0.00,0.00,0.05],[0.00,0.58,-0.40],[0.00,-0.58,-0.40]],
        [[0.00,0.07,0.00],[0.00,-0.55,0.44],[0.00,-0.55,-0.44]]])
    dictio['geom'] = np.array([
             [0.000000, 0.000000, 0.119736],
             [0.000000, 0.761603,-0.478944],
             [0.000000,-0.761603,-0.478944]])
    dictio['RedMass'] = np.array([1.0825,1.0453,1.0810])
    dictio['AtMass'] = np.array([15.9994,1.00794,1.00794])
    dictio['freq'] = np.array([1713.0153, 3727.3731,3849.4717])
    dictio['debug'] = False
    dictio['atomT'] = ['O','H','H']
    return (dictio)

def gaus_dist(sigma):
    '''
    Generate a Gaussian distribution given a sigma
    sigma :: Double
    '''
    u1 = random.uniform(0, 1)
    u2 = random.uniform(0, 1)
    #u1 = 0.5
    #u2 = 0.8
    z = np.sqrt(-2*np.log(u1))*np.sin(2*np.pi*u2)
    return z*sigma

def normal_mode(dictio,label,method):
    '''
    Main driver for initial condition with normal mode sampling. Takes as input a dictionary of inputs.
    dictio :: Dictionary <- inputs
    label :: String <- project name
    in this routine the Boltzmann distribution is calculated for a single initial condition

    degrN :: Int <- number of degrees of freedom
    atomN :: Int <- number of atoms

    some dimensions
    NCMatx :: (degrN, atomN, 3)
    Pcart,vcart :: (atomN,3)

    '''

    # dictionary unpacking
    c = dictio['c']
    beta = dictio['beta']
    freq = dictio['freq']
    geom = dictio['geom']
    atomN = dictio['atomN']
    NCMatx = dictio['NCMatx']
    AtMass = dictio['AtMass']
    RedMass = dictio['RedMass']
    degrN = dictio['degrN']
    debug = dictio['debug']
    atomT = dictio['atomT']
    kb = dictio['kb']
    hbar = dictio['hbar']
    cm_to_SI = dictio['cm_to_SI']
    amu_to_kg = dictio['amu_to_kg']
    joule_to_kcal = dictio['joule_to_kcal']
    bohr_to_m = dictio['bohr_to_m']
    T=dictio['T']
    ##Generate mass matrix from atoms
    mmatrix=mass_matrix(AtMass, atomN)*amu_to_kg
    ##Convert Coordinates to meters for calculations
    xyz_save=geom*bohr_to_m
    xyz_old=np.reshape(xyz_save, (1, 3*atomN))
    ##reshape eigenvector matricices
    NCMatx=np.reshape(NCMatx, (degrN, 3*atomN))
    ##create blank matricies for soon to be created velocities and displacements
    coord_samp=np.zeros((1, 3*atomN), dtype=float)
    ##need to add a test for linearity
    coord_samp_save=np.zeros((degrN, 3*atomN), dtype=float)
    velocity_samp=np.zeros((1, 3*atomN), dtype=float)
    ##This is a test to make sure there were no frequencies included in the hessian calculation from numerical error
    ##Count over normal modes with frequencies greater than 10cm-1 (assumed everything below that is error from the hessian calculation)
    modes=len(freq[freq > 10])
    ##If the number of modes is equal to degrN then proceed assuming that there were no numerical errors, otherwise include all modes but only use those greater than 10cm-1
    if (modes != degrN):
        j=int((3*atomN)-modes)
        COM_modes=j
    else:
        j=0
        COM_modes=0
    ##set variables and begin loop for normal mode analysis or wigner distribution
    Etot=0.0
    i=0
    while (j < degrN):
        ##Convert frequency to SI units
        freqSI=float(freq[j])*cm_to_SI*np.pi*2.0
        ##Normal Mode Sampling for vibrational ground state using a classical boltzmann distribution of energies for each normal mode
        if (method==2):
            ##Generate random number
            rand1=random.uniform(0, 1)
            ##Generate Initial guess of Ei from the cumulative distribution function of normal mode i at temperature T
            Ei=kb*T*(1-rand1)
            Etot=Ei+Etot
            ##Convert Frequencies to joules and calculate Amplitude
            A=np.power((2.0*Ei),0.5)/freqSI
            #Calculate normal coordinate and normal momentum in SI units
            x=A*np.cos(2.0*np.pi*rand1)
            v=-freqSI*A*math.sin(2.0*np.pi*rand1)
        ##Wigner sampling for ground vibrational state
        if (method==3):
            sample=0
            while (sample==0):
                rand1=random.uniform(-1, 1)*np.sqrt(hbar/freqSI)
                rand2=random.uniform(-1, 1)*np.sqrt(hbar*freqSI)
                rand3=random.uniform(0, 1)
                Ei=0.5*(np.power(freqSI*rand1,2)+np.power(rand2,2))
                probability=1.0/(np.pi)*np.exp(-(2.0*Ei)/(hbar*freqSI))
                if (probability > rand3/np.pi):
                    sample=sample+1
                    x=rand1
                    v=rand2
                    Etot=Etot+Ei
        ##Thermal Wigner sampling
        if (method==4):
            sample=0
            while (sample==0):
                alpha=np.tanh(hbar*freqSI/(2.0*kb*T))
                rand1=random.uniform(-1, 1)*np.sqrt(hbar/(freqSI*alpha))
                rand2=random.uniform(-1, 1)*np.sqrt(hbar*freqSI*(1.0/alpha))
                rand3=random.uniform(0, 1)
                Ei=0.5*(np.power(freqSI*rand1,2)+np.power(rand2,2))
                probability=alpha/(np.pi)*np.exp(-2.0*alpha*(np.power(freqSI*rand1,2)+np.power(rand2,2))/(hbar*freqSI))
                if (probability > rand3/(np.pi/alpha)):
                    sample=sample+1
                    x=rand1
                    v=rand2
                    Etot=Etot+Ei
        ##Generate displacements and velocities based on sampling method
        coord_samp=coord_samp+x*NCMatx[j, :]/np.sqrt(mmatrix)
        coord_samp_save[i, :]=x*NCMatx[j, :]/np.sqrt(mmatrix)
        velocity_samp=velocity_samp+v*NCMatx[j, :]/np.sqrt(mmatrix)
        j=j+1
        i=i+1
    ##Add displacement to optimized coordinates and generate cartesian velocities
    vel_new=velocity_samp
    xyz_new=xyz_old+coord_samp
    vel_reshape=np.reshape(vel_new, (atomN, 3))
    xyz_reshape=np.reshape(xyz_new, (atomN, 3))
    ##Check initial total energy based on harmonic oscillator hamiltonian
    E=0.0
    j=COM_modes
    KE=np.power(mmatrix[:]*np.reshape(vel_reshape, (1, 3*atomN)), 2.0)/(2.0*mmatrix[:])
    KE=np.sum(KE)
    while (j < degrN):
        Ei=np.power(float(freq[j])*cm_to_SI*2.0*np.pi, 2.0)*np.power(np.power(mmatrix[:], 0.5)*coord_samp_save[j-COM_modes, :], 2.0)/2.0
        E=E+np.sum(Ei)
        j=j+1
    E=E+KE  
    #This can be uncommented for testing print("initial E, KE, Etot", E-KE, " ", KE, " " , Etot, '\n')
    ##Begin to removal any supurious COM translation or rotation in the molecule and then adjusting the velocities and displacements to conserve energy
    accept=0.0
    Mass=AtMass*amu_to_kg
    while (accept==0):
        ##First calculate the cartesian coordinates with the COM at origin
        com=CenterOfMass(Mass, vel_reshape, atomN)
        com_xyz=CenterOfMass(Mass, xyz_reshape, atomN)
        xyz_reshape_COM=xyz_reshape
        xyz_reshape_COM[:, 0] -= com_xyz[0]
        xyz_reshape_COM[:, 1] -= com_xyz[1]
        xyz_reshape_COM[:, 2] -= com_xyz[2]
        ##Next generate angular momentum, inverse moment of interia matrix, and angular velocity
        Ltot=angular_mo(Mass, xyz_reshape_COM, vel_reshape, atomN)
        invI=inertia(xyz_reshape_COM, Mass, atomN)
        ang_vel=np.dot(invI, Ltot)
        xyz_reshaped=np.reshape(xyz_reshape_COM, (atomN, 3))
        j=0
        while (j < atomN):
            ##Remove spurious rotational velocity
            x=xyz_reshaped[j,0]
            y=xyz_reshaped[j,1]
            z=xyz_reshaped[j,2]
            xyz_corr=np.array([x, y, z])
            ang_corr=angular_vel(Mass, xyz_corr, ang_vel, atomN)
            vel_reshape[j, 0] -= ang_corr[0]
            vel_reshape[j, 1] -= ang_corr[1]
            vel_reshape[j, 2] -= ang_corr[2]
            j=j+1
        ##Remove spurious COM velocity
        com=CenterOfMass(Mass, vel_reshape, atomN)
        vel_reshape[:, 0] -= com[0]
        vel_reshape[:, 1] -= com[1]
        vel_reshape[:, 2] -= com[2]
        ##Calculate harmonic energy and compare to original value
        E=0.0
        j=COM_modes
        KE=np.power(mmatrix[:]*np.reshape(vel_reshape, (1, 3*atomN)), 2.0)/(2.0*mmatrix[:])
        KE=np.sum(KE)
        while (j < degrN):
            Ei=np.power(float(freq[j])*cm_to_SI*2.0*np.pi, 2.0)*np.power(np.power(mmatrix[:], 0.5)*coord_samp_save[j-COM_modes, :], 2.0)/2.0
            E=E+np.sum(Ei)
            j=j+1
        E=E+KE
        #This can be uncommented for testing print("after error removal E, KE, Etot", E-KE, " ", KE, " " , Etot)
        ##If the  value differs by less than 1 percent accept and move to next conndition, otherwise scale the displacements and velocities and try again
        if (abs(E-Etot)/Etot*100 < 1):
            accept=accept+1
            vel_final=vel_reshape
            coord_final=xyz_reshape
        else :
            vel_reshape=vel_reshape*np.power(Etot/E, 0.5)
            xyz_reshape=np.reshape(xyz_old+(coord_samp)*np.power(Etot/E, 0.5), (atomN, 3))
            coord_samp_save=(coord_samp_save)*np.power(Etot/E, 0.5)
    ##Prepare for final coordinate and velocity file writing
    xyz_final=np.reshape(coord_final, (atomN, 3))
    vel_final=np.reshape(vel_final, (atomN, 3))
    ##Convert coordinates to angstroms, velocities to bohr/au
    xyz_final=xyz_final*1.0E10
    vel_final=vel_final*4.57102844e-7
    ##Create final file name
    geomName = label + '.xyz'
    veloName = label + '.velocity.xyz'
    #stringOUT = '\n{}\n{}\n{}\n{}'
    #print(stringOUT.format(veloName, Qcart, geomName, newGeom))
    np.savetxt(veloName,vel_final,fmt='%1.6f')
    # add the atom name in the matrix
    atom_type = np.array(atomT)[:, np.newaxis]
    atom_t_and_geom = np.hstack((atom_type,xyz_final))
    np.savetxt(geomName,atom_t_and_geom,header='{}\n'.format(atomN),comments='',fmt='%s %s %s %s')

def generate_one_boltz(dictio,label):
    '''
    Main driver for initial condition. Takes as input a dictionary of inputs.
    dictio :: Dictionary <- inputs
    label :: String <- project name
    in this routine the Boltzmann distribution is calculated for a single initial condition

    degrN :: Int <- number of degrees of freedom
    atomN :: Int <- number of atoms

    some dimensions
    NCMatx :: (degrN, atomN, 3)
    Pcart,vcart :: (atomN,3)

    '''

    # dictionary unpacking
    c = dictio['c']
    beta = dictio['beta']
    freq = dictio['freq']
    geom = dictio['geom']
    atomN = dictio['atomN']
    NCMatx = dictio['NCMatx']
    AtMass = dictio['AtMass']
    RedMass = dictio['RedMass']
    degrN = dictio['degrN']
    debug = dictio['debug']
    atomT = dictio['atomT']
    kb = dictio['kb']
    hbar = dictio['hbar']
    cm_to_SI = dictio['cm_to_SI']
    amu_to_kg = dictio['amu_to_kg']
    joule_to_kcal = dictio['joule_to_kcal']
    bohr_to_m = dictio['bohr_to_m']
    # DEBUG TIME
    if debug:
        printDict(dictio)
        import pickle
        with open('debug_dictionary_file.pkl', 'wb') as f:
            pickle.dump(dictio, f, pickle.HIGHEST_PROTOCOL)
    # this is lambda :: (degrN)
    lamb = 4 * (np.pi**2) * (freq*100)**2 * c**2
    sigmaP = np.sqrt(RedMass * 1.660537E-27 / beta)
    # change of dimensions
    Pint = (np.array([ gaus_dist(x) for x in sigmaP ]))*6.02214086E21
    # I want to multiply on first axis, both for Pcart and vcart
    # this is why I create and tile the new pint_broadcasted vector
    # it is ugly, there must be a broadcast numpy rule that works better and more efficient
    # but since the degrees of freedom are never gonna be 30.000.000, this loop here is
    # still efficient enough.
    pint_broadcasted = np.empty((degrN,atomN,3))
    for i in range(degrN):
        pint_broadcasted[i] = np.ones((atomN,3)) * Pint[i]
    to_be_summed = NCMatx * pint_broadcasted
    Pcart = np.sum(to_be_summed, axis=0)
    # same multiplication on first axis as with pint
    AtMass_broadcasted = np.empty((atomN,3))
    for i in range(atomN):
        AtMass_broadcasted[i] = np.ones(3) * AtMass[i]
    vcart = Pcart / AtMass_broadcasted
    sigma_q = 1/np.sqrt(beta*lamb)
    # change of dimensions
    Qint = (np.array([ gaus_dist(x) for x in sigma_q ]))*2.4540051E23
    # same multiplication on first axis as with pint and AtMass
    qint_broadcasted = np.empty((degrN,atomN,3))
    for i in range(degrN):
        qint_broadcasted[i] = np.ones((atomN,3)) * Qint[i]
    to_be_summed = NCMatx * qint_broadcasted
    Qcart_temp = np.sum(to_be_summed, axis=0)
    atmass_squared = np.sqrt(AtMass_broadcasted)
    Qcart = Qcart_temp/atmass_squared
    # Pcart is the displacement, added to the main geometry
    # Transformed into ANGSTROM !!
    # in Jan 2019, Dynamix takes geometries in angstrom but velocities in bohr
    newGeom = (Pcart + geom) * 0.529177249
    geomName = label + '.xyz'
    veloName = label + '.velocity.xyz'
    #stringOUT = '\n{}\n{}\n{}\n{}'
    #print(stringOUT.format(veloName, Qcart, geomName, newGeom))
    np.savetxt(veloName,Qcart,fmt='%1.6f')
    # add the atom name in the matrix
    atom_type = np.array(atomT)[:, np.newaxis]
    atom_t_and_geom = np.hstack((atom_type,newGeom))
    np.savetxt(geomName,atom_t_and_geom,header='{}\n'.format(atomN),comments='',fmt='%s %s %s %s')

def parseCL():
    d = '''
This tools is intended to be used as a support to launch molecular dynamics'

usage example:

python3 $MOLCAS/Tools/dynamixtools/dynamixtools.py -t 273 -c 100 -m 1 -i ${Project}.freq.molden

'''
    parser = argparse.ArgumentParser(description=d, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-s", "--seed",
                        dest="seed",
                        required=False,
                        type=int,
                        help="indicate the SEED to use for the generation of randoms")
    parser.add_argument("-l", "--label",
                        dest="label",
                        required=False,
                        type=str,
                        help="label for your project (default is \"geom\")")
    parser.add_argument("-i", "--input",
                        dest="i",
                        required=False,
                        type=str,
                        help="path of the frequency h5 or molden file")
    parser.add_argument("-c", "--condition",
                        dest="condition",
                        required=False,
                        type=int,
                        help="number of initial conditions (default 1)")
    parser.add_argument("-t", "--temperature",
                        dest="temp",
                        required=False,
                        type=float,
                        help="temperature in kelvin for the initial conditions")
    parser.add_argument("-v", "--verbose",
                        dest="debug",
                        required=False,
                        action='store_true',
                        help="more verbose output")
    parser.add_argument("-T", "--TEST",
                        dest="test",
                        required=False,
                        action='store_true',
                        help="keyword use to test the routines")
    parser.add_argument("-D", "--DIGIT",
                        dest="digits",
                        required=False,
                        action='store_true',
                        help="keyword to suppress the counter in the filename (needed for debug)")
    parser.add_argument("-m", "--method",
                        dest="method",
                        required=False,
                        type=int,
                        help=('''\
Keyword to specify the sampling method:
1 Initial conditions based on the molecular vibrational frequencies and energies sampled from a Boltzmann distribution (Default).
2 Thermal normal mode sampling where the cumulitative distribution function for a classical boltzmann distribution at temperature T is used to approximate the energy of each mode.
3 Wigner distribution for the ground vibrational state, n=0.
4 Thermal Wigner distribution for temperature T based on the analytical solution for a canonical ensemble of harmonic oscillators.'''))
    args = parser.parse_args()
    return args

def parseMoldenFreq(fn):
    '''
    This function is the parser for molden freq files.
    '''
    inp = {}
    with open(fn) as f:
        first = f.readline()
        # parse n_freq
        second = f.readline()
        if second.strip() == '[N_FREQ]':
            degrN = int(f.readline().strip())
            inp['degrN'] = degrN
        else:
            sys.exit('This molden format is not recognized N_FREQ')
        # parse frequencies (they should be second)
        freqLabel = f.readline()
        if freqLabel.strip() == '[FREQ]':
            freq = np.empty(degrN)
            for i in range(degrN):
                freq_this_line = f.readline()
                freq[i] = float(freq_this_line.strip())
            inp['freq'] = freq
        else:
            sys.exit('This molden format is not recognized FREQ')
        # parse INT (they should be third)
        intLabel = f.readline()
        if intLabel.strip() == '[INT]':
            intL = np.empty(degrN)
            for i in range(degrN):
                intL_this_line = f.readline()
                intL[i] = float(intL_this_line.strip())
            inp['intL'] = intL
        else:
            sys.exit('This molden format is not recognized INT')
        # parse NATOM
        natomLabel = f.readline()
        if natomLabel.strip() == '[NATOM]':
            atomN = int(f.readline().strip())
            inp['atomN'] = atomN
        else:
            sys.exit('This molden format is not recognized NATOM')
        # parse FR-COORD (they should be fifth field)
        coordLabel = f.readline()
        if coordLabel.strip() == '[FR-COORD]':
            atomT = []
            coordL = np.empty((atomN,3))
            for i in range(atomN):
                geom_this_line = f.readline()
                # H 1.9 1.3 -4.5
                lol = filter(None, geom_this_line.strip('\n').split(' '))
                a,x,y,z = list(lol)
                atomT.append(a)
                coordL[i] = [float(x),float(y),float(z)]
            inp['atomT'] = atomT
            inp['geom'] = coordL
        else:
            sys.exit('This molden format is not recognized FR-COORD')
        # parse FR-NORM-COORD (they should be sixth field)
        norm_coord_Label = f.readline()
        if norm_coord_Label.strip() == '[FR-NORM-COORD]':
            NCMatx = np.empty((degrN,atomN,3))
            for j in range(degrN):
                vibration_line = f.readline()
                for i in range(atomN):
                    geom_this_line = f.readline()
                    evecs = filter(None, geom_this_line.strip('\n').split(' '))
                    x,y,z = list(evecs)
                    NCMatx[j,i] = [float(x),float(y),float(z)]
            inp['NCMatx'] = NCMatx
        else:
            sys.exit('This molden format is not recognized FR-NORM-COORD')
        # parse RMASS (they should be third)
        rmass_Label = f.readline()
        if rmass_Label.strip() == '[RMASS]':
            rmassL = np.empty(degrN)
            for i in range(degrN):
                rmassL_this_line = f.readline()
                rmassL[i] = float(rmassL_this_line.strip())
            inp['RedMass'] = rmassL
        else:
            sys.exit('This molden format is not recognized RMASS')
    # molden file finished, passing the dictionary back
    return(inp)

def parseh5Freq(fn):
    inp = {}
    sys.exit('This feature is still not available')
    return(inp)

def massOf(elem):
    '''
    You get the mass of an element from the label string
    elem :: String
    '''
    dictMass = {'X': 0, 'Ac': 227.028, 'Al': 26.981539, 'Am': 243, 'Sb': 121.757, 'Ar':
            39.948, 'As': 74.92159, 'At': 210, 'Ba': 137.327, 'Bk': 247, 'Be':
            9.012182, 'Bi': 208.98037, 'Bh': 262, 'B': 10.811, 'Br': 79.904,
            'Cd': 112.411, 'Ca': 40.078, 'Cf': 251, 'C': 12.011, 'Ce': 140.115,
            'Cs': 132.90543, 'Cl': 35.4527, 'Cr': 51.9961, 'Co': 58.9332, 'Cu':
            63.546, 'Cm': 247, 'Db': 262, 'Dy': 162.5, 'Es': 252, 'Er': 167.26,
            'Eu': 151.965, 'Fm': 257, 'F': 18.9984032, 'Fr': 223, 'Gd': 157.25,
            'Ga': 69.723, 'Ge': 72.61, 'Au': 196.96654, 'Hf': 178.49, 'Hs':
            265, 'He': 4.002602, 'Ho': 164.93032, 'H': 1.00794, 'In': 114.82,
            'I': 126.90447, 'Ir': 192.22, 'Fe': 55.847, 'Kr': 83.8, 'La':
            138.9055, 'Lr': 262, 'Pb': 207.2, 'Li': 6.941, 'Lu': 174.967, 'Mg':
            24.305, 'Mn': 54.93805, 'Mt': 266, 'Md': 258, 'Hg': 200.59, 'Mo': 95.94,
            'Nd': 144.24, 'Ne': 20.1797, 'Np': 237.048, 'Ni': 58.6934, 'Nb': 92.90638,
            'N': 14.00674, 'No': 259, 'Os': 190.2, 'O': 15.9994, 'Pd': 106.42, 'P':
            30.973762, 'Pt': 195.08, 'Pu': 244, 'Po': 209, 'K': 39.0983, 'Pr':
            140.90765, 'Pm': 145, 'Pa': 231.0359, 'Ra': 226.025, 'Rn': 222,
            'Re': 186.207, 'Rh': 102.9055, 'Rb': 85.4678, 'Ru': 101.07, 'Rf':
            261, 'Sm': 150.36, 'Sc': 44.95591, 'Sg': 263, 'Se': 78.96, 'Si':
            28.0855, 'Ag': 107.8682, 'Na': 22.989768, 'Sr': 87.62, 'S': 32.066,
            'Ta': 180.9479, 'Tc': 98, 'Te': 127.6, 'Tb': 158.92534, 'Tl':
            204.3833, 'Th': 232.0381, 'Tm': 168.93421, 'Sn': 118.71, 'Ti':
            47.88, 'W': 183.85, 'U': 238.0289, 'V': 50.9415, 'Xe': 131.29,
            'Yb': 173.04, 'Y': 88.90585, 'Zn': 65.39, 'Zr': 91.224}
    return(dictMass[elem])

def atomic_masses(atomtype_list):
    '''
    this function takes an atomtype list and return a numpy array with
    atomic masses
    atomtype_list :: [String] <- ['H', 'H', 'O']
    '''
    at_mass = [ massOf(x) for x in atomtype_list ]
    return np.array(at_mass)

def mass_matrix(AtMass, atomN):
    mass=np.zeros((atomN, 1), dtype=float)
    mass_matrix=np.zeros((atomN, 3), dtype=float)
    i=0
    while (i < atomN):
        mass_matrix[i, 0]=AtMass[i]
        mass_matrix[i, 1]=AtMass[i]
        mass_matrix[i, 2]=AtMass[i]
        i=i+1
    mass_matrix=np.reshape(mass_matrix, (1, 3*atomN))
    return mass_matrix

def inertia (xyz,mass,atomN):
    ##I'm not in love with this function and if anyone can do it better please do.
    j=0
    lxx=0
    lyy=0
    lzz=0
    lxy=0
    lxz=0
    lyx=0
    lzx=0
    lyz=0
    lzy=0
    while (j < 3):
        i=0
        while (i < 3):
            h=0
            while (h < atomN):
                stringm=float(mass[h])
                stringx=float(xyz[h, 0])
                stringy=float(xyz[h, 1])
                stringz=float(xyz[h, 2])
                if (i==j) and (i==0) and (j==0):
                    lxx=lxx+stringm*(math.pow(stringy, 2)+math.pow(stringz, 2))
                elif (i==j) and (i==1) and (j==1):
                    lyy=lyy+stringm*(math.pow(stringx, 2)+math.pow(stringz, 2))
                elif (i==j) and (i==2) and (j==2):
                    lzz=lzz+stringm*(math.pow(stringx, 2)+math.pow(stringy, 2))
                elif (i!=j) and (j==0) and (i==1):
                    lxy=lxy+(-stringm*stringx*stringy)
                elif (i!=j) and (j==0) and (i==2):
                    lxz=lxz+(-stringm*stringx*stringz)
                elif (i!=j) and (j==1) and (i==0):
                    lyx=lyx+(-stringm*stringy*stringx)
                elif (i!=j) and (j==1) and (i==2):
                    lyz=lyz+(-stringm*stringy*stringz)
                elif (i!=j) and (j==2) and (i==0):
                    lzx=lzx+(-stringm*stringz*stringx)
                elif (i!=j) and (j==2) and (i==1):
                    lzy=lzy+(-stringm*stringz*stringy)
        
                h=h+1
            i=i+1
        j=j+1
    data1=np.array([[lxx, lxy, lxz], [lyx, lyy, lyz], [lzx, lzy, lzz]])
    invI=np.linalg.inv(data1)
    return invI

def CenterOfMass(mass, vel, atomN):
    h=0
    stringcomx=0
    stringcomy=0
    stringcomz=0
    stringmtot=np.sum(mass[:])
    while (h < atomN):
        stringcomx=stringcomx+float(vel[h, 0])*float(mass[h])
        stringcomy=stringcomy+float(vel[h, 1])*float(mass[h])
        stringcomz=stringcomz+float(vel[h, 2])*float(mass[h])
        h=h+1
    comx=stringcomx/stringmtot
    comy=stringcomy/stringmtot
    comz=stringcomz/stringmtot
    com=np.array([comx, comy, comz])
    return com

def angular_mo (mass, xyz, vel, atomN):
    Lxtot=0
    Lytot=0
    Lztot=0
    h=0
    xyz=np.reshape(xyz, (atomN, 3))
    while (h < atomN):
        px=float(vel[h, 0])*float(mass[h])
        py=float(vel[h, 1])*float(mass[h])
        pz=float(vel[h, 2])*float(mass[h])
        Lx=float(xyz[h, 1])*pz-float(xyz[h, 2])*py
        Ly=float(xyz[h, 2])*px-float(xyz[h, 0])*pz
        Lz=float(xyz[h, 0])*py-float(xyz[h, 1])*px
        Lxtot=Lxtot+Lx
        Lytot=Lytot+Ly
        Lztot=Lztot+Lz
        h=h+1
    Ltot=np.array([Lxtot, Lytot, Lztot])
    return Ltot

def angular_vel (mass, xyz, vel, atomN):
    px=float(vel[0])
    py=float(vel[1])
    pz=float(vel[2])
    wx=-float(xyz[1])*pz+float(xyz[2])*py
    wy=-float(xyz[2])*px+float(xyz[0])*pz
    wz=-float(xyz[0])*py+float(xyz[1])*px
    wtot=np.array([wx, wy, wz])
    return wtot

def main():
    print('')
    args = parseCL()

    if args.test:
        inputs = test_initial_things()
        seedI = 123456789
        random.seed(seedI)
        complete_label = 'test_bolz'
        generate_one_boltz(inputs,complete_label)
        print(inputs)
    else:
        # i is input file from command line
        if args.i:
            fn=args.i
        else:
            # I do not like this termination here, but I still have to figure out how
            # to properly do mutually exclusive argparse keywords.
            # I will keep this exit code here in the meanwhile...
            sys.exit('-i input freq file is a required keyword, --help for help')
        if args.seed:
            seedI = args.seed
            print('seed set to: {}'.format(seedI))
            random.seed(seedI)
        if args.label:
            label = args.label
            print('Project label set to: {}'.format(label))
        else:
            label = 'geom'
        if args.condition:
            number_of_ic = args.condition
        else:
            number_of_ic = 1
        if args.method:
            method = args.method
        else:
            method = 1
        # check if is is molden or h5
        name, ext = os.path.splitext(fn)
        if ext == '.molden':
            inputs = parseMoldenFreq(fn)
        elif ext == '.h5':
            inputs = parseh5Freq(fn)
        else:
            print('You must use freq.molden or .h5 files')
        if args.debug:
            inputs['debug'] = True
        else:
            inputs['debug'] = False
        if args.temp == None:
            print('\nSetting default temperature to 300 K. Use the -t keyword to set a custom temperature.\n')
            inputs['T'] = 300
        else:
            inputs['T'] = args.temp
        ##Conversion Factors and constants
        inputs['kb'] = 1.38064852E-23
        inputs['h'] = 6.62607004E-34
        inputs['hbar'] = 1.0545718E-34
        inputs['beta'] = 1 / (inputs['T'] * inputs['kb'])
        inputs['c'] = 299792458.00
        inputs['cm_to_SI'] = inputs['c'] * 100
        inputs['amu_to_kg'] = 1.66053892173E-27
        inputs['joule_to_kcal'] = 1.439326E+20
        inputs['bohr_to_m'] = 5.2917724900001E-11
        inputs['AtMass'] = atomic_masses(inputs['atomT'])
        for counter in range(number_of_ic):
            if args.digits:
                complete_label = '{}'.format(label)
            else:
                complete_label = '{}{:04}'.format(label,counter)
            if (method==1):
                    generate_one_boltz(inputs,complete_label)
            elif (method==2 or method==3 or method==4):
                    normal_mode(inputs,complete_label,method)
        print('\nThis routine generates geometries in angstrom and velocities in bohr (the format that Molcas requires for a Semiclassical Molecular Dynamics\n')

if __name__ == "__main__":
    main()
