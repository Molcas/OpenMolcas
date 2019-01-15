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
#***********************************************************************

import numpy as np
import random
import argparse
import os
import sys


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
    dictio['degrN'] = 3
    dictio['atomN'] = 3
    dictio['T'] = 298.15
    dictio['kb'] = 1.38064852E-23
    dictio['beta'] = 1 / (dictio['T'] * 1.38064852E-23)
    dictio['c'] = 299792458.00
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
    Generate a gaussian distribution given a sigma
    sigma :: Double
    '''
    u1 = random.uniform(0, 1)
    u2 = random.uniform(0, 1)
    #u1 = 0.5
    #u2 = 0.8
    z = np.sqrt(-2*np.log(u1))*np.sin(2*np.pi*u2)
    return z*sigma


def generate_one_boltz(dictio,label):
    '''
    Main driver for initial condition. Takes as input a dictionary of inputs.
    dictio :: Dictionary <- inputs
    label :: String <- project name
    in this routine the boltzmann distribution is calculated for a single initial condition

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

python3 $MOLCAS/Tools/dynamixtools/dynamixtools.py -t 273 -b 100 -i ${Project}.freq.molden

'''
    parser = argparse.ArgumentParser(description=d)
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
    parser.add_argument("-b", "--boltzmann",
                        dest="bol",
                        required=False,
                        type=int,
                        help="number of initial condition following boltzmann distribution (default 1)")
    parser.add_argument("-t", "--temperature",
                        dest="temp",
                        required=False,
                        type=float,
                        help="temperature in Kelvin for the initial conditions")
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
                    # H 1.9 1.3 -4.5
                    lol = filter(None, geom_this_line.strip('\n').split(' '))
                    x,y,z = list(lol)
                    NCMatx[j,i] = [float(x),float(y),float(z)]
            inp['NCMatx'] = NCMatx
        else:
            sys.exit('This molden format is not recognized FR-COORD')

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

def main():
    print('')
    args = parseCL()
    print(args)

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

            sys.exit('-i input freq file is a required keyword')

        if args.seed:
            seedI = args.seed
            print('seed set to: {}'.format(seedI))
            random.seed(seedI)

        if args.label:
            label = args.label
            print('Project label set to: {}'.format(label))
        else:
            label = 'geom'

        if args.bol:
            # right now this is only boltzmann, we need to rethink this IF when wigner is implemented
            number_of_ic = args.bol
        else:
            number_of_ic = 1

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

        inputs['kb'] = 1.38064852E-23
        inputs['beta'] = 1 / (inputs['T'] * 1.38064852E-23)
        inputs['c'] = 299792458.00
        inputs['AtMass'] = atomic_masses(inputs['atomT'])



        #print('\n\n\n\nAFTER DICTIONARY')
        for counter in range(number_of_ic):
            complete_label = '{}{:04}'.format(label,counter)
            generate_one_boltz(inputs,complete_label)
        print('\nThis routine generates geometries in Angstrom and velocities in Bohr (the format that Molcas requires for a Semiclassical Molecular Dynamics)\n')



if __name__ == "__main__":
    main()

