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
# Copyright (C) 2018, Alessio Valentini, Luis Manuel Frutos            *
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

def giveMeDataJefe():
    '''
    This function returns a dictionary with test inputs instead of reading them
    dictio :: Dictionary
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

    some dimesions
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

    # this is lambda :: (degrN)
    lamb = 4 * (np.pi**2) * (freq*100)**2 * c**2
    sigmaP = np.sqrt(RedMass * 1.660537E-27 / beta)
    # change of dimensions
    Pint = (np.array([ gaus_dist(x) for x in sigmaP ]))*6.02214086E21

    # I want to multiply on first axis, both for Pcart and vcart
    # this is why I am broadcasting/transposing
    pint_transpose = np.broadcast_to(Pint,(degrN,atomN,3)).T
    to_be_summed = NCMatx * pint_transpose
    Pcart = np.sum(to_be_summed, axis=0)
    AtMass_transpose = np.broadcast_to(AtMass,(atomN,3)).T
    vcart = Pcart / AtMass_transpose

    sigma_q = 1/np.sqrt(beta*lamb)
    # change of dimensions
    Qint = (np.array([ gaus_dist(x) for x in sigma_q ]))*2.4540051E23

    # same multiplication on first axis.
    qint_transpose = np.broadcast_to(Qint,(degrN,atomN,3)).T
    to_be_summed = NCMatx * qint_transpose
    Qcart_temp = np.sum(to_be_summed, axis=0)
    atmass_squared = np.sqrt(AtMass_transpose)
    Qcart = Qcart_temp/atmass_squared

    newGeom = Pcart + geom
    geomName = label + '.xyz'
    veloName = label + '.velocity.xyz'
    stringOUT = '\n{}\n{}\n{}\n{}'
    #print(stringOUT.format(veloName, Qcart, geomName, newGeom))
    np.savetxt(veloName,Qcart)
    np.savetxt(geomName,newGeom,header='{}\n'.format(atomN),comments='')

def parseCL():
    parser = argparse.ArgumentParser()
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
    parser.add_argument("-n", "--number",
                        dest="num",
                        required=False,
                        type=int,
                        help="how many initial condition needed? (default 1)")
    parser.add_argument("-t", "--temperature",
                        dest="temp",
                        required=True,
                        type=float,
                        help="temperature in Kelvin for the initial conditions")
    args = parser.parse_args()
    return args

def parseMoldenFreq(fn):
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

def atomic_masses(atomtype_list):
    '''
    this function takes an atomtype list and return a numpy array with
    atomic masses
    atomtype_list :: [String] <- ['H', 'H', 'O']
    '''
    atomMass = {'H': 1.00794, 'O': 15.9994}
    at_mass = [ atomMass[x] for x in atomtype_list ]
    return np.array(at_mass)

def main():
    args = parseCL()
    # i is input file from command line
    if args.i:
        fn=args.i
    else:
        fn=''

    if args.seed:
        seedI = args.seed
        print('seed set to: {}'.format(seedI))
        random.seed(seedI)

    if args.label:
        label = args.label
        print('Project label set to: {}'.format(label))
    else:
        label = 'geom'

    if args.num:
        number_of_ic = args.num
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

    ##################################################
    # take out this line to make it work as intended #
    # inputs = giveMeDataJefe()
    ##################################################

    inputs['T'] = args.temp
    inputs['kb'] = 1.38064852E-23
    inputs['beta'] = 1 / (inputs['T'] * 1.38064852E-23)
    inputs['c'] = 299792458.00
    inputs['AtMass'] = atomic_masses(inputs['atomT'])

    printDict(inputs)

    #print('\n\n\n\nAFTER DICTIONARY')
    for counter in range(number_of_ic):
        complete_label = '{}{:04}'.format(label,counter)
        generate_one_boltz(inputs,complete_label)



if __name__ == "__main__":
    main()

