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

    # some dimesions
    # NCMatx :: (degrN, atomN, 3) 
    # Pcart,vcart :: (atomN,3)

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
    args = parser.parse_args()
    return args

def parseMoldenFreq(fn):
    inp = {}
    return(inp)

def parseh5Freq(fn):
    inp = {}
    return(inp)

def main():
    args = parseCL()
    # i is input file from command line
    if args.i:
        fn=args.i
    else:
        fn='/home/alessio/Molcas-Initial/Test/water.freq.molden'

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
        inputs = giveMeDataJefe()

    ##################################################
    # take out this line to make it work as intended #
    inputs = giveMeDataJefe()
    ##################################################

    printDict(inputs)

    print('\n\n\n\nAFTER DICTIONARY')
    for counter in range(number_of_ic):
        complete_label = '{}{:04}'.format(label,counter)
        generate_one_boltz(inputs,complete_label)



if __name__ == "__main__":
    main()

#      
#      
#      write(2,*) Qcart(2,2)
#      go to 1
#      
#      write(1,*) 'Cartesian geometry:'
#       write(1,*) 'O  ',q0cart(1,1)+Qcart(1,1),q0cart(2,1)+Qcart(2,1),
#     &               q0cart(3,1)+Qcart(3,1)
#       write(1,*) 'H  ',q0cart(1,2)+Qcart(1,2),q0cart(2,3)+Qcart(2,3),
#     &               q0cart(3,2)+Qcart(3,2)
#       write(1,*) 'H  ',q0cart(1,3)+Qcart(1,3),q0cart(2,3)+Qcart(2,3),
#     &               q0cart(3,3)+Qcart(3,3)
#
#       write(1,*) 'Cartesian momentum:'
#       write(1,*) 'O  ',Pcart(1,1),Pcart(2,1),Pcart(3,1)
#       write(1,*) 'H  ',Pcart(1,2),Pcart(2,2),Pcart(3,2)
#       write(1,*) 'H  ',Pcart(1,3),Pcart(2,3),Pcart(3,3)
#      
#      END

#      subroutine GaussDistr (sigma,z)
#      real sigma,z,u1,u2
#C      integer,parameter :: seed = 86456
#C generate an "x" number randomly distributed among (0,1)
#C      x=0.1
#      u1=RAND(0)
#      u2=RAND(0)
#      z=sqrt(-2*log(u1))*sin(2*3.141592654*u2)
#      z=z*sigma
#      return
#      end subroutine

