#!/usr/bin/env python3

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
# Copyright (C) 2023, Jochen Autschbach                                *
#***********************************************************************

'''
A simple convertor that converts hyperfine coupling tensor constants from atomic units (a.u.) to mHZ.
                1 a.u. = (2*)con*gnu mHZ
		Note the factor of 2. 
                For UHF calculations in HF
                the factor of 2 is in needed for the conversion,
                but when it comes to rassi results in OPENMOLCAS, the factor of 
                2 is discarded, since the formula contains that factor of 2
                when calculating HFC.
Input:
        
Output:
'''
import sys

con = 95.5213160024483 #conversion factor not including nuclear g-factor.


# nuclear g-factors (gnu) for specific isotopes. Any isotope not included can be added to the list.
gnu = {"Np": 1.2560,  
       "H" : 5.58569,
       "Cs": 0.737978,
       "Hg": 1.011770,
       "Fr": 0.780,
       "C" : 1.4048,
       "O" : -0.7575,
       "F" : 5.25773,
       "Si": -1.11057,
       "U" : -0.1228,
       "Ti": -0.315392,
       "S" : 0.42921,
       "Ce": 0.314286,
       "Cl": 0.5479157,
       "Er": -0.1618571,
       "159Tb": 1.343,
       "241Am": 0.636
       }

print("""Usage: 1. python au2mhz.py [Isotope] [HFC Value in a.u.]
       2. python au2mhz.py RASSI [Isotope] [HFC Value in a.u.]
       3. python au2mhz.py HF [Isotope] [HFC Value in a.u.]
       versions 1 and 2 are equivalent for RASSI
       check or amend the gnu data list for isotopes 
       """)
print("The converted value in MHz:")
if len(sys.argv) > 1:
   if sys.argv[1] == 'HF': 
        print(float(sys.argv[3])*2.0*con*gnu[sys.argv[2]])
   elif sys.argv[1] == 'RASSI':
        print(float(sys.argv[3])*con*gnu[sys.argv[2]])
   else:
        print(float(sys.argv[2])*con*gnu[sys.argv[1]])

