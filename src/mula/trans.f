************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1996, Niclas Forsberg                                  *
************************************************************************
C!-----------------------------------------------------------------------!
C!
      Subroutine TransEnergy(G01,x_anharm1,harmfreq1,level1,
     &                             G02,x_anharm2,harmfreq2,
     &                              level2,energy,nDim)
C!
C!  Purpose:
C!    Calculate transition energy between two states.
C!
C!  Input:
C!    G01,G02    : Real*8 - (G02-G01) = energy difference between
C!                 the two states.
C!    x_anharm1  : Real*8 two dimensional array - anharmonicity
C!                 constants for ground state.
C!    x_anharm2  : Real*8 two dimensional array - anharmonicity
C!                 constants for excited state.
C!    harmfreq1  : Real*8 array - harmonic frequencies for
C!                 ground state.
C!    harmfreq2  : Real*8 array - harmonic frequencies for
C!                 excited state.
C!    level1     : Integer array - quanta for ground state.
C!    level2     : Integer array - quanta for excited state.
C!
C!  Output:
C!    energy     : Real*8 variable - energy for transition
C!                 between level1 and level2.
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1996.
C!
      Implicit Real*8 ( a-h,o-z )
#include "Constants_mula.fh"
      Real*8    sum,G0,G1,G01,G02
      Real*8    energy

      Real*8 x_anharm1(nDim,nDim),x_anharm2(nDim,nDim)
      Real*8 harmfreq1(nDim),harmfreq2(nDim)
      Integer level1(nDim),level2(nDim)
C!
C!
C!---- Calculate energy for level 1.
      sum = G01
      Do i = 1,nDim
      sum = sum+harmfreq1(i)*(level1(i)+0.5d0)
      Do j = i,nDim
      sum = sum+x_anharm1(i,j)*(level1(i)+0.5d0)*(level1(j)+0.5d0)
      End Do
      End Do
      G0 = sum
C!
C!---- Calculate energy for level 2.
      sum = G02
      Do i = 1,nDim
      sum = sum+harmfreq2(i)*(level2(i)+0.5d0)
      Do j = i,nDim
      sum = sum+x_anharm2(i,j)*(level2(i)+0.5d0)*(level2(j)+0.5d0)
      End Do
      End Do
      G1 = sum
C!
C!---- Calculate energy difference.
      energy = G1-G0
C!
      End
C!
