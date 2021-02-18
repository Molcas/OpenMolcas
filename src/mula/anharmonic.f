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
      Subroutine AnharmonicFreq(x_anharm,harmfreq,anharmfreq,nOsc)
C!
C!  Purpose:
C!    Calculate the anharmonic frequencies.
C!
C!  Input:
C!    x_anharm   : Real*8 two dimensional array - anharmonicity
C!                 constants.
C!    harmfreq   : Real*8 array - harmonic frequencies.
C!
C!  Output:
C!    anharmfreq : Real*8 array - anharmonic frequencies.
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1996.
C!
      Implicit Real*8 ( a-h,o-z )
#include "Constants_mula.fh"
      Real*8 G1,G2
      Real*8 x_anharm (nosc,nosc)
      Real*8 harmfreq (nosc)
      Real*8 anharmfreq (nosc)
#include "WrkSpc.fh"
C!
c       nDim = nOsc
      Call GetMem('level1','Allo','Inte',iplevel1,nOsc)
      Call GetMem('level2','Allo','Inte',iplevel2,nOsc)
C!
      G1 = 0.0d0
      G2 = G1
      do iv=0,nOsc-1
      iWork(iplevel1+iv)=0
      enddo
c       level1 = 0
      Do istate = 1,nOsc
      do iv=0,nOsc-1
      iWork(iplevel2+iv)=0
      enddo
c          level2 = 0
      iWork(iplevel2+istate-1) = 1
c          l_harm=nOsc
      Call TransEnergy(G1,x_anharm,harmfreq,iWork(iplevel1),
     &                     G2,x_anharm,harmfreq,iWork(iplevel2),
     &    anharmfreq(istate),nOsc)
      End Do
C!
      Call GetMem('level1','Free','Inte',iplevel1,nOsc)
      Call GetMem('level2','Free','Inte',iplevel2,nOsc)
C!
      End
