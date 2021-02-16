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
      Subroutine PotDist(F,V,Lambda,PED,NumInt,nOsc)
C!
C!  Purpose:
C!    To give the fractional contributions of the F-matrix to the
C!    potential energy.
C!
C!  Input:
C!    F        : Real*8 two dimensional array -  contains
C!               the force constants expressed in internal
C!    V        : Real*8 two dimensional array  - contains
C!               the eigenvectors of F*G as columns.
C!    Lambda   : Real*8 array - contains the eigenvalues
C!               of F*G.
C!
C!  Output:
C!    PED      : Real*8 three dimensional array - Potential
C!               Energy Distribution for each mode.
C!  Uses:
C!    Linalg
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1996.
C!
c       Use LinAlg
      Implicit Real*8 ( a-h,o-z )
#include "Constants_mula.fh"
      Real*8 F( NumInt,NumInt )
      Real*8 V ( NumInt,NumInt )
      Real*8  Denominator
      Real*8 Lambda( NumInt )
      Real*8 PED ( nOsc,nOsc,nOsc )
C!
C!---- Initialize.
C!
C!---- Calculate Potential Energy Distribution for each mode.
c       PED = 0.0d0
      call dcopy_(nOsc*nOsc*nOsc,0.0d0,0,PED,1)
      Do i = 1,NumInt
      Denominator=max(1.0D-10,Lambda(i))
      Do k = 1,NumInt
      Do l = 1,NumInt
      If ( k.eq.l ) Then
      PED(k,k,i) = (V(k,i)**2*F(k,k))/Denominator
      Else
      PED(k,l,i) = (2*V(k,i)*V(l,i)*F(k,l))/Denominator
      End If
      End Do
      End Do
      End Do
C!
      End
