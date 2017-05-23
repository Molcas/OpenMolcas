************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine GenRadQuad_TA(R,nR,nR_Eff,Alpha)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "debug.fh"
      Real*8 R(2,nR-1), Alpha
*
*---- Last point at infinity is eliminated
*
      If (Debug) Then
         Write (6,*) 'Treutler-Ahlrichs Algorithm'
         Write (6,*) 'Alpha=',Alpha
         Write (6,*) 'nR=',nR
      End If
      Fact=Alpha/Log(Two)
      F6=0.6D00
      Do iR = 1, nR-1
         x       = Two*DBLE(iR)/DBLE(nR)-One
         F1=(One+x)
         F2=Log(Two/(One-x))
         R(1,iR) = Fact * F1**F6 * F2
         R(2,iR) = R(1,iR)**2 * Fact * (
     &             F6*F1**(F6-One) * F2
     &           + F1**F6 * One/(One-x)
     &                                 )
     &           * Two / DBLE(nR)
      End Do
      nR_Eff = nR-1
*
      Return
      End
