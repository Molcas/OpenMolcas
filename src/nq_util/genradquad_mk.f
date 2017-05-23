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
      Subroutine GenRadQuad_MK(R,nR,nR_Eff,m,Alpha,iNQ)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "debug.fh"
      Real*8 R(2,nR-1), Alpha, m
*
*---- Last point at infinity is eliminated
*
      If (Debug) Then
         Write (6,*) 'Log3 Algorithm (Mura-Knowles)'
         Write (6,*) 'Alpha,m=',Alpha,m
         Write (6,*) 'nR=',nR
      End If
      Do iR = 1, nR-1
         x       = DBLE(iR)/DBLE(nR)
         R(1,iR) = - Alpha * log( One - x**m )
         R(2,iR) = R(1,iR)**2 * Alpha * m * x**(m-One)
     &           / ( One - x**m ) / DBLE(nR)
      End Do
      nR_Eff = nR-1
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(iNQ)
      End
