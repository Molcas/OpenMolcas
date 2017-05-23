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
      Subroutine GenRadQuad_MHL(R,nR,nR_Eff,Alpha)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "debug.fh"
      Real*8 R(2,nR-1), Alpha
*
*---- Last point at infinity is eliminated
*
      If (Debug) Then
         Write (6,*) 'EM Algorithm (Murray, Handy, Laming)'
         Write (6,*) 'Alpha=',Alpha
         Write (6,*) 'nR=',nR
      End If
      Do iR = 1, nR-1
         x       = DBLE(iR)/DBLE(nR)
         R(1,iR) = Alpha * (x/(One-x))**2
         R(2,iR) = R(1,iR)**2 * Two* Alpha * x / (One-x)**3 / DBLE(nR)
      End Do
      nR_Eff = nR-1
*
      Return
      End
