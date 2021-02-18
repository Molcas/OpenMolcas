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
      SubRoutine PckInt_mck(abab,nZeta,nab,ab,Zeta)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 abab(nZeta,nab,nab), ab(nZeta), Zeta(nZeta)
*
*--------Integrals
*
*
*
      Do iZeta=1,nZeta
        xMax=0.0d0
        Do iab = 1, nab
            xTest = Abs(abab(iZeta,iab,iab))
            If (xTest.gt.xMax) Then
               xMax = xTest
            End If
        End Do
        ab(iZeta) = Sqrt(xMax)
      End Do
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Zeta)
      End
