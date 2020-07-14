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
      Subroutine Setup_OffAO()
      use Basis_Info, only: nCnttp
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
*
      Do iCnttp = 1, nCnttp
         lComp = 0
         lSh= 0
         kShStr = ipVal(iCnttp)
         kShEnd = ipVal(iCnttp)+nVal_Shells(iCnttp)-1
         Do kSh = kShStr, kShEnd
            If (Prjct(kSh)) Then
               kComp = 2*lSh + 1
            Else
               kComp = (lSh+1)*(lSh+2)/2
            End If
            kOffAO(iCnttp,lSh) = lComp
            If (nBasis_Cntrct(kSh).ne.0.and.nExp(kSh).ne.0)
     &         lComp = lComp + kComp
            lSh = lSh + 1
         End Do
         lOffAO(iCnttp) = lComp
      End Do
*
      Return
      End
