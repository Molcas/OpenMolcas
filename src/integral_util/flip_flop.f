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
       Subroutine Flip_Flop(Primitive)
       Implicit Real*8 (a-h,o-z)
       Logical Primitive
#include "itmax.fh"
#include "info.fh"
#include "Basis_Mode_Parameters.fh"
#include "Basis_Mode.fh"
#include "WrkSpc.fh"
*
      Call IZero(MaxBas,iTabMx+1)
      Call IZero(MaxPrm,iTabMx+1)
      Do iCnttp = 1, nCnttp
         nTest = nVal_Shells(iCnttp)-1
         If (AuxShell(iCnttp) .and.
     &       iCnttp.eq.iCnttp_Dummy) nTest=-1
         Do iCnt = 1, nCntr(iCnttp)
*
           Do 200 iAng=0, iAngMx
               If (iAng.gt.nTest)      Go To 200
               iShll = ipVal(iCnttp) + iAng
               If (nExp(iShll).eq.0)   Go To 200
               If (nBasis_Cntrct(iShll).eq.0) Go To 200
*
*              Decontract only the ordinary basis sets!
*
               If (Primitive.and..Not.AuxShell(iShll)
     &                      .and..Not.FragShell(iShll)) Then
                  ipCff(iShll)=ipCff_Prim(iShll)
                  nBasis(iShll)=nExp(iShll)
               Else
                  ipCff(iShll)=ipCff_Cntrct(iShll)
                  nBasis(iShll)=nBasis_Cntrct(iShll)
               End If
               MaxPrm(iAng) = Max(MaxPrm(iAng),nExp(iShll))
               MaxBas(iAng) = Max(MaxBas(iAng),nBasis(iShll))
*
 200        Continue
         End Do
      End Do
*
      Return
      End
