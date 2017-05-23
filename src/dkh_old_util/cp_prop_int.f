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
      Subroutine Cp_Prop_Int(AInt,nAInt,BInt,nBInt,nBas,nIrrep,Label)
************************************************************************
*                                                                      *
* Object: replace the diagonal blocks of the property integrals.       *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Real*8 AInt(nAInt+4), BInt(nBInt+4)
      Integer nBas(0:nIrrep-1)
*
      iCmp = 1
      iExp = 1
      Do iIrrep = 0, nIrrep-1
         Do 20 jIrrep = 0, iIrrep
            ij = iEor(iIrrep,jIrrep)
            If (iAnd(Label,2**ij).eq.0) Go To 20
            If (iIrrep.eq.jIrrep) Then
*
               Len = nBas(iIrrep)*(nBas(iIrrep)+1)/2
               Do iLen = 0, Len-1
C                 Write (*,*) AInt(iExp+iLen), BInt(iCmp+iLen)
                  AInt(iExp+iLen) = BInt(iCmp+iLen)
               End Do
               iCmp = iCmp + Len
               iExp = iExp + Len
            Else
               Len = nBas(iIrrep)*nBas(jIrrep)
               iExp = iExp + Len
            End If
 20      Continue
      End Do
*
      Return
      End
