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
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      Subroutine CmpInt(XInt,nInt,nBas,nIrrep,Label)
************************************************************************
*                                                                      *
* Object: to remove the offdiagonal nonzero blocks of matrix elements  *
*         for an operator.                                             *
*                                                                      *
*         XInt(1:nInt):array with nonzero elements                     *
*                                                                      *
*         nBas(0:nIrrep-1):number of basis functions in each irrep     *
*                                                                      *
*         Label: symmetry label of the operator for which the          *
*                matrix elements where computed.                       *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    : QEnter                                                  *
*              DCopy  (ESSL)                                           *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
*             University of Lund, SWEDEN                               *
*             March 1991                                               *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Real*8 XInt(nInt+4)
      Integer nBas(0:nIrrep-1)
*
      iCmp = 1
      iExp = 1
      Do 10 iIrrep = 0, nIrrep-1
         Do 20 jIrrep = 0, iIrrep
            ij = iEor(iIrrep,jIrrep)
            If (iAnd(Label,2**ij).eq.0) Go To 20
            If (iIrrep.eq.jIrrep) Then
               Len = nBas(iIrrep)*(nBas(iIrrep)+1)/2
               Do 30 iLen = 0, Len-1
                  XInt(iLen+iCmp) = Xint(iLen+iExp)
 30            Continue
               iCmp = iCmp + Len
               iExp = iExp + Len
            Else
               Len = nBas(iIrrep)*nBas(jIrrep)
               iExp = iExp + Len
            End If
 20      Continue
 10   Continue
      do 40 iadd=0,3
        XInt(iCmp+iadd)=XInt(iExp+iadd)
 40   continue
      nInt = iCmp-1
*
      Return
      End
