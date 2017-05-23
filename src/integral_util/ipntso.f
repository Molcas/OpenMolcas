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
      Integer Function iPntSO(j1,j2,lOper,nbas)
************************************************************************
*                                                                      *
* Object: to compute the offset to an one-electron integral block.     *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             February '91                                             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Integer nbas(0:7)
*
      iPntSO=0
      iSmLbl = lOper
      Do iIrrep = 0, j1
         jMax = iIrrep
         If (iIrrep.eq.j1) jMax = j2-1
         Do 20 jIrrep = 0, jMax
            ij = iEor(iIrrep,jIrrep)
            If (iAnd(iSmLbl,2**ij).eq.0) Go To 20
            If (iIrrep.eq.jIrrep) Then
               iPntSO = iPntSO + nBas(iIrrep)*(nBas(iIrrep)+1)/2
            Else
               iPntSO = iPntSO + nBas(iIrrep)*nBas(jIrrep)
            End If
 20      Continue
      End Do
*
      Return
      End
