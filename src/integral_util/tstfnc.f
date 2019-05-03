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
* Copyright (C) 1990, IBM                                              *
*               1991, Roland Lindh                                     *
************************************************************************
      Logical Function TstFnc(iOper,nIrrep,iCoSet,nCoSet,iChTab,
     &                        iIrrep,iBsFnc,nStab)
************************************************************************
*                                                                      *
* Object: to establish if a function is a basis function of a          *
*         irreducible representation.                                  *
*                                                                      *
* Called from: Input                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              ICopy                                                   *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             September '91                                            *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Integer iOper(0:nIrrep-1), iCoSet(0:7,0:7), iAcc(0:7),
     &          iChTab(0:7,0:7)
#include "print.fh"
#include "real.fh"
      TstFnc = .True.
      Call iCopy(nCoSet,[0],0,iAcc,1)
*
*     Call qEnter('TstFnc')
*
*define _DEBUG_
#ifdef _DEBUG_
      Do i = 0, nCoSet-1
         Write (6,*) (iCoSet(i,j),j=0,nStab-1)
      End Do
      Write (6,*)
      Write (6,*) (iOper(i),i=0,nIrrep-1)
#endif
*
*     Loop over operators
*
      Do 10 i = 0, nIrrep-1
*
*        Find index of the generated center
*
         n = -1
         Do 20 j = 0, nCoSet-1
          If (n.lt.0) Then
            Do 21 k = 0, nStab-1
               If (iOper(i).eq.iCoSet(j,k)) n = j
 21         Continue
          End If
 20      Continue
*
         If (n.lt.0 .or. n.gt.nCoSet-1) Then
            Call WarningMessage(2,'TstFnc: n.lt.0 .or. n.gt.nCoSet-1')
            Write (6,*) ' Coset index',n,' is wrong!'
            Call Abend()
         End If
*
         iCom=iAnd(iOper(i),iBsFnc)
         iAcc(n) = iAcc(n) + iChTab(iIrrep,i)*iPrmt_(iCom)
*
 10   Continue
      Do 30 i = 0, nCoSet-1
         If (iAcc(i).eq.0) TstFnc = .False.
 30   Continue
*
*     Call qExit('TstFnc')
      Return
      End
      Integer Function iPrmt_(iCom)
************************************************************************
*     Returns the phase factor of a basis function under a symmetry    *
*     operation, jOper. iChct contains the information about the       *
*     character of the basis function.                                 *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      iPrmt_= 1
      Do 10 i = 1, 3
         If (iAnd(iCom,2**(i-1)).ne.0) iPrmt_= iPrmt_*(-1)
 10   Continue
      Return
      End
