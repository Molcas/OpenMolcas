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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      SubRoutine Stblzr(iU,nU,iV,nV,iR,iM,nM)
************************************************************************
*                                                                      *
*  Object: to form the proper stabilizer for a pair of entities        *
*          given the respective stabilizers and the operator acting    *
*          on the second entity.                                       *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             June '90                                                 *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Integer iU(nU), iV(nV), iM(8)
      Logical UeqV, RinT_
*
*     See if U and V are the same
*
      UeqV = .True.
      Do 10 i = 1, nV
         If (.Not.RinT_(iU,nU,iV(i))) UeqV=.False.
 10   Continue
      Do 11 i = 1, nU
         If (.Not.RinT_(iV,nV,iU(i))) UeqV=.False.
 11   Continue
*
      If (UeqV) Then
*
*        M is formed as U union RU
*
         Call iCopy(nU,iU,1,iM,1)
         nM = nU
         Do 20 i = 1, nU
            iRU = iEor(iR,iU(i))
            If (.Not.RinT_(iM,nM,iRU)) Then
               nM = nM + 1
               iM(nM) = iRU
            End If
 20      Continue
      Else
*
*        M is formed as U intersection V
*
         Call Inter(iU,nU,iV,nV,iM,nM)
      End If
*
      Return
      End
