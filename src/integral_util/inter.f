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
      SubRoutine Inter(iSet1,nSet1,iSet2,nSet2,iInter,nInter)
************************************************************************
*                                                                      *
* Object : to form the intersection of two sets.                       *
*                                                                      *
* Called from: NucAtt                                                  *
*              TwoEl                                                   *
*                                                                      *
* Calling    : None                                                    *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             February '90                                             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Integer iSet1(0:nSet1-1),iSet2(0:nSet2-1), iInter(0:7)
*
      nInter = 0
      Do 10 i1 = 0, nSet1-1
         Do 20 i2 = 0, nSet2-1
            If (iSet1(i1).eq.iSet2(i2)) Then
               nInter = nInter + 1
               iInter(nInter-1) = iSet1(i1)
               Go To 10
            End If
 20      Continue
 10   Continue
*
      Return
      End
