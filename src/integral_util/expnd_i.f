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
      SubRoutine Expnd_i(Array,n,m)
************************************************************************
*                                                                      *
* Object: to do an in place expansion of a triagularized matrix.       *
*                                                                      *
* Called from: Cntrct                                                  *
*                                                                      *
* Calling    : DCopy  (ESSL)                                           *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             May '90                                                  *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Real*8 Array(m,n*n)
      Call qEnter('Expnd_i')
*
      nij = n*(n+1)/2
      Do 10 i = n, 1, -1
         Do 20 j = n, i+1, -1
            ji = n*(i-1) + j
            ij = n*(j-1) + i
            If (nij.ne.ij) call dcopy_(m,Array(1,nij),1,Array(1,ij),1)
            If (nij.ne.ji) call dcopy_(m,Array(1,nij),1,Array(1,ji),1)
            nij = nij - 1
 20      Continue
         ii = n*(i-1) + i
         If (nij.ne.ii) call dcopy_(m,Array(1,nij),1,Array(1,ii),1)
         nij = nij - 1
 10   Continue
*
      Call qExit('Expnd_i')
      Return
      End
