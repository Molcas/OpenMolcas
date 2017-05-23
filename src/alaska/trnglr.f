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
* Copyright (C) 1990,1992, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      SubRoutine Trnglr(Array,m,n)
************************************************************************
*                                                                      *
* Object: to do an in place expansion of a triangularized matrix.      *
*                                                                      *
* Called from: Tcrtnc                                                  *
*                                                                      *
* Calling    : DCopy  (ESSL)                                           *
*              DaXpY  (ESSL)                                           *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             May '90                                                  *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             Modified to inplace triangularization, January '92.      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
      Real*8 Array(m,n*n)
*     Call qEnter('Trnglr')
*
*     Observe that the desymmetrization will not yield a symmetric
*     result. In order to apply the triangularization we will have
*     to symmetrize the matrix first.
*
      Do 100 i = 1, n
         Do 200 j = 1, i-1
            mji = n*(i-1)+j
            mij = n*(j-1)+i
            Call DaXpY_(m,One,Array(1,mij),1,
     &                       Array(1,mji),1)
 200     Continue
 100  Continue
*
      Do 10 i = 1, n
         Do 20 j = 1, i
            mji = n*(i-1) + j
            nij = i*(i-1)/2 + j
            If (nij.ne.mji) call dcopy_(m,Array(1,mji),1,
     &                                   Array(1,nij),1)
*
 20      Continue
 10   Continue
*
*     Call qExit('Trnglr')
      Return
      End
