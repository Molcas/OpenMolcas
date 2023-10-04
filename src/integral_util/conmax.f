!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************
      SubRoutine ConMax(A,nPrim,mPrim,B,nCont,C,mCont)
!***********************************************************************
!                                                                      *
! Object: to find the largest element in the contraction matrix  for   *
!         each primitive index.                                        *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             July '91                                                 *
!***********************************************************************
      Implicit None
      Integer, Intent(In):: nPrim,nCont,mPrim,mCont
      Real*8, Intent(In):: B(nPrim,nCont), C(mPrim,mCont)
      Real*8, Intent(Out):: A(nPrim,mPrim)
!
      Integer iPrim, jPrim
      Real*8, External:: DDot_
      Real*8 Temp
      Do iPrim = 1, nPrim
         Temp=DDot_(nCont,B(iPrim,1),nPrim,B(iPrim,1),nPrim)
         Do jPrim = 1, mPrim
            A(iPrim,jPrim) = Temp
         End Do
      End Do
!
      Do jPrim = 1, mPrim
         Temp=DDot_(mCont,C(jPrim,1),mPrim,C(jPrim,1),mPrim)
         Do iPrim = 1, nPrim
            A(iPrim,jPrim) = Sqrt(A(iPrim,jPrim)*Temp)
         End Do
      End Do
!
      Return
      End SubRoutine ConMax
