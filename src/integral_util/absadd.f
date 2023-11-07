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
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************
      SubRoutine AbsAdd(nVec,Vec1,inc1,Vec2,inc2)
!***********************************************************************
!                                                                      *
! Object: to add the absolute values of a vector to another vector.    *
!         The square root due to the Cauchy-Schwatz equation.          *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             May '92.                                                 *
!***********************************************************************
      Implicit None
      Integer, Intent(In):: nVec, inc1, inc2
      Real*8, Intent(InOut):: Vec1(1+(nVec-1)*Inc1)
      Real*8, Intent(InOut):: Vec2(1+(nVec-1)*Inc2)

      Integer  iVec
!
      Do iVec = 1, nVec
         Vec2(1+(iVec-1)*Inc2) = Vec2(1+(iVec-1)*Inc2) +
     &                  Sqrt(Abs(Vec1(1+(iVec-1)*Inc1)))
      End Do
!
      Return
      End SubRoutine AbsAdd
