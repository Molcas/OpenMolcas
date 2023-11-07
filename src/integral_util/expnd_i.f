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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************
      SubRoutine Expnd_i(Array,n,m)
!***********************************************************************
!                                                                      *
! Object: to do an in place expansion of a triagularized matrix.       *
!                                                                      *
! Called from: Cntrct                                                  *
!                                                                      *
! Calling    : DCopy  (ESSL)                                           *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             May '90                                                  *
!***********************************************************************
      Implicit None
      Integer m, n
      Real*8 Array(m,n*n)

      Integer nij, i, j, ij, ji, ii
!
      nij = n*(n+1)/2
      Do 10 i = n, 1, -1
         Do 20 j = n, i+1, -1
            ji = n*(i-1) + j
            ij = n*(j-1) + i
            If (nij.ne.ij) Array(:,ij)=Array(:,nij)
            If (nij.ne.ji) Array(:,ji)=Array(:,nij)
            nij = nij - 1
 20      Continue
         ii = n*(i-1) + i
         If (nij.ne.ii) Array(:,ii)=Array(:,nij)
         nij = nij - 1
 10   Continue
!
      Return
      End SubRoutine Expnd_i
