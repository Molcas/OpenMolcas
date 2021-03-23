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
! Copyright (C) 1990,1992, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

subroutine Trnglr(Array,m,n)
!***********************************************************************
!                                                                      *
! Object: to do an in place expansion of a triangularized matrix.      *
!                                                                      *
! Called from: Tcrtnc                                                  *
!                                                                      *
! Calling    : DCopy  (ESSL)                                           *
!              DaXpY  (ESSL)                                           *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             May '90                                                  *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified to inplace triangularization, January '92.      *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "print.fh"
real*8 Array(m,n*n)

! Observe that the desymmetrization will not yield a symmetric
! result. In order to apply the triangularization we will have
! to symmetrize the matrix first.

do i=1,n
  do j=1,i-1
    mji = n*(i-1)+j
    mij = n*(j-1)+i
    call DaXpY_(m,One,Array(1,mij),1,Array(1,mji),1)
  end do
end do

do i=1,n
  do j=1,i
    mji = n*(i-1)+j
    nij = i*(i-1)/2+j
    if (nij /= mji) call dcopy_(m,Array(1,mji),1,Array(1,nij),1)
  end do
end do

return

end subroutine Trnglr
