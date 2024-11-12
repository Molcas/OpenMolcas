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

subroutine Binom(n,i,iBin)
!***********************************************************************
!                                                                      *
! Object: to compute the binomial factor                               *
!                                                                      *
! Called from: TraPAB                                                  *
!                                                                      *
! Calling    : None                                                    *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!***********************************************************************

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: n, i
integer(kind=iwp), intent(out) :: iBin
integer(kind=iwp) :: iDen, j, Num

Num = 1
iDen = 1
do j=1,i
  Num = Num*(n-j+1)
  iDen = iDen*j
end do
iBin = Num/iDen

return

end subroutine Binom
