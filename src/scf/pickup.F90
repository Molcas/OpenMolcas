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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!***********************************************************************

subroutine PickUp(Tri,Vec,n)
!***********************************************************************
!                                                                      *
!     purpose: Pick up diagonal elements from triangular matrix and    *
!              put it to a vector                                      *
!                                                                      *
!     input:                                                           *
!       Tri     : triangular matrix                                    *
!       n       : dimension                                            *
!                                                                      *
!     output:                                                          *
!       Vec     : vector                                               *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: Tri(nTri_Elem(n))
real(kind=wp), intent(out) :: Vec(n)
integer(kind=iwp) :: i, ij

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

ij = 0
do i=1,n
  ij = ij+i
  Vec(i) = Tri(ij)
end do

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
return

end subroutine PickUp
