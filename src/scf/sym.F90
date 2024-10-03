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

subroutine Sym(A,S,n)
!***********************************************************************
!                                                                      *
!     purpose: Symmetrize matrix                                       *
!                                                                      *
!     input:                                                           *
!       A       : input matrix (square)                                *
!       n       : dimension                                            *
!                                                                      *
!     output:                                                          *
!       S       : Symmetrized matrix (triangular)                      *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: A(n,n)
real(kind=wp), intent(out) :: S(nTri_Elem(n))
integer(kind=iwp) :: i, ij, j

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

ij = 0
do i=1,n
  do j=1,i
    ij = ij+1
    S(ij) = (A(i,j)+A(j,i))*Half
  end do
end do

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
return

end subroutine Sym
