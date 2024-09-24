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

implicit none
integer n
real*8 A(n,n), S(n*(n+1)/2)
integer i, j, ij

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

ij = 0
do i=1,n
  do j=1,i
    ij = ij+1
    S(ij) = (A(i,j)+A(j,i))/2.0d+00
  end do
end do

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
return

end subroutine Sym
