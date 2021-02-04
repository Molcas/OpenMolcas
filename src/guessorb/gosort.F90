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
!***********************************************************************
!                                                                      *
!     purpose: Sort the set of eigenvalues and eigenvectors            *
!                                                                      *
!                                                                      *
!     input:                                                           *
!       EVal    : the set of eigenvalues in random order               *
!       EVec    : the set of eigenvectors in random order              *
!       n,nB    : dimensions                                           *
!                                                                      *
!     output:                                                          *
!       EVal    : sorted set of eigenvalues                            *
!       EVec    : sorted set of eigenvectors                           *
!                                                                      *
!     called from: DCore, NewOrb, IvoGen                               *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
!     University of Lund, Sweden, 1992                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

subroutine goSort(EVal,EVec,n,nB)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n, nB
real(kind=wp), intent(inout) :: Eval(n), EVec(nB,n)
integer(kind=iwp) :: i, j, k, l
real(kind=wp) :: Swap

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

do i=1,n-1
  k = i
  do j=i+1,n
    if (EVal(j) < EVal(k)) k = j
  end do
  if (k /= i) then
    Swap = EVal(k)
    EVal(k) = EVal(i)
    EVal(i) = Swap
    do l=1,nB
      Swap = EVec(l,k)
      EVec(l,k) = -EVec(l,i)
      EVec(l,i) = Swap
    end do
  end if
end do

end subroutine goSort
