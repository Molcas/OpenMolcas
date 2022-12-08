!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Fold_tMat(nSym,nBas,A,B)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Fold up symmetry blocked matrix A (UT-storage)                   *
!     into triangular matrices with scaled (by two) off-diagonals      *
!     The input and output matrices may be                             *
!     identical.                                                       *
!                                                                      *
!     calling arguments:                                               *
!     nSym    : input, integer                                         *
!               number of symmetry blocks                              *
!     nBas    : input, array of integers                               *
!               matrix dimension per symmetry block                    *
!     A       : input, array of real*8                                 *
!               Unfolded input matrix (UT-storage)                     *
!     B       : output, array of real*8                                *
!               Folded output matrix  (UT-storage)                     *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem
use Constants, only: Two
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym)
real(kind=wp), intent(in) :: A(*)
real(kind=wp), intent(_OUT_) :: B(*)
integer(kind=iwp) :: i, iOff, iSym, j

!************************************
iOff = 0

do iSym=1,nSym

  do j=1,nBas(iSym)

    do i=j+1,nBas(iSym)

      B(iOff+iTri(i,j)) = two*A(iOff+iTri(i,j))

    end do

    B(iOff+nTri_Elem(j)) = A(iOff+nTri_Elem(j))

  end do

  iOff = iOff+nTri_Elem(nBas(iSym))

end do

!ij = 0
!do isym=1,nsym
!  do i=1,nbas(isym)
!    do j=1,i-1
!      ij = ij+1
!      B(ij) = Two*A(ij)
!    end do
!    ij = ij+1
!    B(ij) = A(ij)
!  end do
!end do

return

end subroutine Fold_tMat
