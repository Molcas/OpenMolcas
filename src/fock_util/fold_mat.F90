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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************

subroutine Fold_Mat(nSym,nBas,A,B)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Fold up symmetry blocked matrix A stored as SQUARE blocks        *
!     into TRIANGULAR matrices. The matrix A can be non-symmetric.     *
!     The input and output matrices CANNOT BE IDENTICAL !!!            *
!                                                                      *
!     calling arguments:                                               *
!     nSym    : input, integer                                         *
!               number of symmetry blocks                              *
!     nBas    : input, array of integers                               *
!               matrix dimension per symmetry block                    *
!     A       : input, array of real*8  (square storage)               *
!               Unfolded input matrix                                  *
!     B       : output, array of real*8 (triangular storage)           *
!               Folded output matrix                                   *
!                                                                      *
!----------------------------------------------------------------------*
!     Author:   F. Aquilante                                           *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym)
real(kind=wp), intent(in) :: A(*)
real(kind=wp), intent(_OUT_) :: B(*)
integer(kind=iwp) :: i, iOff1, iOff2, iSym, j

!************************************
iOff1 = 0
iOff2 = 0

do iSym=1,nSym

  do j=1,nBas(iSym)

    B(iOff1+nTri_Elem(j)) = A(iOff2+nBas(iSym)*(j-1)+j)

    do i=j+1,nBas(iSym)

      B(iOff1+iTri(i,j)) = A(iOff2+nBas(iSym)*(j-1)+i)+A(iOff2+nBas(iSym)*(i-1)+j)

    end do

  end do

  iOff1 = iOff1+nTri_Elem(nBas(iSym))
  iOff2 = iOff2+nBas(iSym)**2

end do

return

end subroutine Fold_Mat
