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
! Copyright (C) 1996, Markus P. Fuelscher                              *
!***********************************************************************

subroutine Fold(nSym,nBas,A,B)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Fold up symmetry blocked matrix A stored as square blocks        *
!     into triangular matrices and scale the off-diagonal elements     *
!     by a factor of 2 (two). The input and output matrices may be     *
!     identical.                                                       *
!                                                                      *
!     calling arguments:                                               *
!     nSym    : input, integer                                         *
!               number of symmetry blocks                              *
!     nBas    : input, array of integers                               *
!               matrix dimension per symmetry block                    *
!     A       : input, array of real*8                                 *
!               Unfolded input matrix                                  *
!     B       : output, array of real*8                                *
!               Folded output matrix                                   *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1996                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Constants, only: Two
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym)
real(kind=wp), intent(in) :: A(*)
real(kind=wp), intent(_OUT_) :: B(*)
integer(kind=iwp) :: iBas, iOff1, iOff2, iSym, jBas, mBas

iOff1 = 0
iOff2 = 0
do iSym=1,nSym
  mBas = nBas(iSym)
  do iBas=1,mBas
    do jBas=1,iBas-1
      B(iOff2+jBas) = Two*A(iOff1+jBas)
    end do
    B(iOff2+iBas) = A(iOff1+iBas)
    iOff1 = iOff1+mBas
    iOff2 = iOff2+iBas
  end do
end do

return

end subroutine Fold
