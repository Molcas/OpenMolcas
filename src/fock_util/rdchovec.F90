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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************

!  RdChoVec
!
!> @brief
!>   Read dense Cholesky vectors in SO basis
!> @author T. B. Pedersen
!>
!> @details
!> Read dense SO Cholesky vectors \f$ L^J_{pq} \f$ where
!> \f$ J \f$ = \p iVec1, \p iVec1+1, ..., \p iVec1+nVec-1 and indices
!> \f$ p \ge q \f$.
!>
!> @note
!> Cholesky vectors must have been sorted into the
!> dense representation (by the Cholesky utility)
!> and the storage file must be open.
!>
!> @param[out] Vec   Array containing the vectors
!> @param[in]  nDim  Vector dimension
!> @param[in]  nVec  Number of vectors
!> @param[in]  iVec1 Index of first vector
!> @param[in]  Lunit Logical unit number
!***********************************************************************
subroutine RdChoVec(Vec,nDim,nVec,iVec1,Lunit)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDim, nVec, iVec1, Lunit
real(kind=wp), intent(out) :: Vec(nDim,nVec)
integer(kind=iwp) :: iOpt, iAdr, nTot

iOpt = 2
iAdr = nDim*(iVec1-1)+1
nTot = nDim*nVec
call dDaFile(Lunit,iOpt,Vec,nTot,iAdr)

end subroutine RdChoVec
