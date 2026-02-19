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
! Copyright (C) 2000, Markus P. Fuelscher                              *
!***********************************************************************

subroutine ClnMO(CMO)
!***********************************************************************
!                                                                      *
!     In order to preserve symmetry of orbitals which belong           *
!     to a point group of higher order than those allowed in           *
!     MOLCAS, it is sometimes necessary to enforce the MO coefficients *
!     to remain zero. This subroutine does that job in every           *
!     iteration.                                                       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 2000                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!     This subroutine replaces the one written earlier by L. Serrano   *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 2000                                 *
!                                                                      *
!***********************************************************************

use general_data, only: CleanMask, NSYM, NBAS
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: CMO(*)
integer(kind=iwp) :: i, ij, iSym, j

! Body

ij = 0
do iSym=1,nSym
  do i=1,nBas(iSym)
    do j=1,nBas(iSym)
      ij = ij+1
      if (CleanMask(ij) == 1) CMO(ij) = Zero
    end do
  end do
end do

end subroutine ClnMO
