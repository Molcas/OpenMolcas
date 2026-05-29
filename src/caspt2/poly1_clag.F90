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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine POLY1_CLag(NCONF,NLEV,CI,CLag,RDMEIG)
! PER-AAKE MALMQUIST, 92-12-07
! THIS PROGRAM CALCULATES THE 1-EL DENSITY
! MATRIX FOR A CASSCF WAVE FUNCTION.

use caspt2_module, only: cLab10, iAdr10, MxCI
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NCONF, NLEV
real(kind=wp), intent(in) :: CI(NCONF), RDMEIG(NLEV**2)
real(kind=wp), intent(inout) :: CLag(NCONF)
real(kind=wp), allocatable :: SGM1(:)

if (NLEV > 0) then
  call MMA_ALLOCATE(SGM1,MXCI,LABEL='SGM1')
  call DENS1_RPT2_CLag(CI,NCONF,SGM1,MXCI,CLag,NCONF,RDMEIG,nLev)
end if
!return !! for test purpose

! REINITIALIZE USE OF DMAT.
! The fields IADR10 and CLAB10 are kept in caspt2_module
! CLAB10 replaces older field called LABEL.
IADR10(:,1) = -1
IADR10(:,2) = 0
CLAB10(:) = '   EMPTY'
IADR10(1,1) = 0
! HENCEFORTH, THE CALL PUT(NSIZE,LABEL,ARRAY) WILL ENTER AN
! ARRAY ON LUDMAT AND UPDATE THE TOC.
if (NLEV > 0) call MMA_DEALLOCATE(SGM1)

return

end subroutine POLY1_CLag
