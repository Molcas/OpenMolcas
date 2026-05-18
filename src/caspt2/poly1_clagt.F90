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

subroutine POLY1_CLagT(CI1,CI2,CLag1,CLag2,RDMEIG,Scal)
! PER-AAKE MALMQUIST, 92-12-07
! THIS PROGRAM CALCULATES THE 1-EL DENSITY
! MATRIX FOR A CASSCF WAVE FUNCTION.

use sguga, only: SGS
use stdalloc, only: mma_allocate, mma_deallocate
use definitions, only: wp, iwp
use caspt2_module, only: nConf
use caspt2_module, only: MxCI, iAdr10, cLab10

implicit none
real(kind=wp), intent(in) :: CI1(NCONF), CI2(NCONF), RDMEIG(*), Scal
real(kind=wp), intent(inout) :: CLag1(*), CLag2(*)
real(kind=wp), allocatable :: SGM1(:)
integer(kind=iwp) :: nLev, I

nLev = SGS%nLev

if (NLEV > 0) then
  call mma_allocate(SGM1,MXCI,Label='SGM1')
  call DENS1T_RPT2_CLag(CI1,CI2,SGM1,CLag1,CLag2,RDMEIG,Scal,nLev)
end if

! REINITIALIZE USE OF DMAT.
! The fields IADR10 and CLAB10 are kept in caspt2_module
! CLAB10 replaces older field called LABEL.
do I=1,64
  IADR10(I,1) = -1
  IADR10(I,2) = 0
  CLAB10(I) = '   EMPTY'
end do
IADR10(1,1) = 0
! HENCEFORTH, THE CALL PUT(NSIZE,LABEL,ARRAY) WILL ENTER AN
! ARRAY ON LUDMAT AND UPDATE THE TOC.
if (NLEV > 0) call mma_deallocate(SGM1)

return

end subroutine POLY1_CLagT
