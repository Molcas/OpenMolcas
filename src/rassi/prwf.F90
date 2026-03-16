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

subroutine PRWF(SGS,CIS,ISYCI,CI,CITHR)

use definitions, only: iwp, wp
use gugx, only: SGStruct, CIStruct
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
integer(kind=iwp) ISYCI
real(kind=wp), intent(in) :: CI(*), CITHR
integer(kind=iwp), allocatable :: ICS(:)
integer(kind=iwp) NLEV, NMIDV

NLEV = SGS%nLev
NMIDV = CIS%nMidV

call mma_allocate(ICS,NLEV,Label='ICS')
call PRWF1(SGS,CIS,NLEV,NMIDV,SGS%ISM,ICS,CIS%NOCSF,CIS%IOCSF,CIS%NOW,CIS%IOW,ISYCI,CI,CITHR)
call mma_deallocate(ICS)

end subroutine PRWF
