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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************

!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*
subroutine MKMAW(SGS)

use gugx, only: SGStruct
use stdalloc, only: mma_allocate
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
type(SGStruct), intent(inout) :: SGS
integer(kind=iwp) :: IC, ID, ISUM, IU, IV

call mma_allocate(SGS%MAW,[1,SGS%nVert],[0,3],Label='SGS%MAW')

! COPY LOWER PART OF DIRECT ARC WEIGHT TABLE INTO MAW:
SGS%MAW(SGS%MVSta:SGS%nVert,0:3) = SGS%DAW(SGS%MVSta:SGS%nVert,0:3)
! COPY UPPER PART OF REVERSE ARC WEIGHT TABLE INTO MAW. HOWEVER,
!    NOTE THAT THE MAW TABLE IS ACCESSED BY THE UPPER VERTEX.
SGS%MAW(1:SGS%MVSta-1,0:3) = 0
do IU=1,SGS%MVSta-1
  do IC=0,3
    ID = SGS%Down(IU,IC)
    if (ID /= 0) SGS%MAW(IU,IC) = SGS%RAW(ID,IC)
  end do
end do
! FINALLY, ADD AN OFFSET TO ARCS LEADING TO MIDLEVELS:
ISUM = 1
do IV=SGS%MVSta,SGS%MVEnd
  do IC=0,3
    IU = SGS%Up(IV,IC)
    if (IU == 0) cycle
    SGS%MAW(IU,IC) = ISUM+SGS%MAW(IU,IC)
  end do
  ISUM = ISUM+SGS%RAW(IV,4)
end do
do IV=SGS%MVSta,SGS%MVEnd
  do IC=0,3
    if (SGS%Down(IV,IC) == 0) cycle
    SGS%MAW(IV,IC) = ISUM+SGS%MAW(IV,IC)
  end do
  ISUM = ISUM+SGS%DAW(IV,4)
end do
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' THE MODIFIED ARC WEIGHT TABLE IN MKMAW:'
do IV=1,SGS%nVert
  write(u6,'(1X,I4,5X,5(1X,I6))') IV,SGS%MAW(IV,0:3)
end do
write(u6,*)
#endif

end subroutine MKMAW
