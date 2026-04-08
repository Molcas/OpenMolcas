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
! Copyright (C) 1991,1997, Jeppe Olsen                                 *
!               2015, Lasse Kragh Soerensen                            *
!***********************************************************************

subroutine RASSG3(CB,SB,NBATS,LBATS,I1BATS,IBATS,LUC,LUHC,I_AM_OUT,N_ELIMINATED_BATCHES)
! Direct RAS routine employing combined MOC/n-1 resolution method
!
! Jeppe Olsen   Winter of 1991
!               May 1997 : Connected to SBLOCK
!
! Lasse Soerensen October 2015
!                 Do not calculate unwanted batches for highly
!                 excited states.

use lucia_data, only: IDISK
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: CB(*)
real(kind=wp), intent(_OUT_) :: SB(*)
integer(kind=iwp), intent(in) :: NBATS, LBATS(*), I1BATS(*), IBATS(8,*), LUC, LUHC, I_AM_OUT(*), N_ELIMINATED_BATCHES
integer(kind=iwp) :: DUM(1), I, I_AM_NOT_WANTED, IEND, ILEN, IOFF, ISBLK, ISBOFF, ISTA, JBATS, NSB
integer(kind=iwp), allocatable :: SBOFF(:), SBSIZ(:)

#ifdef _DEBUGPRINT_
write(u6,*) ' ================='
write(u6,*) ' RASSG3 speaking :'
write(u6,*) ' ================='
write(u6,*) ' RASSG3 : NBATS = ',NBATS
#endif

!SVC: Compute offsets of a sigma batch in the sigma array.
!     The batches used inside sblock(s) use a batch size corresponding
!     to the 'expanded form' as computed inside part_civ2. This is
!     stored inside 7th element of IBATS. Later, the size that needs to
!     be actually written to disc uses the 'packed form', stored inside
!     the 8th element of IBATS. This also computes the total size NSB.

call mma_allocate(SBSIZ,NBATS,Label='SBSIZ')
call mma_allocate(SBOFF,NBATS,Label='SBOFF')

NSB = 0
do JBATS=1,NBATS
  ISTA = I1BATS(JBATS)
  IEND = I1BATS(JBATS)+LBATS(JBATS)-1
  SBSIZ(JBATS) = sum(IBATS(7,ISTA:IEND))
  NSB = NSB+SBSIZ(JBATS)
end do
SBOFF(1) = 1
do JBATS=2,NBATS
  SBOFF(JBATS) = SBOFF(JBATS-1)+SBSIZ(JBATS-1)
end do

!SVC: the entire sigma array is zeroed here, because each process will
!     zero only its own sigma blocks, and we need to do a global sum
!     operations later to combine blocks before writing.
SB(1:NSB) = Zero

do JBATS=1,NBATS

  ! Lasse addition start
  ! MGD here we try to remove the whole batch if possible
  ! later, we will put to zero individual blocks in case
  ! the batch had a mix of maximum and non-maximum occupation
  ! and thus survived this test

  I_AM_NOT_WANTED = 0
  do ISBLK=I1BATS(JBATS),I1BATS(JBATS)+LBATS(JBATS)-1
    I_AM_NOT_WANTED = 0
    do I=1,N_ELIMINATED_BATCHES
      if (I_AM_OUT(I) == ISBLK) then
        I_AM_NOT_WANTED = 1
        exit
      end if
    end do
    if (I_AM_NOT_WANTED == 0) exit
  end do
  if (I_AM_NOT_WANTED == 1) cycle

  ! Lasse addition end


  ISBOFF = SBOFF(JBATS)
  ! Obtain sigma for batch of blocks
  call SBLOCK(LBATS(JBATS),IBATS(1,I1BATS(JBATS)),1,CB,SB(ISBOFF),LUC,0,0,0,0,0)

end do

call GADSUM(SB,NSB)
!SVC: Write sigma array to disk here, after sum reduction.
!     The writing is done in consecutive blocks, but since I don't know
!     if this block structure is used internally, I didn't optimize this.
if (LUHC > 0) IDISK(LUHC) = 0
do JBATS=1,NBATS
  ISBOFF = SBOFF(JBATS)
  do ISBLK=I1BATS(JBATS),I1BATS(JBATS)+LBATS(JBATS)-1
    IOFF = IBATS(6,ISBLK)
    ILEN = IBATS(8,ISBLK)
    DUM(1) = ILEN
    call ITODS(DUM,1,-1,LUHC)
    !MGD zero afterwards since it is easier
    I_AM_NOT_WANTED = 0
    do I=1,N_ELIMINATED_BATCHES
      if (I_AM_OUT(I) == ISBLK) then
        I_AM_NOT_WANTED = 1
        exit
      end if
    end do
    if (I_AM_NOT_WANTED == 1) SB(ISBOFF-1+IOFF:ISBOFF-1+IOFF+ILEN-1) = Zero

    call TODSC(SB(ISBOFF-1+IOFF),ILEN,-1,LUHC)
  end do
end do

call mma_deallocate(SBSIZ)
call mma_deallocate(SBOFF)

DUM(1) = -1
call ITODS(DUM,1,-1,LUHC)

#ifdef _DEBUGPRINT_
write(u6,*) ' Final S-vector on disc'
call WRTVCD(SB,LUHC,1,-1)
#endif

end subroutine RASSG3
