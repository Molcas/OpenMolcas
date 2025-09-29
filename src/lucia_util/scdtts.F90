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
! Copyright (C) 1995, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine SCDTTS(BLOCKS,IBLOCK,NBLOCK,NSMST,NSASO,NSBSO,IDC)
! Scale batch of
! blocks between determinant and combination form
!
! dets to combs
!
! The blocks are assumed to be in packed form !!
!
! Jeppe Olsen, August 1995

use Index_Functions, only: nTri_Elem
use Constants, only: Two, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
real(kind=wp), intent(inout) :: BLOCKS(*)
integer(kind=iwp), intent(in) :: NBLOCK, IBLOCK(8,NBLOCK), NSMST, NSASO(NSMST,*), NSBSO(NSMST,*), IDC
integer(kind=iwp) :: IASM, IATP, IBSM, IBTP, IOFFP, IPACK, JBLOCK, NELMNT, NIA, NIB
real(kind=wp), parameter :: SQ2 = sqrt(Two), SQ2I = sqrt(Half)

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' ======================='
write(u6,*) ' Information from SCDTTS'
write(u6,*) ' ======================='
write(u6,*) ' Input vector'
call WRTTTS(BLOCKS,IBLOCK,NBLOCK,NSMST,NSASO,NSBSO,2)
#endif

do JBLOCK=1,NBLOCK

  IATP = IBLOCK(1,JBLOCK)
  IBTP = IBLOCK(2,JBLOCK)
  IASM = IBLOCK(3,JBLOCK)
  IBSM = IBLOCK(4,JBLOCK)
  IOFFP = IBLOCK(6,JBLOCK)
  if (IBLOCK(1,JBLOCK) > 0) then
    ! Is this block diagonal in packed form
    if ((IASM == IBSM) .and. (IATP == IBTP)) then
      IPACK = 1
    else
      IPACK = 0
    end if
    NIA = NSASO(IASM,IATP)
    NIB = NSBSO(IBSM,IBTP)
    if (IPACK == 1) then
      NELMNT = nTri_Elem(NIA)
    else
      NELMNT = NIA*NIB
    end if
    ! Ms combinations
    if (IDC == 2) then
      BLOCKS(IOFFP:IOFFP+NELMNT-1) = SQ2*BLOCKS(IOFFP:IOFFP+NELMNT-1)
      if (IPACK == 1) call SCLDIA(BLOCKS(IOFFP),SQ2I,NIA,1)
    end if

  end if
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' Output vector'
call WRTTTS(BLOCKS,IBLOCK,NBLOCK,NSMST,NSASO,NSBSO,2)
#endif

end subroutine SCDTTS
