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
subroutine RFTTS(BLOCKSI,BLOCKSO,IBLOCK,NBLOCK,NSMST,NSASO,NSBSO,IDC)
! Reformat between determinant and combination form of matrices.
! No scaling is performed.
!
! dets to combs
!
! Combination storage mode is defined BY IDC
!
! Jeppe Olsen, August 1995

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: BLOCKSI(*)
real(kind=wp), intent(_OUT_) :: BLOCKSO(*)
integer(kind=iwp), intent(in) :: NBLOCK, IBLOCK(8,NBLOCK), NSMST, NSASO(NSMST,*), NSBSO(NSMST,*), IDC
integer(kind=iwp) :: IASM, IATP, IBSM, IBTP, IOFFI, IOFFO, IPACK, JBLOCK, LENGTH, NELMNT, NIA, NIB

LENGTH = 0

#ifdef _DEBUGPRINT_
write(u6,*) ' Information from RFTTS'
write(u6,*) ' ======================'
write(u6,*) ' Input vector'
call WRTTTS(BLOCKSI,IBLOCK,NBLOCK,NSMST,NSASO,NSBSO,1)
#endif

do JBLOCK=1,NBLOCK

  IATP = IBLOCK(1,JBLOCK)
  IBTP = IBLOCK(2,JBLOCK)
  IASM = IBLOCK(3,JBLOCK)
  IBSM = IBLOCK(4,JBLOCK)
  if (IBLOCK(1,JBLOCK) > 0) then

    IOFFI = IBLOCK(5,JBLOCK)
    IOFFO = IBLOCK(6,JBLOCK)
    ! Is this block diagonal in packed form
    if ((IDC == 2) .and. (IASM == IBSM) .and. (IATP == IBTP)) then
      IPACK = 1
    else
      IPACK = 0
    end if
    NIA = NSASO(IASM,IATP)
    NIB = NSBSO(IBSM,IBTP)
    ! Number of elements in output block
    if (IPACK == 1) then
      NELMNT = nTri_Elem(NIA)
    else
      NELMNT = NIA*NIB
    end if
    !write(u6,*) ' JBLOCK, NELMNT = ',JBLOCK,NELMNT
    !write(u6,*) ' RFTTS : IATP IBTP IASM IBSM ',IATP,IBTP,IASM,IBSM
    !write(u6,*) ' RFTTS : NIA NIB IOFFI,IOFFO',NIA,NIB,IOFFI,IOFFO

    if (IPACK == 0) then
      ! Just copy
      BLOCKSO(IOFFO:IOFFO+NELMNT-1) = BLOCKSI(IOFFI:IOFFI+NELMNT-1)
    else
      ! unpacked => packed
      !    TRIPK31(AUTPAK,APAK,MATDIM,NDIM)
      call TRIPK31(BLOCKSI(IOFFI),BLOCKSO(IOFFO),NIA,NIA)
    end if
    LENGTH = LENGTH+NELMNT
  end if
end do

BLOCKSI(1:LENGTH) = BLOCKSO(1:LENGTH)

#ifdef _DEBUGPRINT_
write(u6,*) ' Information from RFTTS'
write(u6,*) ' ======================'
write(u6,*) ' Output vector'
call WRTTTS(BLOCKSO,IBLOCK,NBLOCK,NSMST,NSASO,NSBSO,2)
#endif

end subroutine RFTTS
