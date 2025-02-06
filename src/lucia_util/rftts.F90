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

subroutine RFTTS(BLOCKSI,BLOCKSO,IBLOCK,NBLOCK,ICOPY,NSMST,NSASO,NSBSO,IDC,PS,IWAY,IPRNT)
! Reformat between determinant and combination form of
! matrices. No scaling is performed .
!
! IWAY = 1 : dets to combs
! IWAY = 2 : combs to dets
!
! Combination storage mode is defined BY IDC
!
! Jeppe Olsen, August 1995

use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: BLOCKSI(*), BLOCKSO(*), PS
integer(kind=iwp) :: NBLOCK, IBLOCK(8,NBLOCK), ICOPY, NSMST, NSASO(NSMST,*), NSBSO(NSMST,*), IDC, IWAY, IPRNT
integer(kind=iwp) :: IASM, IATP, IBSM, IBTP, IOFFI, IOFFO, IPACK, ISCI, ISCO, JBLOCK, LENGTH, NELMNT, NIA, NIB, NTEST

NTEST = 0
NTEST = max(NTEST,IPRNT)

LENGTH = 0
if (IWAY == 1) then
  ISCI = 1
  ISCO = 2
else
  ISCI = 2
  ISCO = 1
end if

if (NTEST > 10) then
  write(u6,*) ' Information from RFTTS'
  write(u6,*) ' ======================'
  write(u6,*) ' Input vector'
  call WRTTTS(BLOCKSI,IBLOCK,NBLOCK,NSMST,NSASO,NSBSO,ISCI)
end if

do JBLOCK=1,NBLOCK

  IATP = IBLOCK(1,JBLOCK)
  IBTP = IBLOCK(2,JBLOCK)
  IASM = IBLOCK(3,JBLOCK)
  IBSM = IBLOCK(4,JBLOCK)
  if (IBLOCK(1,JBLOCK) > 0) then

    if (IWAY == 1) then
      IOFFI = IBLOCK(5,JBLOCK)
      IOFFO = IBLOCK(6,JBLOCK)
    else
      IOFFO = IBLOCK(5,JBLOCK)
      IOFFI = IBLOCK(6,JBLOCK)
    end if
    ! Is this block diagonal in packed form
    if ((IDC == 2) .and. (IASM == IBSM) .and. (IATP == IBTP)) then
      IPACK = 1
    else
      IPACK = 0
    end if
    NIA = NSASO(IASM,IATP)
    NIB = NSBSO(IBSM,IBTP)
    ! Number of elements in output block
    if ((IPACK == 1) .and. (ISCO == 2)) then
      NELMNT = NIA*(NIA+1)/2
    else
      NELMNT = NIA*NIB
    end if
    !write(u6,*) ' JBLOCK, NELMNT = ',JBLOCK,NELMNT
    !write(u6,*) ' RFTTS : IATP IBTP IASM IBSM ',IATP,IBTP,IASM,IBSM
    !write(u6,*) ' RFTTS : NIA NIB IOFFI,IOFFO',NIA,NIB,IOFFI,IOFFO

    if (IPACK == 0) then
      ! Just copy
      call COPVEC(BLOCKSI(IOFFI),BLOCKSO(IOFFO),NELMNT)
    else
      if (IWAY == 1) then
        ! unpacked => packed
        !    TRIPK3(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
        call TRIPK3(BLOCKSI(IOFFI),BLOCKSO(IOFFO),1,NIA,NIA,PS)
      else
        ! Packed => unpacked
        call TRIPK3(BLOCKSO(IOFFO),BLOCKSI(IOFFI),2,NIA,NIA,PS)
      end if
    end if
    LENGTH = LENGTH+NELMNT
  end if
end do

if (ICOPY /= 0) call COPVEC(BLOCKSO,BLOCKSI,LENGTH)

if (NTEST > 10) then
  write(u6,*) ' Information from RFTTS'
  write(u6,*) ' ======================'
  write(u6,*) ' Output vector'
  call WRTTTS(BLOCKSO,IBLOCK,NBLOCK,NSMST,NSASO,NSBSO,ISCO)
end if

end subroutine RFTTS
