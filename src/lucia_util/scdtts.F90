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

subroutine SCDTTS(BLOCKS,IBLOCK,NBLOCK,NSMST,NSASO,NSBSO,IDC,IWAY,IPRNT)
! Scale batch of
! blocks between determinant and combination form
!
! IWAY = 1 : dets to combs
! IWAY = 2 : combs to dets
!
! The blocks are assumed to be in packed form !!
!
!. Jeppe Olsen, August 1995

implicit real*8(A-H,O-Z)
! General input
dimension NSASO(NSMST,*), NSBSO(NSMST,*)
dimension BLOCKS(*)
integer IBLOCK(8,NBLOCK)
!logical DIAGBL

NTEST = 0
NTEST = max(NTEST,IPRNT)
if (NTEST > 10) then
  write(6,*)
  write(6,*) ' ======================='
  write(6,*) ' Information from SCDTTS'
  write(6,*) ' ======================='
  write(6,*) ' Input vector'
  call WRTTTS(BLOCKS,IBLOCK,NBLOCK,NSMST,NSASO,NSBSO,2)
end if

SQ2 = sqrt(2.0d0)
SQ2I = 1.0d0/SQ2

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
      NELMNT = NIA*(NIA+1)/2
    else
      NELMNT = NIA*NIB
    end if
    ! Ms combinations
    if (IDC == 2) then
      if (IWAY == 1) then
        FACTOR = SQ2
      else
        FACTOR = SQ2I
      end if
      call SCALVE(BLOCKS(IOFFP),FACTOR,NELMNT)
      if (IPACK == 1) then
        FACTOR = 1.0d0/FACTOR
        call SCLDIA(BLOCKS(IOFFP),FACTOR,NIA,1)
      end if
    end if

  end if
end do

if (NTEST >= 10) then
  write(6,*) ' Output vector'
  call WRTTTS(BLOCKS,IBLOCK,NBLOCK,NSMST,NSASO,NSBSO,2)
end if

end subroutine SCDTTS
