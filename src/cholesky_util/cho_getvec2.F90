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

subroutine CHO_GETVEC2(CHOVEC,LENVEC,NUMVEC,IVEC1,ISYM,SCR,LSCR)
!
! Purpose: read Cholesky vectors IVEC=IVEC1,....,IVEC1+NUMVEC-1
!          of symmetry ISYM from file. The vectors are returned
!          in the "current" reduced set. This routine attempts
!          to minimize gather/scatter operations along the lines
!          of cho_getvec1. However, in this version, buffering
!          is (hopefully) further improved by reading vectors
!          across reduced sets.
!
! NOTE: the scratch array SCR(LSCR) is used to read vectors from
!       disk and should not be smaller than NNBSTR(ISYM,1)+1,
!       preferably more.

use Cholesky, only: Cho_AdrVec, InfVec, iScr, nnBstR, nSys_call
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LENVEC, NUMVEC, IVEC1, ISYM, LSCR
real(kind=wp), intent(inout) :: CHOVEC(LENVEC,NUMVEC), SCR(LSCR)
integer(kind=iwp) :: IAB, ILOC, IMAPC, IOFF(0:1), IREDC, IVEC2, JRED, JRED1, JRED2, JVEC1, JVEC2, KJUNK, KOFF, KSCR, KVEC, KVEC1, &
                     LEFT, LNUM, LRED, LVEC, LVEC1, MUSED, NVRD
character(len=*), parameter :: SECNAM = 'CHO_GETVEC2'

! Some initializations.
! ---------------------

ILOC = 3

IVEC2 = IVEC1+NUMVEC-1

KJUNK = 1
KSCR = KJUNK+1
LEFT = LSCR-KSCR+1
if (LEFT < 1) call CHO_QUIT('Insufficient scratch space in '//SECNAM,101)

SCR(KJUNK) = Zero
IOFF(0) = KJUNK

! Start buffer batch loop.
! ------------------------

KVEC1 = 1
JVEC1 = IVEC1
IREDC = -1
IMAPC = -1
do while (JVEC1 <= IVEC2)

  ! Read as many vectors as fit into scratch space.
  ! -----------------------------------------------

  JRED1 = INFVEC(JVEC1,2,ISYM)
  NVRD = 0
  MUSED = 0
  call CHO_VECRD(SCR(KSCR),LEFT,JVEC1,IVEC2,ISYM,NVRD,IREDC,MUSED)
  if (CHO_ADRVEC == 1) then
    NSYS_CALL = NSYS_CALL+1
  else if (CHO_ADRVEC == 2) then
    NSYS_CALL = NSYS_CALL+NVRD
  else
    call CHO_QUIT('CHO_ADRVEC error in '//SECNAM,102)
  end if

  ! Quit if no vectors were read.
  ! -----------------------------

  if (NVRD < 1) call CHO_QUIT('Insufficient scratch space for read in '//SECNAM,101)

  ! Loop over reduced sets in scratch space.
  ! ----------------------------------------

  JVEC2 = JVEC1+NVRD-1
  JRED2 = INFVEC(JVEC2,2,ISYM)
  LVEC1 = JVEC1
  IOFF(1) = KSCR-1
  do JRED=JRED1,JRED2

    ! Count vectors read from this reduced set.
    ! -----------------------------------------

    LNUM = 0
    LVEC = LVEC1-1
    do while (LVEC < JVEC2)
      LVEC = LVEC+1
      LRED = INFVEC(LVEC,2,ISYM)
      if (LRED == JRED) then
        LNUM = LNUM+1 ! increase counter
      else
        LVEC = JVEC2 ! break loop
      end if
    end do

    if (LNUM > 0) then

      ! Read index arrays for this reduced set (if needed).
      ! ---------------------------------------------------

      if (JRED /= IREDC) then
        call CHO_GETRED(JRED,ILOC,.false.)
        call CHO_SETREDIND(ILOC)
        IREDC = JRED
      end if

      ! Set up rs-to-rs map (if needed).
      ! --------------------------------

      if (JRED /= IMAPC) then
        call CHO_RS2RS(ISCR,size(ISCR),2,3,JRED,ISYM)
        IMAPC = JRED
      end if

      ! Copy vectors to result array.
      ! -----------------------------

      do LVEC=1,LNUM
        KVEC = KVEC1+LVEC-1
        do IAB=1,NNBSTR(ISYM,2)
          KOFF = IOFF(min(ISCR(IAB),1))+ISCR(IAB)
          CHOVEC(IAB,KVEC) = SCR(KOFF)
        end do
        IOFF(1) = IOFF(1)+NNBSTR(ISYM,3)
      end do

      ! Update local vector counters.
      ! -----------------------------

      KVEC1 = KVEC1+LNUM
      LVEC1 = LVEC1+LNUM

    end if

  end do

  ! Update global vector counter.
  ! -----------------------------

  JVEC1 = JVEC1+NVRD

end do

end subroutine CHO_GETVEC2
