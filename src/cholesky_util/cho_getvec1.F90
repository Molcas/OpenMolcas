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

subroutine CHO_GETVEC1(CHOVEC,LENVEC,NUMVEC,IVEC1,ISYM,SCR,LSCR)
!
! Purpose: read Cholesky vectors IVEC=IVEC1,....,IVEC1+NUMVEC-1
!          of symmetry ISYM from file. The vectors are returned
!          in the "current" reduced set. This routine attempts
!          to minimize gather/scatter operations and uses batched
!          reading to (hopefully) improve buffering.
!
! NOTE: the scratch array SCR(LSCR) is used to read vectors from
!       disk and should not be smaller than NNBSTR(ISYM,1)+1.

use Cholesky, only: Cho_AdrVec, InfVec, iScr, LuPri, nnBstR, nSys_call
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LENVEC, NUMVEC, IVEC1, ISYM, LSCR
real(kind=wp), intent(inout) :: CHOVEC(LENVEC,NUMVEC), SCR(LSCR)
integer(kind=iwp) :: IAB, IBATCH, IBVEC1, ILOC, IOFF(0:1), IRED, IRED1, IRED2, IREDC, IVEC2, JNUM, JNUM_RD, JRED, JVEC, JVEC1, &
                     JVEC2, JVEC_END, KBVEC1, KJUNK, KOFF, KSCR, KVEC, KVEC1, LEFT, LTOT, MINL, MUSED, NBATCH, NUMV, NVEC
logical(kind=iwp), parameter :: LOCDBG = .false.
character(len=*), parameter :: SECNAM = 'CHO_GETVEC1'

! Some initializations.
! ---------------------

ILOC = 3
IVEC2 = IVEC1+NUMVEC-1

KJUNK = 1
KSCR = KJUNK+1
LEFT = LSCR-KSCR+1
if (LEFT <= 0) call CHO_QUIT('Insufficient scratch space in '//SECNAM,101)

SCR(KJUNK) = Zero
IOFF(0) = KJUNK

! Get reduced sets of first and last vector.
! ------------------------------------------

IRED1 = INFVEC(IVEC1,2,ISYM)
IRED2 = INFVEC(IVEC2,2,ISYM)

! Loop through reduced sets to be read.
! -------------------------------------

KVEC1 = 1
JVEC1 = IVEC1
do IRED=IRED1,IRED2

  ! Count vectors in this reduced set.
  ! ----------------------------------

  JNUM = 0
  JVEC = JVEC1-1
  do while (JVEC < IVEC2)
    JVEC = JVEC+1
    JRED = INFVEC(JVEC,2,ISYM)
    if (JRED == IRED) then
      JNUM = JNUM+1 ! increase counter
    else
      JVEC = IVEC2 ! break while loop
    end if
  end do

  ! Skip if this reduced set is empty.
  ! ----------------------------------

  if (JNUM /= 0) then

    ! Check vector range.
    ! -------------------

    if (LOCDBG) then
      JVEC2 = JVEC1+JNUM-1
      if (JVEC2 > IVEC2) then
        write(LUPRI,*) SECNAM,': IRED  = ',IRED
        write(LUPRI,*) SECNAM,': JNUM  = ',JNUM
        write(LUPRI,*) SECNAM,': JVEC1 = ',JVEC1
        write(LUPRI,*) SECNAM,': JVEC2 = ',JVEC2
        call CHO_QUIT('Vector index error in '//SECNAM,103)
      end if
    end if

    ! Read reduced set index arrays.
    ! ------------------------------

    call CHO_GETRED(IRED,ILOC,.false.)
    call CHO_SETREDIND(ILOC)

    ! If reduced sets are identical, simply read the vectors
    ! directly into CHOVEC array and go to next reduced set.
    ! ------------------------------------------------------

    if (NNBSTR(ISYM,3) == NNBSTR(ISYM,2)) then
      !if (CHO_ADRVEC == 1) then
      !  IOPT = 2
      !  IADR = INFVEC(JVEC1,3,ISYM)
      !  LTOT = NNBSTR(ISYM,2)*JNUM
      !  call DDAFILE(LUCHO(ISYM),IOPT,CHOVEC(1,KVEC1),LTOT,IADR)
      !  NSYS_CALL = NSYS_CALL+1
      !else if (CHO_ADRVEC == 2) then
      !  IOPT = 2
      !  LTOT = NNBSTR(ISYM,2)
      !  do KK=1,JNUM
      !    IADR = INFVEC(JVEC1+KK-1,3,ISYM)
      !    call DDAFILE(LUCHO(ISYM),IOPT,CHOVEC(1,KVEC1+KK-1),LTOT,IADR)
      !    NSYS_CALL = NSYS_CALL+1
      !  end do
      !else
      !  call CHO_QUIT('[1] CHO_ADRVEC error in '//SECNAM,102)
      !end if
      !-tbp: replaced above to make use of buffer in cho_vecrd.
      LTOT = NNBSTR(ISYM,2)*JNUM
      JVEC_END = JVEC1+JNUM-1
      JNUM_RD = 0
      IREDC = IRED
      MUSED = 0
      call CHO_VECRD(CHOVEC(1,KVEC1),LTOT,JVEC1,JVEC_END,ISYM,JNUM_RD,IREDC,MUSED)
      if (JNUM_RD /= JNUM) call CHO_QUIT('Logical error [RD1] in '//SECNAM,103)
      NSYS_CALL = NSYS_CALL+1
    else

      ! Set up batch over vectors to be read.
      ! -------------------------------------

      MINL = NNBSTR(ISYM,3)
      if (MINL < 1) then
        NVEC = 0
      else
        NVEC = min(LEFT/MINL,JNUM)
      end if
      if (NVEC < 1) then
        write(LUPRI,*) SECNAM,': insufficient scratch space:'
        write(LUPRI,*) 'LEFT = ',LEFT
        write(LUPRI,*) 'JNUM = ',JNUM
        write(LUPRI,*) 'MINL = ',MINL
        write(LUPRI,*) 'NVEC = ',NVEC
        write(LUPRI,*) 'Input:'
        write(LUPRI,*) 'IVEC1  = ',IVEC1
        write(LUPRI,*) 'NUMVEC = ',NUMVEC
        write(LUPRI,*) 'LENVEC = ',LENVEC
        write(LUPRI,*) 'ISYM   = ',ISYM
        call CHO_QUIT('Insufficient scratch space in '//SECNAM,104)
        NBATCH = 0 ! to avoid compiler warnings
      else
        NBATCH = (JNUM-1)/NVEC+1
      end if

      ! Set up mapping between reduced sets.
      ! ------------------------------------

      call CHO_RS2RS(ISCR,size(ISCR),2,3,IRED,ISYM)

      ! Start batch loop.
      ! -----------------

      do IBATCH=1,NBATCH

        if (IBATCH == NBATCH) then
          NUMV = JNUM-NVEC*(NBATCH-1)
        else
          NUMV = NVEC
        end if
        IBVEC1 = JVEC1+NVEC*(IBATCH-1)
        KBVEC1 = KVEC1+NVEC*(IBATCH-1)

        ! Read vectors.
        ! -------------

        !if (CHO_ADRVEC == 1) then
        !   IOPT = 2
        !   LTOT = NNBSTR(ISYM,3)*NUMV
        !   IADR = INFVEC(IBVEC1,3,ISYM)
        !   call DDAFILE(LUCHO(ISYM),IOPT,SCR(KSCR),LTOT,IADR)
        !   NSYS_CALL = NSYS_CALL+1
        !else if (CHO_ADRVEC == 2) then
        !   IOPT = 2
        !   LTOT = NNBSTR(ISYM,3)
        !   KTRG = KSCR
        !   do KK=1,NUMV
        !      IADR = INFVEC(IBVEC1+KK-1,3,ISYM)
        !      call DDAFILE(LUCHO(ISYM),IOPT,SCR(KTRG),LTOT,IADR)
        !      KTRG = KTRG+NNBSTR(ISYM,3)
        !      NSYS_CALL = NSYS_CALL+1
        !   end do
        !else
        !   call CHO_QUIT('[2] CHO_ADRVEC error in '//SECNAM,102)
        !end if
        !-tbp: replaced above to make use of buffer in cho_vecrd.
        LTOT = NNBSTR(ISYM,3)*NUMV
        JVEC_END = IBVEC1+NUMV-1
        JNUM_RD = 0
        IREDC = IRED
        MUSED = 0
        call CHO_VECRD(SCR(KSCR),LTOT,IBVEC1,JVEC_END,ISYM,JNUM_RD,IREDC,MUSED)
        if (JNUM_RD /= NUMV) call CHO_QUIT('Logical error [RD2] in '//SECNAM,103)
        if (CHO_ADRVEC == 1) then
          NSYS_CALL = NSYS_CALL+1
        else if (CHO_ADRVEC == 2) then
          NSYS_CALL = NSYS_CALL+NUMV
        else
          call CHO_QUIT('CHO_ADRVEC error in '//SECNAM,102)
        end if

        ! Copy vectors into result array.
        ! -------------------------------

        do JVEC=1,NUMV
          KVEC = KBVEC1+JVEC-1
          IOFF(1) = IOFF(0)+NNBSTR(ISYM,3)*(JVEC-1)
          do IAB=1,NNBSTR(ISYM,2)
            KOFF = IOFF(min(ISCR(IAB),1))+ISCR(IAB)
            CHOVEC(IAB,KVEC) = SCR(KOFF)
          end do
        end do

      end do
    end if
  end if

  ! Set next vector to be treated.
  ! ------------------------------

  KVEC1 = KVEC1+JNUM
  JVEC1 = JVEC1+JNUM

end do

end subroutine CHO_GETVEC1
