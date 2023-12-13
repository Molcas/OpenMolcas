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

subroutine CHO_VECRD1(SCR,LSCR,JVEC1,IVEC2,ISYM,JNUM,IREDC,MUSED,DOREAD)
!
! Purpose: read as many vectors as fit into SCR array starting
!          at vector JVEC1 and reading at most until vector IVEC2.
!          On exit, JNUM is the number of vectors read.
!          On entry as well as exit, IREDC identifies the reduced
!          set stored in core (at position "3"; use -1 if none
!          or unkown). If DOREAD=.false. no vectors are actually
!          read in, but JNUM and MUSED are returned as appropriate.
!          Thus, array SCR is not referenced for DOREAD=.false.
!
! NOTE: if no vectors can be read, JNUM=0 and MUSED=0 are returned,
!       but execution is NOT stopped here!!!

use Cholesky, only: Cho_AdrVec, InfVec, LuCho, LuPri, nDimRS, nnBstR
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LSCR, JVEC1, IVEC2, ISYM
real(kind=wp), intent(inout) :: SCR(LSCR)
integer(kind=iwp), intent(out) :: JNUM, MUSED
integer(kind=iwp), intent(inout) :: IREDC
logical(kind=iwp), intent(in) :: DOREAD
integer(kind=iwp) :: IADR, ILOC, IOPT, IVEC, JADR, JRED, JVEC, KOFFV, KSCR, LENR, LTOT, NTST
real(kind=wp) :: XNRM
logical(kind=iwp) :: FULL
logical(kind=iwp), parameter :: LOCDBG = .false.
character(len=*), parameter :: SECNAM = 'CHO_VECRD1'
real(kind=wp), external :: ddot_

JRED = 0  ! fix compiler warning

if (CHO_ADRVEC == 1) then  ! WA addressing

  ! Count how many vectors can be read.
  ! -----------------------------------

  JNUM = 0
  LTOT = 0
  JVEC = JVEC1-1
  FULL = LTOT >= LSCR
  if (.not. allocated(NDIMRS)) then
    ILOC = 3 ! use scratch location in reduced index arrays
    do while ((JVEC < IVEC2) .and. (.not. FULL))
      JVEC = JVEC+1
      JRED = INFVEC(JVEC,2,ISYM)
      if (JRED /= IREDC) then
        call CHO_GETRED(JRED,ILOC,.false.)
        call CHO_SETREDIND(ILOC)
        IREDC = JRED
      end if
      LTOT = LTOT+NNBSTR(ISYM,ILOC)
      if (LTOT > LSCR) then
        JVEC = JVEC-1
        LTOT = LTOT-NNBSTR(ISYM,ILOC)
        FULL = .true.
      else
        JNUM = JNUM+1
      end if
    end do
  else
    do while ((JVEC < IVEC2) .and. (.not. FULL))
      JVEC = JVEC+1
      JRED = INFVEC(JVEC,2,ISYM)
      LTOT = LTOT+NDIMRS(ISYM,JRED)
      if (LTOT > LSCR) then
        JVEC = JVEC-1
        LTOT = LTOT-NDIMRS(ISYM,JRED)
        FULL = .true.
      else
        JNUM = JNUM+1
      end if
    end do
  end if

  ! Read vectors (if any).
  ! ----------------------

  if (DOREAD .and. (LTOT > 0)) then
    IOPT = 2
    IADR = INFVEC(JVEC1,3,ISYM)
    call DDAFILE(LUCHO(ISYM),IOPT,SCR,LTOT,IADR)
  end if

else if (CHO_ADRVEC == 2) then ! DA adressing

  ! Read as many vectors as can be read, one at a time.
  ! ---------------------------------------------------

  JNUM = 0
  LTOT = 0
  KSCR = 1
  JVEC = JVEC1-1
  FULL = LTOT >= LSCR
  if (.not. allocated(NDIMRS)) then
    ILOC = 3 ! use scratch location in reduced index arrays
    do while ((JVEC < IVEC2) .and. (.not. FULL))
      JVEC = JVEC+1
      JRED = INFVEC(JVEC,2,ISYM)
      if (JRED /= IREDC) then
        call CHO_GETRED(JRED,ILOC,.false.)
        call CHO_SETREDIND(ILOC)
        IREDC = JRED
      end if
      LTOT = LTOT+NNBSTR(ISYM,ILOC)
      if (LTOT > LSCR) then
        JVEC = JVEC-1
        LTOT = LTOT-NNBSTR(ISYM,ILOC)
        FULL = .true.
      else
        JNUM = JNUM+1
        if (DOREAD) then
          IOPT = 2
          LENR = NNBSTR(ISYM,ILOC)
          IADR = INFVEC(JVEC,3,ISYM)
          call DDAFILE(LUCHO(ISYM),IOPT,SCR(KSCR),LENR,IADR)
          KSCR = KSCR+NNBSTR(ISYM,ILOC)
        end if
      end if
    end do
  else
    do while ((JVEC < IVEC2) .and. (.not. FULL))
      JVEC = JVEC+1
      JRED = INFVEC(JVEC,2,ISYM)
      LTOT = LTOT+NDIMRS(ISYM,JRED)
      if (LTOT > LSCR) then
        JVEC = JVEC-1
        LTOT = LTOT-NDIMRS(ISYM,JRED)
        FULL = .true.
      else
        JNUM = JNUM+1
        if (DOREAD) then
          IOPT = 2
          LENR = NDIMRS(ISYM,JRED)
          IADR = INFVEC(JVEC,3,ISYM)
          call DDAFILE(LUCHO(ISYM),IOPT,SCR(KSCR),LENR,IADR)
          KSCR = KSCR+NDIMRS(ISYM,JRED)
        end if
      end if
    end do
  end if

else  ! unknown file addressing

  call CHO_QUIT('CHO_ADRVEC error in '//SECNAM,102)
  LTOT = 0 ! dummy assignment to avoid compiler warnings

end if

! Return total memory used for read.
! ----------------------------------

MUSED = LTOT

! Debug print.
! ------------

if (LOCDBG) then
  write(LUPRI,*)
  write(LUPRI,*) SECNAM,':'
  write(LUPRI,*) 'Vector addressing: ',CHO_ADRVEC
  write(LUPRI,*) 'DOREAD: ',DOREAD
  if (JNUM < 1) then
    if (DOREAD) then
      write(LUPRI,*) 'No vectors read!'
    else
      write(LUPRI,*) 'No vectors can be read!'
    end if
  else
    if (DOREAD) then
      write(LUPRI,*) 'Vectors ',JVEC1,' to ',JVEC1+JNUM-1,' of symmetry ',ISYM,' read from unit ',LUCHO(ISYM)
      if (allocated(NDIMRS)) then
        KOFFV = 1
        do IVEC=1,JNUM
          JVEC = JVEC1+IVEC-1
          JADR = INFVEC(JVEC,3,ISYM)
          JRED = INFVEC(JVEC,2,ISYM)
          XNRM = sqrt(DDOT_(NDIMRS(ISYM,JRED),SCR(KOFFV),1,SCR(KOFFV),1))
          write(LUPRI,*) 'Vector:',JVEC,' address: ',JADR,' norm: ',XNRM
          KOFFV = KOFFV+NDIMRS(ISYM,JRED)
        end do
        NTST = KOFFV-1
        if (NTST /= MUSED) call CHO_QUIT('Vector dimension error in '//SECNAM,104)
      end if
    else
      write(LUPRI,*) 'Vectors ',JVEC1,' to ',JVEC1+JNUM-1,' of symmetry ',ISYM,' can be read'
    end if
  end if
  call XFLUSH(LUPRI)
end if

end subroutine CHO_VECRD1
