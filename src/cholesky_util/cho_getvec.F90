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

subroutine CHO_GETVEC(CHOVEC,LENVEC,NUMVEC,IVEC1,ISYM,SCR,LSCR)
!
! Purpose: read Cholesky vectors IVEC=IVEC1,....,IVEC1+NUMVEC-1
!          of symmetry ISYM from file. The vectors are returned
!          in the "current" reduced set. The algorithm used for
!          reading is taken from input (via Cholesky module).

use Cholesky, only: CHO_IOVEC, InfVec, iScr, LuCho, LuPri, MaxVec, nnBstR, nSym
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LENVEC, NUMVEC, IVEC1, ISYM, LSCR
real(kind=wp), intent(inout) :: CHOVEC(LENVEC,NUMVEC), SCR(LSCR)
integer(kind=iwp) :: IFAIL, IVEC, IVEC2, JADR, JVEC
real(kind=wp) :: XNRM
logical(kind=iwp), parameter :: LOCDBG = .false.
character(len=*), parameter :: SECNAM = 'CHO_GETVEC'
real(kind=wp), external :: ddot_

! Return if no vectors requested.
! -------------------------------

if (NUMVEC < 1) then
  if (LOCDBG) then
    write(LUPRI,*) SECNAM,': WARNING: no vectors in this call!'
    write(LUPRI,*) SECNAM,': NUMVEC = ',NUMVEC
  end if
  return
end if

! Check vector dimension: should be identical to current reduced
! set. Check also symmetry and vector index.
! --------------------------------------------------------------

if (LOCDBG) then
  IFAIL = 0
  if ((LENVEC /= NNBSTR(ISYM,2)) .or. (LENVEC < 1)) then
    write(LUPRI,*) SECNAM,': illegal vector dimension:'
    write(LUPRI,*) SECNAM,': LENVEC = ',LENVEC
    IFAIL = IFAIL+1
  end if
  if ((ISYM < 1) .or. (ISYM > NSYM)) then
    write(LUPRI,*) SECNAM,': illegal symmetry input:'
    write(LUPRI,*) SECNAM,': ISYM = ',ISYM
    IFAIL = IFAIL+1
  end if
  IVEC2 = IVEC1+NUMVEC-1
  if ((IVEC1 < 1) .or. (IVEC1 > MAXVEC) .or. (IVEC2 < 1) .or. (IVEC2 > MAXVEC)) then
    write(LUPRI,*) SECNAM,': illegal vector index:'
    write(LUPRI,*) SECNAM,': IVEC1,IVEC2 = ',IVEC1,IVEC2
    IFAIL = IFAIL+1
  else
    if (INFVEC(IVEC1,3,ISYM) < 0) then
      write(LUPRI,*) SECNAM,': illegal first vector address:'
      write(LUPRI,*) SECNAM,': address: ',INFVEC(IVEC1,3,ISYM)
      IFAIL = IFAIL+1
    end if
    if (INFVEC(IVEC2,3,ISYM) < 0) then
      write(LUPRI,*) SECNAM,': illegal last vector address:'
      write(LUPRI,*) SECNAM,': address: ',INFVEC(IVEC2,3,ISYM)
      IFAIL = IFAIL+1
    end if
  end if
  if ((CHO_IOVEC == 1) .or. (CHO_IOVEC == 2) .or. (CHO_IOVEC == 3) .or. (CHO_IOVEC == 4)) then
    if (size(ISCR) < NNBSTR(ISYM,2)) then
      write(LUPRI,*) SECNAM,': insufficient iscratch:'
      write(LUPRI,*) SECNAM,': SIZE(ISCR) = ',size(ISCR)
      IFAIL = IFAIL+1
    end if
  end if
  if (IFAIL /= 0) then
    write(LUPRI,*) SECNAM,': unable to continue!'
    call CHO_QUIT('Error in '//SECNAM,104)
  end if
end if

! Call reading routine.
! ---------------------

if (CHO_IOVEC == 1) then
  call CHO_GETVEC1(CHOVEC,LENVEC,NUMVEC,IVEC1,ISYM,SCR,LSCR)
else if ((CHO_IOVEC == 2) .or. (CHO_IOVEC == 3) .or. (CHO_IOVEC == 4)) then
  call CHO_GETVEC2(CHOVEC,LENVEC,NUMVEC,IVEC1,ISYM,SCR,LSCR)
else
  call CHO_GETVEC0(CHOVEC,LENVEC,NUMVEC,IVEC1,ISYM,SCR,LSCR)
end if

! Debug print.
! ------------

if (LOCDBG) then
  call XFLUSH(LUPRI)
  write(LUPRI,*)
  write(LUPRI,*) SECNAM,':'
  write(LUPRI,*) 'Vectors ',IVEC1,' to ',IVEC1+NUMVEC-1,' of symmetry ',ISYM,' read from unit ',LUCHO(ISYM)
  write(LUPRI,*) 'Vector dimension: ',LENVEC,' (current reduced set)'
  write(LUPRI,*) 'Algorithm: ',CHO_IOVEC
  do IVEC=1,NUMVEC
    JVEC = IVEC1+IVEC-1
    JADR = INFVEC(JVEC,3,ISYM)
    XNRM = sqrt(DDOT_(LENVEC,CHOVEC(:,IVEC),1,CHOVEC(:,IVEC),1))
    write(LUPRI,*) 'Vector:',JVEC,' address: ',JADR,' norm: ',XNRM
  end do
  call XFLUSH(LUPRI)
end if

end subroutine CHO_GETVEC
