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

subroutine CHO_PUTVEC(CHOVEC,LENVEC,NUMVEC,IVEC1,ISYM)
!
! Purpose: write Cholesky vectors IVEC=IVEC1,...,IVEC1+NUMVEC-1
!          of symmetry ISYM to file.

use Cholesky, only: Cho_AdrVec, Cho_Real_Par, InfVec, LuCho, LuPri, MaxVec, nnBstR, nSym
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: LENVEC, NUMVEC, IVEC1, ISYM
real(kind=wp), intent(_IN_) :: CHOVEC(LENVEC,NUMVEC)
integer(kind=iwp) :: IADR, IADR2, IOPT, IVEC, IVEC2, JADR, JVEC, LTOT
real(kind=wp) :: XNRM
logical(kind=iwp) :: CHK_OVERFLOW
logical(kind=iwp), parameter :: LOCDBG = .false.
character(len=*), parameter :: SECNAM = 'CHO_PUTVEC'
real(kind=wp), external :: ddot_

CHK_OVERFLOW = .not. Cho_Real_Par

! Return if no vectors.
! ---------------------

if (NUMVEC < 1) then
  if (LOCDBG) then
    write(LUPRI,*) SECNAM,': WARNING: no vectors in this call!'
    write(LUPRI,*) SECNAM,': NUMVEC = ',NUMVEC
  end if
  return
end if

! Check vector dimension: should be the same as current reduced set.
! ------------------------------------------------------------------

if (LENVEC /= NNBSTR(ISYM,2)) then
  call CHO_QUIT('Illegal vector dimension in '//SECNAM,104)
end if
if (LENVEC < 1) then
  if (LOCDBG) then
    write(LUPRI,*) SECNAM,': WARNING: negative vector dimension'
    write(LUPRI,*) SECNAM,': LENVEC = ',LENVEC
  end if
  return
end if

! Check symmetry.
! ---------------

if ((ISYM < 1) .or. (ISYM > NSYM)) then
  write(LUPRI,*) SECNAM,': symmetry out of bounds'
  write(LUPRI,*) 'ISYM = ',ISYM
  call CHO_QUIT('Symmetry out of bounds in '//SECNAM,104)
end if

! Check vector index.
! -------------------

IVEC2 = IVEC1+NUMVEC-1
if ((IVEC1 < 1) .or. (IVEC1 > MAXVEC) .or. (IVEC2 < 1) .or. (IVEC2 > MAXVEC)) then
  write(LUPRI,*) SECNAM,': vector index out of bounds'
  write(LUPRI,*) 'IVEC1 = ',IVEC1,' IVEC2 = ',IVEC2
  write(LUPRI,*) '...must be between 1 and ',MAXVEC
  call CHO_QUIT('Vector index out of bounds in '//SECNAM,104)
end if

! Check for overflow for WA file addressing.
! ------------------------------------------

if (CHK_OVERFLOW .and. (CHO_ADRVEC == 1)) then
  IADR2 = INFVEC(IVEC2,4,ISYM)
  if (INFVEC(IVEC1,4,ISYM) < 0) then
    write(LUPRI,*) 'Error in ',SECNAM,':'
    write(LUPRI,*) 'Illegal disk address for first vector: ',INFVEC(IVEC1,4,ISYM)
    if (INFVEC(IVEC1,4,ISYM) < -1) write(LUPRI,*) '....is it an overflow?'
    write(LUPRI,*) 'IVEC1 = ',IVEC1,' ISYM = ',ISYM
    call CHO_QUIT('Illegal disk address in '//SECNAM,104)
  else if (IADR2 < INFVEC(IVEC1,4,ISYM)) then
    write(LUPRI,*) 'Error in ',SECNAM,':'
    write(LUPRI,*) 'Illegal disk address for last vector: ',IADR2
    if (IADR2 < -1) write(LUPRI,*) '....is it an overflow?'
    write(LUPRI,*) 'IVEC2 = ',IVEC2,' ISYM = ',ISYM
    call CHO_QUIT('Illegal disk address in '//SECNAM,104)
  end if
end if

! Call the low-level I/O routines.
! CHO_ADRVEC=1: WA files.
! CHO_ADRVEC=2: DA files.
! Set (next) disk addresses.
! --------------------------------

if (CHO_ADRVEC == 1) then
  IOPT = 1  ! synchronous write option
  LTOT = LENVEC*NUMVEC
  IADR = INFVEC(IVEC1,3,ISYM)
  call DDAFILE(LUCHO(ISYM),IOPT,CHOVEC,LTOT,IADR)
  do IVEC=1,NUMVEC-1
    JVEC = IVEC1+IVEC-1
    INFVEC(JVEC+1,3,ISYM) = INFVEC(JVEC,3,ISYM)+LENVEC
  end do
  IVEC = NUMVEC
  JVEC = IVEC1+IVEC-1
  IADR = INFVEC(JVEC,3,ISYM)
  if (JVEC < MAXVEC) INFVEC(JVEC+1,3,ISYM) = INFVEC(JVEC,3,ISYM)+LENVEC
else if (CHO_ADRVEC == 2) then
  IOPT = 1  ! synchronous write option
  LTOT = LENVEC
  do IVEC=1,NUMVEC-1
    JVEC = IVEC1+IVEC-1
    IADR = INFVEC(JVEC,3,ISYM)
    call DDAFILE(LUCHO(ISYM),IOPT,CHOVEC(1,IVEC),LTOT,IADR)
    INFVEC(JVEC+1,3,ISYM) = IADR
  end do
  IVEC = NUMVEC
  JVEC = IVEC1+IVEC-1
  IADR = INFVEC(JVEC,3,ISYM)
  call DDAFILE(LUCHO(ISYM),IOPT,CHOVEC(1,IVEC),LTOT,IADR)
  if (JVEC < MAXVEC) INFVEC(JVEC+1,3,ISYM) = IADR
else
  call CHO_QUIT('CHO_ADRVEC out of bounds in '//SECNAM,102)
end if

! Debug stuff.
! ------------

if (LOCDBG) then
  write(LUPRI,*)
  write(LUPRI,*) SECNAM,':'
  write(LUPRI,*) 'Vectors ',IVEC1,' to ',IVEC1+NUMVEC-1,' of symmetry ',ISYM,' written to unit ',LUCHO(ISYM)
  write(LUPRI,*) 'Vector dimension: ',LENVEC
  do IVEC=1,NUMVEC
    JVEC = IVEC1+IVEC-1
    JADR = INFVEC(JVEC,3,ISYM)
    XNRM = sqrt(DDOT_(LENVEC,CHOVEC(1,IVEC),1,CHOVEC(1,IVEC),1))
    write(LUPRI,*) 'Vector:',JVEC,' address: ',JADR,' norm: ',XNRM
  end do
end if

end subroutine CHO_PUTVEC
