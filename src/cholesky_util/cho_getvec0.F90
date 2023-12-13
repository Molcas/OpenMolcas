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

subroutine CHO_GETVEC0(CHOVEC,LENVEC,NUMVEC,IVEC1,ISYM,SCR,LSCR)
!=======================================================================
!==== DEPRECATED - USE CHO_X_GETVEC OR CHO_X_VECRD INSTEAD =============
!=======================================================================
!
! Purpose: read Cholesky vectors IVEC=IVEC1,...,IVEC1+NUMVEC-1
!          of symmetry ISYM from file. The vectors are returned
!          in the "current" reduced set.
!
! NOTE: array SCR(LSCR) is used for storing the vectors in the
!       red. set from disk and for a full first red. set vector.
!       Thus, to be certain that enough memory is available,
!       use LSCR = 2 x dimension of first reduced set.

use Cholesky, only: iiBstR, IndRed, InfVec, LuPri, nnBstR, nSys_call
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LENVEC, NUMVEC, IVEC1, ISYM, LSCR
real(kind=wp), intent(out) :: CHOVEC(LENVEC,NUMVEC)
real(kind=wp), intent(inout) :: SCR(LSCR)
integer(kind=iwp) :: IAB, ILOC, IRED, IREDC, IVEC, JAB, JNUM, JRED, JVEC, KAB, KEND1, KREAD, KRED1, KREDU, LSCR1, MUSED
real(kind=wp) :: XNRM
logical(kind=iwp), parameter :: LOCDBG = .false.
character(len=*), parameter :: SECNAM = 'CHO_GETVEC0'
real(kind=wp), external :: ddot_

! Initialize output array.
! ------------------------

CHOVEC(:,:) = Zero

! Read reduced set index arrays for first vector.
! -----------------------------------------------

IRED = INFVEC(IVEC1,2,ISYM)
ILOC = 3
call CHO_GETRED(IRED,ILOC,.false.)
call CHO_SETREDIND(3)
KRED1 = 1
KREAD = KRED1+NNBSTR(ISYM,1)
KEND1 = KREAD+NNBSTR(ISYM,3)
LSCR1 = LSCR-KEND1+1
if (LSCR1 < 0) then
  write(LUPRI,*) 'Insufficient scratch space in ',SECNAM
  write(LUPRI,*) 'Available: ',LSCR,'   Need: ',KEND1-1
  write(LUPRI,*) '- needed for RED1: ',NNBSTR(ISYM,1)
  write(LUPRI,*) '- needed for READ: ',NNBSTR(ISYM,3)
  call CHO_QUIT('[1] Insufficient scratch space in '//SECNAM,102)
end if

! Read vectors and re-order into current reduced set via reduced set 1.
! NOTE: if the read vectors are already in red. set 1, don't resort.
! ---------------------------------------------------------------------

do JVEC=1,NUMVEC

  IVEC = IVEC1+JVEC-1
  JRED = INFVEC(IVEC,2,ISYM)
  if (JRED /= IRED) then   ! read new reduced set
    call CHO_GETRED(JRED,ILOC,.false.)
    call CHO_SETREDIND(3)
    KEND1 = KREAD+NNBSTR(ISYM,3)
    LSCR1 = LSCR-KEND1+1
    if (LSCR1 < 0) then
      write(LUPRI,*) 'Insufficient scratch space in ',SECNAM
      write(LUPRI,*) 'Available: ',LSCR,'   Need: ',KEND1-1
      write(LUPRI,*) '- needed for RED1: ',NNBSTR(ISYM,1)
      write(LUPRI,*) '- needed for READ: ',NNBSTR(ISYM,3)
      call CHO_QUIT('[2] Insufficient scratch space in '//SECNAM,102)
    end if
    IRED = JRED
  end if

  !IOPT = 2
  !IADR = INFVEC(IVEC,3,ISYM)
  !CALL DDAFILE(LUCHO(ISYM),IOPT,SCR(KREAD),NNBSTR(ISYM,3),IADR)
  !-tbp: replaced above code to make use of buffer through cho_vecrd.
  JNUM = 0
  IREDC = IRED
  MUSED = 0
  call CHO_VECRD(SCR(KREAD),NNBSTR(ISYM,3),IVEC,IVEC,ISYM,JNUM,IREDC,MUSED)
  if (JNUM /= 1) call CHO_QUIT('Logical error in '//SECNAM,103)
  NSYS_CALL = NSYS_CALL+1
  if (LOCDBG) then
    XNRM = sqrt(DDOT_(NNBSTR(ISYM,3),SCR(KREAD),1,SCR(KREAD),1))
    write(LUPRI,*) SECNAM,': Vector:',IVEC,' address: ',INFVEC(IVEC,3,ISYM),' norm: ',XNRM,' sym. ',ISYM,' red. set: ',IRED, &
                   ' dim.: ',NNBSTR(ISYM,3)
  end if

  SCR(KRED1:KRED1+NNBSTR(ISYM,1)-1) = Zero
  if (IRED > 1) then
    do JAB=1,NNBSTR(ISYM,3)   ! sort into rs1 ordering
      KAB = IIBSTR(ISYM,3)+JAB
      IAB = INDRED(KAB,3)-IIBSTR(ISYM,1)
      SCR(KRED1+IAB-1) = SCR(KREAD+JAB-1)
    end do
    KREDU = KRED1   ! point rs2 sort to red1 resort
  else if (IRED == 1) then
    KREDU = KREAD   ! point rs2 sort to read (already red1)
  else
    write(LUPRI,*) SECNAM,': ERROR: IRED is negative: ',IRED
    call CHO_QUIT('Reduced set error in '//SECNAM,104)
    KREDU = -999999 ! just to avoid compiler warnings
  end if

  do JAB=1,NNBSTR(ISYM,2)  ! sort into in rs2 ordering
    KAB = IIBSTR(ISYM,2)+JAB
    IAB = INDRED(KAB,2)-IIBSTR(ISYM,1)
    CHOVEC(JAB,JVEC) = SCR(KREDU+IAB-1)
  end do

end do

end subroutine CHO_GETVEC0
