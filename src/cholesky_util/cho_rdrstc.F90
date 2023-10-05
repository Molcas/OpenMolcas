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

subroutine CHO_RDRSTC(IFAIL)
!
! Purpose: read decomposition restart info and store in common
!          block. If IFAIL != 0 on exit, some error occurred and,
!          most likely, some of the restart info is not
!          defined/initialized.
!
! NB!!!! the restart files MUST be open on entry...

use Cholesky, only: InfRed, InfVec, LuPri, LuRst, MaxRed, MaxVec, nSym, NumCho, XCho_AdrVec, XDamp, XNBAS, XnnShl, XnPass, &
                    XnShell, XnSym, XScDiag, XSpan, XThrCom, XThrDiag, XThrNeg, XTooNeg, XWarNeg
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: IFAIL
integer(kind=iwp), parameter :: LSCR = 8
integer(kind=iwp) :: IADR, IOPT, ISYM, J, JSCR(LSCR), NRD
real(kind=wp) :: DSCR(LSCR)
character(len=*), parameter :: SECNAM = 'CHO_RDRSTC'

! Set return code.
! ----------------

IFAIL = 0

! Read molecular info.
! --------------------

IADR = 0

IOPT = 2
NRD = 4
call IDAFILE(LURST,IOPT,JSCR,NRD,IADR)
XNSYM = JSCR(1)
XNSHELL = JSCR(2)
XNNSHL = JSCR(3)
if ((XNSYM < 1) .or. (XNSYM > 8)) then
  write(LUPRI,'(A,A,I10)') SECNAM,': #irreps from restart file: ',XNSYM
  IFAIL = 1
  call Finish_this()
  return
else
  IOPT = 2
  call IDAFILE(LURST,IOPT,XNBAS,XNSYM,IADR)
end if

! Read decomposition configuration info.
! --------------------------------------

IOPT = 2
NRD = 2
call IDAFILE(LURST,IOPT,JSCR,NRD,IADR)
if (JSCR(1) == 0) then
  XSCDIAG = .false.
else if (JSCR(1) == 1) then
  XSCDIAG = .true.
else
  write(LUPRI,'(A,A,I10)') SECNAM,': integer flag for screening not recognized:',JSCR(1)
  IFAIL = 2
  call Finish_this()
  return
end if
XCHO_ADRVEC = JSCR(2)

IOPT = 2
NRD = 8
call DDAFILE(LURST,IOPT,DSCR,NRD,IADR)
XTHRCOM = DSCR(1)
XTHRDIAG = DSCR(2)
XDAMP(1) = DSCR(3)
XDAMP(2) = DSCR(4)
XSPAN = DSCR(5)
XTHRNEG = DSCR(6)
XWARNEG = DSCR(7)
XTOONEG = DSCR(8)

! Read decomposition info.
! ------------------------

IOPT = 2
NRD = 1
call IDAFILE(LURST,IOPT,JSCR,NRD,IADR)
XNPASS = JSCR(1)
if ((XNPASS < 1) .or. (XNPASS > MAXRED)) then
  write(LUPRI,'(A,A,I10)') SECNAM,': #reduced sets in restart:',XNPASS
  IFAIL = 3
  call Finish_this()
  return
else
  IOPT = 2
  INFRED(:) = 0
  call IDAFILE(LURST,IOPT,INFRED,XNPASS,IADR)
  if (INFRED(1) /= 0) then
    write(LUPRI,'(A,A,I10)') SECNAM,': disk address of 1st reduced set:',INFRED(1)
    IFAIL = 4
    call Finish_this()
    return
  end if
end if

do ISYM=1,NSYM
  IOPT = 2
  NRD = 1
  call IDAFILE(LURST,IOPT,JSCR,NRD,IADR)
  NUMCHO(ISYM) = JSCR(1)
  if ((NUMCHO(ISYM) < 0) .or. (NUMCHO(ISYM) > MAXVEC)) then
    write(LUPRI,'(A,A,I2,A,I10)') SECNAM,': #Cholesky vectors (sym.',ISYM,'): ',NUMCHO(ISYM)
    IFAIL = 5
    call Finish_this()
    return
  else if (NUMCHO(ISYM) == 0) then
    INFVEC(:,:,ISYM) = 0
  else
    INFVEC(:,:,ISYM) = 0
    do J=1,size(INFVEC,2)
      IOPT = 2
      call IDAFILE(LURST,IOPT,INFVEC(:,J,ISYM),NUMCHO(ISYM),IADR)
    end do
  end if
end do

call Finish_this()

contains

! failures jump to this point
subroutine Finish_this()

  if (IFAIL /= 0) write(LUPRI,'(A,A)') SECNAM,': refusing to read more restart info!'

end subroutine Finish_this

end subroutine CHO_RDRSTC
