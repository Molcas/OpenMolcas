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

subroutine CHO_WRRSTC(IPASS)
!
! Purpose: write decomposition restart info for integral pass IPASS.
!
! NB!!!  The restart files are assumed open on entry.

use Cholesky, only: CHO_ADRVEC, Damp, InfRed, InfVec, IntMap, LuMap, LuRst, nBas, nnShl, nShell, nSym, NumCho, SCDIAG, Span, &
                    ThrCom, ThrDiag, ThrNeg, TOONEG, WARNEG
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IPASS
integer(kind=iwp), parameter :: LSCR = 10
integer(kind=iwp) :: IADR, IOPT, ISYM, J, JADR, JSCR(LSCR), NDIM, NTOT, NWR
real(kind=wp) :: DSCR(LSCR)

! Start address on file.
! ----------------------

IADR = 0

! Write molecular and configuration info.
! ---------------------------------------

IOPT = 1
NWR = 4
JSCR(1) = NSYM
JSCR(2) = NSHELL
JSCR(3) = NNSHL
JSCR(4) = 0
call IDAFILE(LURST,IOPT,JSCR,NWR,IADR)

IOPT = 1
NWR = NSYM
JSCR(1:NSYM) = NBAS(1:NSYM)
call IDAFILE(LURST,IOPT,JSCR,NWR,IADR)

IOPT = 1
NWR = 2
if (SCDIAG) then
  JSCR(1) = 1
else
  JSCR(1) = 0
end if
JSCR(2) = CHO_ADRVEC
call IDAFILE(LURST,IOPT,JSCR,NWR,IADR)

IOPT = 1
NWR = 8
DSCR(1) = THRCOM
DSCR(2) = THRDIAG
DSCR(3) = DAMP(1)
DSCR(4) = DAMP(2)
DSCR(5) = SPAN
DSCR(6) = THRNEG
DSCR(7) = WARNEG
DSCR(8) = TOONEG
call DDAFILE(LURST,IOPT,DSCR,NWR,IADR)

! Write vector info.
! ------------------

IOPT = 1
NWR = 1
JSCR(1) = IPASS
call IDAFILE(LURST,IOPT,JSCR,NWR,IADR)

IOPT = 1
call IDAFILE(LURST,IOPT,INFRED,IPASS,IADR)

do ISYM=1,NSYM
  IOPT = 1
  NWR = 1
  JSCR(1) = NUMCHO(ISYM)
  call IDAFILE(LURST,IOPT,JSCR,NWR,IADR)
  if (NUMCHO(ISYM) > 0) then
    do J=1,size(INFVEC,2)
      IOPT = 1
      NTOT = NUMCHO(ISYM)
      call IDAFILE(LURST,IOPT,InfVec(:,J,ISYM),NTOT,IADR)
    end do
  end if
end do

! Write integral shell pair map to disk.
! --------------------------------------

IOPT = 1
NDIM = 0
if (allocated(IntMap)) then
  NDIM = size(IntMap)
  JADR = 0
  call IDAFILE(LUMAP,IOPT,INTMAP,NDIM,JADR)
end if

end subroutine CHO_WRRSTC
