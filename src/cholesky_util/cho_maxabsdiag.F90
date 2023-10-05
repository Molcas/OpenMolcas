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

subroutine CHO_MAXABSDIAG(DIAG,IRED,DGMAX)
!
! Purpose: set max. abs. DIAG (reduced set IRED) in each symmetry, and
!          return the global max. abs. in DGMAX.

use Cholesky, only: Cho_1Center, DiaMax, DiaMaxT, iiBstR, IndRed, LuPri, nnBstR, nSym
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Diag(*)
integer(kind=iwp), intent(in) :: IRED
real(kind=wp), intent(out) :: DGMax
integer(kind=iwp) :: AB, AB1, AB2, IAB, ISYM
#ifdef _DEBUGPRINT_
#define _DBG_ .true.
#else
#define _DBG_ .false.
#endif
logical(kind=iwp), parameter :: LOCDBG = _DBG_
character(len=*), parameter :: SECNAM = 'CHO_MAXABSDIAG'

if (CHO_1CENTER) then ! specialization for 1-center approximation
  call CHO_MAXABSDIAG_1C(DIAG,IRED,DGMAX)
  return
end if

if (IRED == 1) then
  do ISYM=1,NSYM
    if (NNBSTR(ISYM,IRED) < 1) then
      DIAMAX(ISYM) = Zero
    else
      DIAMAX(ISYM) = abs(DIAG(IIBSTR(ISYM,IRED)+1))
      AB1 = IIBSTR(ISYM,IRED)+2
      AB2 = IIBSTR(ISYM,IRED)+NNBSTR(ISYM,IRED)
      do AB=AB1,AB2
        DIAMAX(ISYM) = max(DIAMAX(ISYM),abs(DIAG(AB)))
      end do
    end if
    DIAMAXT(ISYM) = DIAMAX(ISYM)
  end do
else if ((IRED == 2) .or. (IRED == 3)) then
  do ISYM=1,NSYM
    if (NNBSTR(ISYM,IRED) < 1) then
      DIAMAX(ISYM) = Zero
    else
      AB = INDRED(IIBSTR(ISYM,IRED)+1,IRED)
      DIAMAX(ISYM) = abs(DIAG(AB))
      AB1 = IIBSTR(ISYM,IRED)+2
      AB2 = IIBSTR(ISYM,IRED)+NNBSTR(ISYM,IRED)
      do IAB=AB1,AB2
        AB = INDRED(IAB,IRED)
        DIAMAX(ISYM) = max(DIAMAX(ISYM),abs(DIAG(AB)))
      end do
    end if
    if (NNBSTR(ISYM,1) < 1) then
      DIAMAXT(ISYM) = Zero
    else
      DIAMAXT(ISYM) = abs(DIAG(IIBSTR(ISYM,1)+1))
      AB1 = IIBSTR(ISYM,1)+2
      AB2 = IIBSTR(ISYM,1)+NNBSTR(ISYM,1)
      do AB=AB1,AB2
        DIAMAXT(ISYM) = max(DIAMAXT(ISYM),abs(DIAG(AB)))
      end do
    end if
  end do
else
  write(LUPRI,*) SECNAM,': unknown reduced set, IRED = ',IRED
  call CHO_QUIT('Unknown reduced set in '//SECNAM,104)
end if

DGMAX = DIAMAX(1)
do ISYM=2,NSYM
  DGMAX = max(DGMAX,DIAMAX(ISYM))
end do

if (LOCDBG) then
  write(LUPRI,*) SECNAM,': in reduced set ',IRED,':'
  write(LUPRI,*) 'DIAMAX  = ',(DIAMAX(ISYM),ISYM=1,NSYM)
  write(LUPRI,*) 'DIAMAXT = ',(DIAMAXT(ISYM),ISYM=1,NSYM)
  write(LUPRI,*) 'DGMAX   = ',DGMAX
end if

end subroutine CHO_MAXABSDIAG
