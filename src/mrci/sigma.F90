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

subroutine SIGMA(SGM,AREF,CI,INTSY,INDX,BMN,IBMN,BIAC2,BICA2,BFIN3,ISAB,AC1,AC2,BFIN4,ABIJ,AIBJ,AJBI,ASCR1,BSCR1,FSCR1,FSEC,FOCK, &
                 ASCR2,BSCR2,FSCR2,DBK)

use mrci_global, only: ESHIFT, FIJKL, GFAC, ICPF, IFIRST, IREFX, IREST, ITER, NCONF, NREF, POTNUC
use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: SGM(*), BMN(*), BIAC2(*), BICA2(*), BFIN3(*), AC1(*), AC2(*), BFIN4(*), ASCR1(*), BSCR1(*), &
                                FSCR1(*), ASCR2(*), BSCR2(*), FSCR2(*), DBK(*)
real(kind=wp), intent(in) :: AREF(*)
real(kind=wp), intent(inout) :: CI(*), ABIJ(*), AIBJ(*), AJBI(*), FSEC(*), FOCK(*)
integer(kind=iwp), intent(in) :: INTSY(*), INDX(*), ISAB(*)
integer(kind=iwp), intent(_OUT_) :: IBMN(*)
integer(kind=iwp) :: ICSF, IREF, KTYP
real(kind=wp) :: GINV, SQG, SQGP

SGM(1:NCONF) = Zero

call CSFTRA(' CSF',CI,AREF)
SQGP = One
SQG = One
if (ICPF == 1) then
  SQGP = sqrt(GFAC)
  SQG = One/SQGP
  do IREF=1,NREF
    ICSF = IREFX(IREF)
    CI(ICSF) = SQGP*CI(ICSF)
  end do
end if

call DIAGC(INTSY,CI,SGM)
if ((IFIRST == 0) .and. ((IREST == 1) .or. (ITER > 1))) then
  call ABCI(INTSY,INDX,CI,SGM,BMN,IBMN,BIAC2,BICA2,BFIN3)
  call ABCD(INTSY,INDX,ISAB,CI,SGM,AC1,AC2,BFIN4)
end if
call IJKL(INTSY,INDX,CI,SGM,FIJKL)

call FAIBJ(INTSY,INDX,CI,SGM,ABIJ,AIBJ,AJBI,ASCR1,BSCR1,FSCR1,FSEC)

if (ITER > 0) then
  KTYP = 1
  ! Switch KTYP=1 means AI is actually handling AIJK integrals.
  call AI_MRCI(INTSY,INDX,CI,SGM,FOCK,ASCR2,BSCR2,FSCR2,DBK,KTYP)
end if
call FIJ(INTSY,INDX,CI,SGM,FOCK,ASCR2,BSCR2,FSCR2,DBK)

SGM(1:NCONF) = SGM(1:NCONF)+(POTNUC-ESHIFT)*CI(1:NCONF)
if (ICPF == 1) then
  GINV = One/GFAC
  SGM(1:NCONF) = GINV*SGM(1:NCONF)
  do IREF=1,NREF
    ICSF = IREFX(IREF)
    CI(ICSF) = SQG*CI(ICSF)
    SGM(ICSF) = SQGP*SGM(ICSF)
  end do
end if

call CSFTRA('MCSF',SGM,AREF)

return

end subroutine SIGMA
