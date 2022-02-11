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
                 BFIN5,ASCR2,BSCR2,FSCR2,DBK,ICSPCK)
!PAM04 BFIN5,ASCR2,BSCR2,FSCR2,DBK,CSPCK)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: SGM(*), AREF(*), CI(*), BMN(*), BIAC2(*), BICA2(*), BFIN3(*), AC1(*), AC2(*), BFIN4(*), ABIJ(*), AIBJ(*), &
                 AJBI(*), ASCR1(*), BSCR1(*), FSCR1(*), FSEC(*), FOCK(*), BFIN5(*), ASCR2(*), BSCR2(*), FSCR2(*), DBK(*)
integer(kind=iwp) :: INTSY(*), INDX(*), IBMN(*), ISAB(*), ICSPCK(*)
#include "WrkSpc.fh"
#include "mrci.fh"
integer(kind=iwp) :: ICSF, IREF, KTYP
real(kind=wp) :: GINV, SQG, SQGP

call DCOPY_(NCONF,[Zero],0,SGM,1)

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
  call ABCI_MRCI(INTSY,INDX,CI,SGM,BMN,IBMN,BIAC2,BICA2,BFIN3)
  call ABCD_MRCI(INTSY,INDX,ISAB,CI,SGM,AC1,AC2,BFIN4)
end if
!PAM04 call IJKL(INTSY,INDX,CI,SGM,FIJKL)
call IJKL(INTSY,INDX,CI,SGM,WORK(LFIJKL))

!PAM04 call FAIBJ(INTSY,INDX,CI,SGM,ABIJ,AIBJ,AJBI,BFIN1,BFIN1,ASCR1,BSCR1,FSCR1,FSEC)
call FAIBJ(INTSY,INDX,CI,SGM,ABIJ,AIBJ,AJBI,ASCR1,BSCR1,FSCR1,FSEC)

if (ITER > 0) then
  KTYP = 1
  ! Switch KTYP=1 means AI is actually handling AIJK integrals.
  !PAM04 call AI(INTSY,INDX,CI,SGM,FOCK,BFIN5,BFIN5,ASCR2,BSCR2,FSCR2,DBK,KTYP)
  call AI_MRCI(INTSY,INDX,CI,SGM,FOCK,ASCR2,BSCR2,FSCR2,DBK,KTYP)
end if
!PAM04 call FIJ(CSPCK,INTSY,INDX,CI,SGM,FOCK,ASCR2,BSCR2,FSCR2,DBK)
call FIJ_MRCI(ICSPCK,INTSY,INDX,CI,SGM,FOCK,ASCR2,BSCR2,FSCR2,DBK)

call DAXPY_(NCONF,POTNUC-ESHIFT,CI,1,SGM,1)
if (ICPF == 1) then
  GINV = One/GFAC
  call DSCAL_(NCONF,GINV,SGM,1)
  do IREF=1,NREF
    ICSF = IREFX(IREF)
    CI(ICSF) = SQG*CI(ICSF)
    SGM(ICSF) = SQGP*SGM(ICSF)
  end do
end if

call CSFTRA('MCSF',SGM,AREF)

return
! Avoid unused argument warnings
if (.false.) call Unused_real_array(BFIN5)

end subroutine SIGMA
