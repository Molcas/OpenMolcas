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

subroutine MKTDAO(CMO,TDMO,TDAO,SCR)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
#include "mrci.fh"
real(kind=wp) :: CMO(NCMO), TDMO(NBAST,NBAST), TDAO(NBAST,NBAST), SCR(NBMAX,NBMAX)
integer(kind=iwp) :: I, IEB, IECMO, IEO, II, ISB, ISCMO, ISCO, ISO, ISYM, J, JJ, NB, NBCO, NBD, NBF, NCO, ND, NF, NO

! REORDER TDMO (USE TDAO AS TEMPORARY STORAGE):
call FZERO(TDAO,NBAST**2)
do I=1,NORBT
  II = ICH(I)
  if (II <= 0) goto 10
  do J=1,NORBT
    JJ = ICH(J)
    if (JJ <= 0) goto 5
    TDAO(I,J) = TDMO(II,JJ)
5   continue
  end do
10 continue
end do
call DCOPY_(NBAST**2,TDAO,1,TDMO,1)
call FZERO(TDAO,NBAST**2)
IECMO = 0
IEO = 0
IEB = 0
do ISYM=1,NSYM
  ISO = IEO+1
  ISB = IEB+1
  NO = NORB(ISYM)
  NB = NBAS(ISYM)
  IEO = IEO+NO
  IEB = IEB+NB
  ! ORBITALS PRE-FROZEN IN MOTRA, OR FROZEN IN MRCI:
  NF = NFMO(ISYM)+NFRO(ISYM)
  NBF = NB*NF
  ISCMO = IECMO+1
  IECMO = IECMO+NBF
  ! ORBITALS EXPLICITLY USED IN CI:
  NCO = NO-NFRO(ISYM)-NDEL(ISYM)
  ISCO = ISO+NFRO(ISYM)
  NBCO = NB*NCO
  ISCMO = IECMO+1
  IECMO = IECMO+NBCO
  if (NCO > 0) then
    call DGEMM_('N','N',NB,NCO,NCO,One,CMO(ISCMO),NB,TDMO(ISCO,ISCO),NBAST,Zero,SCR,NB)
    call DGEMM_('N','T',NB,NB,NCO,One,SCR,NB,CMO(ISCMO),NB,One,TDAO(ISB,ISB),NBAST)
  end if
  ! ORBITALS PRE-DELETED IN MOTRA OR DELETED IN MRCI:
  ND = NDMO(ISYM)+NDEL(ISYM)
  NBD = NB*ND
  ISCMO = IECMO+1
  IECMO = IECMO+NBD
end do

return

end subroutine MKTDAO
