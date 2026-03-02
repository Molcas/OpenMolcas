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
Subroutine mkTraCI(nTORB,TORB,STSYM,nConf,CI)
use caspt2_module, only: nSym, nIsh, nAsh, nRas1, nRas2, nRas3, nSsh, nOrb, nBas, nAES

use definitions, only: iwp, wp
integer(kind=iwp), intent(in):: nTORB,STSYM,nCONF
real(kind=wp), intent(inout):: TORB(nTORB)
real(kind=wp), intent(inout):: CI(nConf)

integer(kind=iwp) ITOEND,ISYM,NI,NA,NR1,NR2,NR3,NS,NO,NB,ITOSTA,ITO,iSTART

ITOEND=0
DO ISYM=1,NSYM
   NI=NISH(ISYM)
   NA=NASH(ISYM)
   NR1=NRAS1(ISYM)
   NR2=NRAS2(ISYM)
   NR3=NRAS3(ISYM)
   NS=NSSH(ISYM)
   NO=NORB(ISYM)
   NB=NBAS(ISYM)
   ITOSTA=ITOEND+1
   ITOEND=ITOEND+NI**2+NR1**2+NR2**2+NR3**2+NS**2

   ITO=ITOSTA+NI**2
   IF (NA<=0) CYCLE
   IF (NR1>0) THEN
      ISTART=NAES(ISYM)+1
      CALL TRACI_RPT2(ISTART,NR1,TORB(ITO),STSYM,NCONF,CI)
   END IF
   ITO=ITO+NR1**2
   IF (NR2>0) THEN
      ISTART=NAES(ISYM)+NR1+1
      CALL TRACI_RPT2(ISTART,NR2,TORB(ITO),STSYM,NCONF,CI)
   END IF
   ITO=ITO+NR2**2
   IF(NR3>0) THEN
      ISTART=NAES(ISYM)+NR1+NR2+1
      CALL TRACI_RPT2(ISTART,NR3,TORB(ITO),STSYM,NCONF,CI)
   END IF
END DO
End Subroutine mkTraCI

