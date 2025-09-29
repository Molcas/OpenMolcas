!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2001, Jeppe Olsen                                      *
!               2001, Giovanni Li Manni                                *
!               2001, Dongxia Ma                                       *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine NSTPTP_GAS(NGAS,ISPGRP,NSTSGP,NSMST,NSTSSPGP,IGRP,MXNSTR,NSMCLS,NSMCLSE,NSMCLSE1)
! Find number of strings per symmetry for the supergroup defined
! by the groups of ISPGRP. The obtained number of strings per sym
! is stored in NSTSSPGP(*,IGRP)
!
! Jeppe Olsen, Giovanni Li Manni, Dongxia Ma
! Dicember 2011 - the old version is too slow for many GAS spaces
!                 (new version simpler and quicker)
!
! Also delivered:
!
! NSMCLS : MAX Number of symmetry classes for given supergroup,
!          i.e. number of combinations of symmetries of groups
!          containing strings
! NSMCLSE : Number of symmetry classes for given supergroup
!           obtained by restricting allowed symmetries in
!           a given group by a max and min.
! NSMCLSE1 : As NSMCLSE, but the symmetry of the last active
!            orbital space where there is more than one symmetry
!            is left out

use Symmetry_Info, only: Mul
use lucia_data, only: MXPNGAS, MXPNSMST
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NGAS, ISPGRP(NGAS), NSMST, NSTSGP(NSMST,*), IGRP
integer(kind=iwp), intent(inout) :: NSTSSPGP(NSMST,IGRP)
integer(kind=iwp), intent(out) :: MXNSTR, NSMCLS, NSMCLSE, NSMCLSE1
integer(kind=iwp) :: IGAS, ISM, ISM1(MXPNSMST), ISM2(MXPNSMST), ISM_IGAS, ISM_IGASM1, ISYM, MNSM(MXPNGAS), MSM1(MXPNSMST), &
                     MSM2(MXPNSMST), MXSM(MXPNGAS)

#ifdef _DEBUGPRINT_
write(u6,*) ' ======================'
write(u6,*) ' NSTPTP_GAS is speaking'
write(u6,*) ' ======================'

write(u6,*) ' Supergroup in action'
call IWRTMA(ISPGRP,1,NGAS,1,NGAS)
#endif

! The NSMCLS* parameters

! Max and min for each GASpace

do IGAS=1,NGAS
  MXSM(IGAS) = 1
  do ISYM=1,NSMST
    if (NSTSGP(ISYM,ISPGRP(IGAS)) /= 0) MXSM(IGAS) = ISYM
  end do
  MNSM(IGAS) = NSMST
  do ISYM=NSMST,1,-1
    if (NSTSGP(ISYM,ISPGRP(IGAS)) /= 0) MNSM(IGAS) = ISYM
  end do
end do
! NSMCLSE
NSMCLSE = product(MXSM(1:NGAS)-MNSM(1:NGAS)+1)
! NSMCLSE1
NSMCLSE1 = NSMCLSE
do IGAS=1,NGAS
  ! In ISM1, the number of strings per symmetry for the first
  ! IGAS-1 spaces are given, obtain in ISM2 the number of strings per sym
  ! for the first IGAS spaces
  ! Also: in MSM1, MSM2, counts the number of nontrivial combinations per sym
  if (IGAS == 1) then
    ! ISM1: The number of strings per symmetry for zero electrons
    ISM1(1:NSMST) = 0
    ISM1(1) = 1
    MSM1(1:NSMST) = 0
    MSM1(1) = 1
  else
    ! copy from the ISM2 obtained for preceding IGAS
    ISM1(1:NSMST) = ISM2(1:NSMST)
    MSM1(1:NSMST) = MSM2(1:NSMST)
  end if
  ISM2(1:NSMST) = 0
  MSM2(1:NSMST) = 0
  do ISM_IGASM1=1,NSMST
    do ISM_IGAS=1,NSMST
      ISM = Mul(ISM_IGASM1,ISM_IGAS)
      ISM2(ISM) = ISM2(ISM)+ISM1(ISM_IGASM1)*NSTSGP(ISM_IGAS,ISPGRP(IGAS))
      if (ISM1(ISM_IGASM1)*NSTSGP(ISM_IGAS,ISPGRP(IGAS)) /= 0) MSM2(ISM) = MSM2(ISM)+MSM1(ISM_IGASM1)
    end do
  end do
end do !loop over IGAS
NSTSSPGP(:,IGRP) = ISM2(1:NSMST)

MXNSTR = max(0,maxval(NSTSSPGP(:,IGRP)))
NSMCLS = max(0,maxval(MSM2(1:NSMST)))

#ifdef _DEBUGPRINT_
write(u6,*) ' Number of strings per symmetry for supergroup',IGRP
call IWRTMA10(NSTSSPGP(:,IGRP),1,NSMST,1,NSMST)
write(u6,*) ' Largest number of strings of given sym ',MXNSTR

write(u6,'(A,3(2X,I8))') ' NSMCLS,NSMCLSE,NSMCLSE1=',NSMCLS,NSMCLSE,NSMCLSE1
#endif

end subroutine NSTPTP_GAS
