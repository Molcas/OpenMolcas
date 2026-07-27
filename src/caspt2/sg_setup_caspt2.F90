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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
<<<<<<< HEAD
SUBROUTINE SG_SETUP_CASPT2()

use Molcas, only: MxLev
use fciqmc_interface, only: DoFCIQMC
use stdalloc, only: mma_allocate
use RefWfn, only: L2Act, Level
use sguga, only: SGS, CIS, EXS, SG_Init, SG_Init_Simple
use caspt2_module, only: DMRG, DoCumulant, iSCF, iSpin, nActEl,        &
                         nEle3, nHole1, nRas1, nRas2, nRas3, nSym, STSym
use caspt2_module, only: nAsh, MxCI
use definitions, only: iwp
IMPLICIT NONE

INTEGER(kind=iwp) I,IT,nLEV,ILEV,ISYM,ISM(MxLev)

nLEV=0
DO ISYM=1,NSYM
   DO IT=1,NASH(ISYM)
      nLEV=nLEV+1
      ILEV=LEVEL(nLEV)
      ISM(ILEV)=ISYM
   END DO
END DO

if ((.NOT.DoCumulant) .and. (nactel.gt.0) .and. (iscf.eq.0)            &
      .and. (.not. DoFCIQMC) .and. (.not. DMRG)) Then

   call SG_Init(nSym,nActEl,iSpin,                                     &
               SGS,CIS,EXS,                                            &
               nHole1,nEle3,                                           &
               nRas1,nRas2,nRas3,                                      &
               xLevel=Level,xL2Act=L2Act,                              &
               xnLev=nLev,xNSM=ISM)

else

   call SG_Init_Simple(nSym,nActEl,iSpin,SGS,CIS,                      &
                       xLevel=Level,xL2Act=L2Act,                      &
                       xnLev=nLev,xNSM=ISM,Do_MkSGuga=.FALSE.)
   SGS%iSpin=0
   SGS%nActEl=0

! INITIALIZE SPLIT-GRAPH GUGA DATA SETS:
   Call mma_allocate(CIS%NCSF,SGS%nSym,Label='CIS%NCSF')
   CIS%NCSF(:)=0
   CIS%NCSF(STSYM)=1
   Call mma_allocate(EXS%ICoup,[1,3],[1,1],Label='EXS%ICoup')
   Call mma_allocate(EXS%VTab,[1,1],Label='EXS%VTab')
endif

MXCI=1
DO I=1,NSYM
   MXCI=MAX(MXCI,CIS%NCSF(I))
END DO
=======

subroutine SG_SETUP_CASPT2()

use Molcas, only: MxLev
use caspt2_qmc_interface, only: DoFCIQMC
use RefWfn, only: L2Act, Level
use sguga, only: SG_Init, SG_Init_Simple
use sguga_states, only: CIS, EXS, SGS
use general_data, only: iSpin, nActel, nAsh, nEle3=>nElec3, nHole1, STSym, nLev
use caspt2_module, only: DMRG, DoCumulant, iSCF, MxCI, nRas1, nRas2, nRas3, nSym
use stdalloc, only: mma_allocate
use rasdef, only: nRas,nRasEl,nRsPrt

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ISM(MxLev), ISYM, IT, nRs1T
integer(kind=iwp), parameter :: istate=1

nLEV = 0
do ISYM=1,NSYM
  do IT=1,NASH(ISYM)
    nLEV = nLEV+1
    ISM(LEVEL(nLEV)) = ISYM
  end do
end do

If (nHole1+nEle3/=0) Then
   nRsPrt=3
   nRas(:,1)=nRas1(:)
   nRas(:,2)=nRas2(:)
   nRas(:,3)=nRas3(:)
   nRs1T=Sum(nRas1(1:nSym))
   nRasEl(1)=2*nRs1T-nHole1
   nRasEl(2)=nActel-nEle3
   nRasEl(3)=nActel
Else
   nRsPrt=1
   nRas(:,1)=nRas2(:)
   nRasEl(1)=nActel
End If

if ((.not. DoCumulant) .and. (nactel > 0) .and. (iscf == 0) .and. (.not. DoFCIQMC) .and. (.not. DMRG)) then

  call SG_Init(nSym,nActEl,iSpin,SGS(istate),CIS(istate),                             &
               nRas,nRasEl,nRsPrt,EXS(istate),                                &
               xLevel=Level,xL2Act=L2Act,xnLev=nLev,xNSM=ISM)

else

  call SG_Init_Simple(nSym,nActEl,iSpin,SGS(istate),CIS(istate),                       &
                     nRas,nRasEl,nRsPrt,                               &
                     xLevel=Level,xL2Act=L2Act,xnLev=nLev,             &
                     xNSM=ISM,Do_MkSGuga=.false.)
  SGS(istate)%iSpin = 0
  SGS(istate)%nActEl = 0

  ! INITIALIZE SPLIT-GRAPH GUGA DATA SETS:
  call mma_allocate(CIS(istate)%NCSF,SGS(istate)%nSym,Label='CIS%NCSF')
  CIS(istate)%NCSF(:) = 0
  CIS(istate)%NCSF(STSYM) = 1
  call mma_allocate(EXS(istate)%ICoup,[1,3],[1,1],Label='EXS%ICoup')
  call mma_allocate(EXS(istate)%VTab,[1,1],Label='EXS%VTab')
end if

MXCI = maxval(CIS(istate)%NCSF(1:NSYM))
>>>>>>> upstream-openmolcas/master

! NOTE: AT THIS POINT, WE HAVE ALLOCATED MEMORY SPACE FOR SGUGA USE:
! MVL,MVR,NOW,IOW,NOCP,IOCP,NOCSF,IOCSF,ICASE,ICOUP,VTAB,TMP.

<<<<<<< HEAD
END SUBROUTINE SG_SETUP_CASPT2
=======
end subroutine SG_SETUP_CASPT2
>>>>>>> upstream-openmolcas/master
