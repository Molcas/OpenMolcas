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

subroutine SG_SETUP_CASPT2()

use Molcas, only: MxLev
use caspt2_qmc_interface, only: DoFCIQMC
use RefWfn, only: L2Act, Level
use sguga, only: SG_Init, SG_Init_Simple
use sguga, only: CIS
use general_data, only: iSpin, nActel, nAsh, nEle3=>nElec3, nHole1, STSym, nLev
use caspt2_module, only: DMRG, DoCumulant, iSCF, MxCI, nRas1, nRas2, nRas3, nSym
use stdalloc, only: mma_allocate
use general_data, only: nRas,nRasEl,nRsPrt

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ISM(MxLev), ISYM, IT, nRs1T
integer(kind=iwp), parameter :: iState=1

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

  call SG_Init(iState,nSym,nActEl,iSpin,             &
               nRas,nRasEl,nRsPrt,                         &
               xLevel=Level,xL2Act=L2Act,xnLev=nLev,xNSM=ISM)

else

  call SG_Init_Simple(iState,nSym,nActEl,iSpin,       &
                     nRas,nRasEl,nRsPrt,                   &
                     xLevel=Level,xL2Act=L2Act,xnLev=nLev,             &
                     xNSM=ISM,Do_MkSGuga=.false.)
  CIS(iState)%NCSF(STSYM) = 1
end if

MXCI = maxval(CIS(iState)%NCSF(1:NSYM))

! NOTE: AT THIS POINT, WE HAVE ALLOCATED MEMORY SPACE FOR SGUGA USE:
! MVL,MVR,NOW,IOW,NOCP,IOCP,NOCSF,IOCSF,ICASE,ICOUP,VTAB,TMP.

end subroutine SG_SETUP_CASPT2
