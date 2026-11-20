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
use CI_Info, only: CI_St
use CI_States, only: CI_Setup_RAS
use Orbital_Spaces, only: nAsh
use CI_Input, only: iSpin, nActel, nElec3, nHole1, nLev, STSym
use caspt2_module, only: DMRG, DoCumulant, HZERO, iSCF, MxCI, nRas1, nRas2, nRas3, nSym
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ISM(MxLev), ISYM, IT
integer(kind=iwp), parameter :: istate = 1

nLEV = 0
do ISYM=1,NSYM
  do IT=1,NASH(ISYM)
    nLEV = nLEV+1
    ISM(LEVEL(nLEV)) = ISYM
  end do
end do

! NEVPT2 applies E_pq to the reference also for a single determinant (MKNEVE4).
if ((.not. DoCumulant) .and. (nactel > 0) .and. ((iscf == 0) .or. (HZERO == 'DYALL')) .and. (.not. DoFCIQMC) .and. (.not. DMRG)) &
  then

  call CI_Setup_RAS(iState,nSym,nActEl,iSpin,STSym,nRas1(1:nSym),nRas2(1:nSym),nRas3(1:nSym),nHole1,nElec3,.false.,nLev, &
                    Level(1:nLev),L2Act(1:nLev),ISM(1:nLev),'FULL')

else

  call CI_Setup_RAS(iState,nSym,nActEl,iSpin,STSym,nRas1(1:nSym),nRas2(1:nSym),nRas3(1:nSym),nHole1,nElec3,.false.,nLev, &
                    Level(1:nLev),L2Act(1:nLev),ISM(1:nLev),'NONE',LSM_Single=STSYM)
end if

MXCI = maxval(CI_St(iState)%nCSF(1:NSYM))

! NOTE: AT THIS POINT, WE HAVE ALLOCATED MEMORY SPACE FOR SGUGA USE:
! MVL,MVR,NOW,IOW,NOCP,IOCP,NOCSF,IOCSF,ICASE,ICOUP,VTAB,TMP.

end subroutine SG_SETUP_CASPT2
