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

subroutine SG_setup_RASSI(nSym,nActEl,iSpin,iState,nRs1,nRs2,nRs3,nHole1,nElec3)

use Molcas, only: MxLev
use CI_States, only: CI_Setup_RAS
use rassi_aux, only: Level
use rassi_data, only: NASH
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nSym, nActEl, iSpin, iState, nRs1(nSym), nRs2(nSym), nRs3(nSym), nHole1, nElec3
integer(kind=iwp) :: ILEV, iq, ISM(MxLev), ISYM, IT, L2Act(MxLev), nLev

nLev = 0
do ISYM=1,NSYM
  do IT=1,NASH(ISYM)
    nLev = nLev+1
    ILEV = Level(nLev)
    ISM(ILEV) = ISYM
  end do
end do

L2Act(1:MxLev) = [(iq,iq=1,MxLev)]

! RASSI always uses the three RAS partitions
call CI_Setup_RAS(iState,nSym,nActEl,iSpin,nRs1,nRs2,nRs3,nHole1,nElec3,.true.,nLev,Level(1:nLev),L2Act(1:nLev),ISM(1:nLev),'FULL')

end subroutine SG_setup_RASSI
