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

!> The CI state of MCLR (CI_St(iState), CI_Info) for the whole run: set up
!> once after the reference wave function has been read (RdJobIph,
!> RdJobIph_td; Mode 'GRAPH': graph and walks, for the reorderings of
!> SG2SymG), freed at the end of MCLR.
subroutine SG_Setup_MCLR(Mode)

use molcas, only: MxLev
use CI_States, only: CI_Setup_RAS
use Orbital_Spaces, only: nSym
use CI_Input, only: iSpin, nActEl, nElec3, nHole1, nRS1, nRS2, nRS3, STSym
use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: Mode
integer(kind=iwp) :: iBas, iq, ISM(MxLev), iSym, Level(MxLev), nLev
integer(kind=iwp), parameter :: iState = 1

nLev = 0
do iSym=1,nSym
  do iBas=1,nRs1(iSym)
    nLev = nLev+1
    ISM(nLev) = iSym
  end do
end do
do iSym=1,nSym
  do iBas=1,nRs2(iSym)
    nLev = nLev+1
    ISM(nLev) = iSym
  end do
end do
do iSym=1,nSym
  do iBas=1,nRs3(iSym)
    nLev = nLev+1
    ISM(nLev) = iSym
  end do
end do

Level(1:MxLev)=[(iq,iq=1,MxLev)]

call CI_Setup_RAS(istate,nSym,nActEl,iSpin,STSym,nRs1(1:nSym),nRs2(1:nSym),nRs3(1:nSym),nHole1,nElec3,.false.,nLev,Level(1:nLev), &
                  Level(1:nLev),ISM(1:nLev),Mode)

end subroutine SG_Setup_MCLR
