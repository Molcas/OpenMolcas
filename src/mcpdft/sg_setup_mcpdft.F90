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
! Copyright (C) 2026, Roland Lindh                                     *
!***********************************************************************

!> The CI state of MC-PDFT (CI_St(iState), CI_Info) for the whole run, set
!> up after the reference wave function has been read (Proc_InpX) and freed
!> at the end of MC-PDFT. MC-PDFT does no CI itself, so the state has no
!> graph (Mode 'NONE'); it holds the state's spin, electrons, symmetry and
!> levels. The levels are the active orbitals in symmetry-blocked order.
subroutine SG_Setup_MCPDFT()

use Molcas, only: MxLev
use CI_States, only: CI_Setup_RAS
use Orbital_Spaces, only: nAsh, nSym
use CI_Input, only: iSpin, nActEl, nElec3, nHole1, nRs1, nRs2, nRs3, STSym
use Definitions, only: iwp

implicit none
integer(kind=iwp), parameter :: iState = 1
integer(kind=iwp) :: iq, ISM(MxLev), iSym, it, Level(MxLev), nLev

nLev = 0
do iSym=1,nSym
  do it=1,nAsh(iSym)
    nLev = nLev+1
    ISM(nLev) = iSym
  end do
end do
Level(:) = [(iq,iq=1,MxLev)]

call CI_Setup_RAS(iState,nSym,nActEl,iSpin,STSym,nRs1(1:nSym),nRs2(1:nSym),nRs3(1:nSym),nHole1,nElec3,.false.,nLev, &
                  Level(1:nLev),Level(1:nLev),ISM(1:nLev),'NONE')

end subroutine SG_Setup_MCPDFT
