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
! Copyright (C) 2014,2026, Ignacio Fdez. Galvan                        *
!***********************************************************************

! June 2026: Merge with TSHop/TSHinit functionality

!define _DEBUGPRINT_
subroutine Track_State(OVLP)

use rassi_aux, only: ipglob
use Cntrl, only: HOP, NJOB, NSTAT, NSTATE, TRACK
use rassi_global_arrays, only: HAM
use Constants, only: Zero, Quart
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: ovlp(nstate,nstate)
integer(kind=iwp) :: initState, iState, maxHop, newState, nHop
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: j
#endif
real(kind=wp) :: E0, Em, Ep, MaxOv, ThisOv
logical(kind=iwp) :: lAllowHop, lHop, lMaxHop
real(kind=wp), parameter :: Ethr = 0.03_wp

! Check that there are 2 JOB files, with the same number of states
if (nJob /= 2) call SysAbendMsg('Track_State','The number of JOB files should be 2.','')
if (nStat(2) /= nStat(1)) call SysAbendMsg('Track_State','The number of states in the JOB files should be the same.','')

! Get the root to track
call Get_iScalar('Relax CASSCF root',initState)

! Find the root in the first JOB (current states) that has maximum
! overlap with the tracked root in the second JOB (previous states)
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) 'OVERLAP MATRIX FOR THE ORIGINAL STATES:'
write(u6,*)
do iState=1,nState
  write(u6,'(5(1X,F15.8))') (Ovlp(iState,j),j=1,iState)
end do
#endif
if (IPGLOB >= 2) then
  write(u6,*)
  write(u6,*) 'Initial root: ',initState
  write(u6,*) 'Overlaps with current states:'
end if
MaxOv = Zero
newState = 0
do iState=1,nStat(1)
  ThisOv = Ovlp(iState,initState+nStat(1))
  if (IPGLOB >= 2) write(u6,'(I5,1X,F15.8)') iState,ThisOv
  if (TRACK .and. (abs(ThisOv) > MaxOv)) then
    MaxOv = abs(ThisOv)
    newState = iState
  end if
end do

! The HOP keyword is somewhat special:
! - Check only the immediately upper or lower states
! - Prefer hop down if both possible
! - Hard-coded threshold of 0.25 overlap
! - Hard-coded 0.03 energy difference threshold
if (HOP) then
  E0 = HAM(initState,initState)
  if (IPGLOB >= 2) then
    write(u6,'(//,6X,A)') repeat('*',100)
    write(u6,'(6X,A,98X,A)') '*','*'
    write(u6,'(6X,A,36X,A,37X,A)') '*',' Surface hopping section ','*'
    write(u6,'(6X,A,98X,A)') '*','*'
    write(u6,'(6X,A,//)') repeat('*',100)
    write(u6,'(6X,A,6X,ES15.6,A,/)') 'The initial state energy is:',E0,' a.u.'
  end if
  newState = initState
  if (initState < nStat(1)) then
    Ep = HAM(initState+1,initState+1)
    if (IPGLOB >= 2) then
      write(u6,'(6X,A,I2)') 'There is an upper root, which is: ',initState+1
      write(u6,'(6X,A,6X,ES15.6,A)') 'Its energy is:',Ep,' a.u.'
      write(u6,'(6X,A,12X,ES15.6,A,/)') 'Ediff = ',E0-Ep,' a.u.'
    end if
    if ((abs(E0-Ep) < Ethr) .and. (abs(Ovlp(initState+1,initState+nStat(1)))) >= Quart) newState = initState+1
  else if (IPGLOB >= 2) then
    write(u6,'(6X,A)') 'There is no upper root'
  end if
  if (initState > 1) then
    Em = HAM(initState-1,initState-1)
    if (IPGLOB >= 2) then
      write(u6,'(6X,A,I2)') 'There is a lower root, which is: ',initState-1
      write(u6,'(6X,A,6X,ES15.6,A)') 'Its energy is:',Em,' a.u.'
      write(u6,'(6X,A,12X,ES15.6,A,/)') 'Ediff = ',E0-Em,' a.u.'
    end if
    if ((abs(E0-Em) < Ethr) .and. (abs(Ovlp(initState-1,initState+nStat(1))) >= Quart)) newState = initState-1
  else if (IPGLOB >= 2) then
    write(u6,'(6X,A)') 'There is no lower root'
  end if
  if (newState /= initState) then
    write(u6,'(/,6X,3A)') '+',repeat('-',78),'+'
    write(u6,'(6X,A1,T86,A1)') '|','|'
    write(u6,'(6X,A1,T35,A,T86,A1)') '|','A HOP event is detected!','|'
    write(u6,'(6X,A1,T86,A1)') '|','|'
    write(u6,'(6X,A1,T32,2(A,I3,4X),T86,A1)') '|','From state:',initState,'To state:',newState,'|'
    ! Check if the number of Hops is limited:
    call qpg_iScalar('MaxHops',lMaxHop)
    if (lMaxHop) then
      call Get_iScalar('MaxHops',maxHop)
      if (maxHop < 1) lMaxHop = .false.
    end if
    if (lMaxHop) then
      call qpg_iScalar('Number of Hops',lHop)
      if (lHop) then
        call Get_iScalar('Number of Hops',nHop)
      else
        nHop = 0
      end if
      if (maxHop <= nHop) then
        lAllowHop = .false.
        write(u6,'(6X,A1,T40,A,T86,A1)') '|','maxHop > nHop','|'
        write(u6,'(6X,A1,T31,A,T86,A1)') '|','This surface HOP is not allowed','|'
        write(u6,'(6X,A1,T24,A,T86,A1)') '|','because the number of allowed Hops is exceeded','|'
      else
        lAllowHop = .true.
      end if
    else
      lAllowHop = .true.
    end if
    ! Bottom of the printed box
    write(u6,'(6X,A1,T86,A1)') '|','|'
    write(u6,'(6X,3A,/)') '+',repeat('-',78),'+'
    ! Set the number of Hops
    if (lAllowHop) then
      nHop = nHop+1
      call Put_iScalar('Number of Hops',nHop)
    else
      newState = initState
    end if
  end if
end if

if (IPGLOB >= 2) write(u6,*) 'New root: ',newState

! If no state is found, something wrong has happened
if (newState == 0) call SysAbendMsg('Track_State','No overlaps!','')

! Store the new state for geometry optimization
if (newState /= initState) then
  call Put_iScalar('Relax CASSCF root',newState)
  call Put_iScalar('Relax Original root',newState)
  call Put_iScalar('NumGradRoot',newState)
end if

end subroutine Track_State
