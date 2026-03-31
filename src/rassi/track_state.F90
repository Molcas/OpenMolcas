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

!define _DEBUGPRINT_
subroutine Track_State(OVLP)

use rassi_aux, only: ipglob
use Cntrl, only: NJOB, NSTAT, NSTATE
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: ovlp(nstate,nstate)
integer(kind=iwp) :: initState, iState, newState
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: j
#endif
real(kind=wp) :: MaxOv, ThisOv

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
  if (abs(ThisOv) > MaxOv) then
    MaxOv = abs(ThisOv)
    newState = iState
  end if
end do
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
