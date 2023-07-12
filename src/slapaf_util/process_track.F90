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
! Copyright (C) 2016, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine Process_Track()

use Slapaf_Info, only: Request_RASSI, RootMap
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
#include "print.fh"
integer(kind=iwp) :: i, MaxId(2), nOv, nRoots
logical(kind=iwp) :: Done, Found
character(len=8) :: Method
integer(kind=iwp), allocatable :: OldMap(:), RootIdx(:)
real(kind=wp), allocatable :: Overlaps(:,:), Ovlp(:)

call Get_cArray('Relax Method',Method,8)
if ((Method /= 'CASSCF') .and. (Method /= 'RASSCF') .and. (Method /= 'CASSCFSA') .and. (Method /= 'RASSCFSA') .and. &
    (Method /= 'CASPT2') .and. (Method /= 'RASPT2')) then
  call WarningMessage(2,'Error in Process_Track')
  write(u6,*) '***************** ERROR ********************'
  write(u6,*) ' The TRACK keyword can only be used with'
  write(u6,*) ' states computed by the RASSCF or CASPT2'
  write(u6,*) ' programs.'
  write(u6,*) '********************************************'
  call Quit_OnUserError()
end if

! Find the number of roots, and whether state overlaps are available

nRoots = 1
call Qpg_iScalar('Number of roots',Found)
if (Found) call Get_iScalar('Number of roots',nRoots)
call qpg_dArray('State Overlaps',Found,nOv)
if (Found .and. (nOv == (2*nRoots)**2)) then
  Request_RASSI = .false.
else
  Request_RASSI = .true.
end if

! Make sure that the root mapping is only done once per iteration

call Qpg_iScalar('Track Done',Done)
call Get_lScalar('Track Done',Done)
if (Request_RASSI .or. Done) return

! Modify the root map according to the overlaps:
!
!   RootMap(i) = N
!
! where i is the original root number (at iter=1), and N is the
! root number at the current iteration.

call mma_allocate(Ovlp,nOv)
call mma_allocate(Overlaps,nRoots,nRoots)
call get_dArray('State Overlaps',Ovlp,nOv)
do i=1,nRoots
  Overlaps(:,i) = Ovlp((2*i-1)*nRoots+1:2*i*nRoots)
end do
call mma_deallocate(Ovlp)
if (nPrint(1) >= 5) then
  call RecPrt('Overlaps with previous states','',Overlaps,nRoots,nRoots)
  call mma_allocate(OldMap,nRoots)
  OldMap(:) = RootMap(1:nRoots)
end if
call mma_allocate(RootIdx,nRoots)
do i=1,nRoots
  MaxId = maxloc(abs(Overlaps))
  RootIdx(MaxId(1)) = MaxId(2)
  Overlaps(MaxId(1),:) = Zero
  Overlaps(:,MaxId(2)) = Zero
end do
do i=1,nRoots
  RootMap(i) = RootIdx(RootMap(i))
end do
call Put_iArray('Root Mapping',RootMap,nRoots)
if (nPrint(1) >= 5) then
  write(u6,*)
  write(u6,100) 'Root map'
  write(u6,*)
  write(u6,100) 'Original  Prev.  This'
  write(u6,100) '  root    iter.  iter.'
  write(u6,100) '----------------------'
  do i=1,nRoots
    write(u6,101) i,OldMap(i),RootMap(i)
  end do
  write(u6,100) '----------------------'
  write(u6,*)
  call mma_deallocate(OldMap)
end if

call mma_deallocate(RootIdx)
call mma_deallocate(Overlaps)

! Update the RunFile for automatic gradients

call Qpg_iScalar('Relax CASSCF root',Found)
if (Found) then
  call Get_iScalar('Relax CASSCF root',i)
  if (RootMap(i) /= i) then
    call Put_iScalar('Relax CASSCF root',RootMap(i))
    call Put_iScalar('Relax Original root',RootMap(i))
  end if
end if
call Qpg_iScalar('NumGradRoot',Found)
if (Found) then
  call Get_iScalar('NumGradRoot',i)
  if (RootMap(i) /= i) call Put_iScalar('NumGradRoot',RootMap(i))
end if
call Put_lScalar('Track Done',.true.)

return

100 format(3X,A)
101 format(3X,I6,1X,I6,1X,I6)

end subroutine Process_Track
