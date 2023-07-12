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
! Copyright (C) 2015,2016, Ignacio Fdez. Galvan                        *
!***********************************************************************

subroutine Process_Gradients()

use Slapaf_Info, only: ApproxNADC, Energy, Energy0, Gx, Gx0, iState, iter, NAC, NADC, Request_Alaska, RootMap, TwoRunFiles
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: Columbus, i, nRoots, nsAtom, RC
real(kind=wp) :: E0, E1
logical(kind=iwp) :: Exists, Found
real(kind=wp), allocatable :: Ener(:), Grads(:,:,:)
integer(kind=iwp), external :: Read_Grad

Request_Alaska = .false.
nsAtom = size(Gx,2)
call mma_Allocate(Grads,3,nsAtom,3)

! First check that all the needed gradients are available
! and stop to compute them if they are not.
! For a two-RunFile job, this behaves as if it is one state.

if (TwoRunFiles) iState(:) = 0

if (iState(1) /= 0) iState(1) = RootMap(iState(1))
if (iState(2) /= 0) iState(2) = RootMap(iState(2))
i = max(iState(1),iState(2))
iState(2) = min(iState(1),iState(2))
iState(1) = i

if ((iState(1) /= 0) .and. (iState(2) /= 0)) then

  ! Two states

  do i=2,1,-1
    RC = Read_Grad(Grads(:,:,i),3*nsAtom,iState(i),0,0)
    if (RC == 0) then
      Request_Alaska = .true.
      call Put_iScalar('Relax CASSCF root',iState(i))
      call Put_iScalar('NumGradRoot',iState(i))
      iState(1) = iState(i)
      iState(2) = 0
      exit
    end if
  end do
  if ((.not. Request_Alaska) .and. NADC) then
    RC = Read_Grad(Grads(:,:,3),3*nsAtom,0,iState(1),iState(2))
    if (RC == 0) Request_Alaska = .true.
  end if

else

  ! One state

  iState(1) = 0
  iState(2) = 0
  call Qpg_iScalar('Relax CASSCF root',Exists)
  if (Exists) call Get_iScalar('Relax CASSCF root',iState(1))
  if (iState(1) == 0) iState(1) = 1
  RC = Read_Grad(Grads(:,:,1),3*nsAtom,iState(1),0,0)
  if (RC == 0) Request_Alaska = .true.
end if

if (Request_Alaska) then
  call mma_Deallocate(Grads)
  NADC = .false.
  return
end if

! Once the gradients are read, convert them to forces
! and store them as appropriate.

! Energy and gradient of the first (higher) state

nRoots = 1
call Qpg_iScalar('Number of roots',Found)
if (Found) call Get_iScalar('Number of roots',nRoots)
call mma_Allocate(Ener,nRoots)
call Get_dArray('Last energies',Ener,nRoots)
if (max(iState(1),iState(2)) > nRoots) then
  call WarningMessage(2,'Too few energies in RUNFILE')
  call Abend()
end if
Energy(iter) = Ener(iState(1))
E1 = Ener(iState(1))
Gx(:,:,iter) = -Grads(:,:,1)

! For a two-RunFile job, read the second energy
! and gradient from RUNFILE2

if (TwoRunFiles) then
  call NameRun('RUNFILE2')
  iState(2) = 0
  call Qpg_iScalar('Relax CASSCF root',Exists)
  if (Exists) call Get_iScalar('Relax CASSCF root',iState(2))
  if (iState(2) == 0) iState(2) = 1
  nRoots = 1
  call Qpg_iScalar('Number of roots',Found)
  if (Found) call Get_iScalar('Number of roots',nRoots)
  call mma_Deallocate(Ener)
  call mma_Allocate(Ener,nRoots)
  call Get_dArray('Last energies',Ener,nRoots)
  call Get_dArray('GRAD',Grads(:,:,2),3*nsAtom)
  call NameRun('#Pop')
  RC = -1
end if

if (iState(2) > 0) then

  ! With two states the Lagrangian is different!
  ! We will optimize on the average energy and have a constraint
  ! for the energy difference. We change the energies and gradients
  ! on the fly.

  E0 = Ener(iState(2))
  Energy(iter) = (E1+E0)*Half
  Energy0(iter) = E1-E0
  Gx(:,:,iter) = Half*(Gx(:,:,iter)-Grads(:,:,2))
  Gx0(:,:,iter) = Grads(:,:,2)-Grads(:,:,1)
  call Get_iScalar('Columbus',Columbus)

  ! In case of a true conical intersection there is also coupling.

  if (NADC) then
    call Get_iScalar('Columbus',Columbus)
    if (Columbus /= 1) then
      NAC(:,:,iter) = Grads(:,:,3)

      ! If the coupling derivative vector could not be calculated,
      ! use an approximate one.

      if (RC < 0) then
        ApproxNADC = .true.
        call Branching_Plane_Update(Gx,Gx0,NAC(:,:,iter),3*nsAtom,iter)
      end if
    end if
  end if
end if

call mma_Deallocate(Ener)
call mma_Deallocate(Grads)

return

end subroutine Process_Gradients
