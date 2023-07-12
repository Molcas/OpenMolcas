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

subroutine Init2()

use Slapaf_Info, only: Coor, Cx, DipM, dqInt, dqInt_Aux, Energy, Energy0, Get_Slapaf, Grd, Gx, Gx0, iter, lOld_Implicit, MaxItr, &
                       mTROld, NAC, NADC, qInt, RefGeo, TwoRunFiles
use Kriging_Mod, only: nSet
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: Columbus, i, iDummy, iMode, i_N, iRoot, j, Length, mLambda, nData, nDip, nGrad, nqInt, nQQ, nRoots
real(kind=wp) :: Columbus_Energy(2), E0, Temp, Value_l
logical(kind=iwp) :: Exist_2, Found, Is_Roots_Set, lMMGrd
real(kind=wp), allocatable :: DMs(:,:), MMGrd(:,:), Tmp(:)

!                                                                      *
!***********************************************************************
!                                                                      *
! Dummy-add the TS and saddle constraints, so that the arrays have
! a large enough size, these constraints will be actually added or
! not later on.

call Merge_Constraints('UDC','TSC','purge.Udc',mLambda,iDummy)
call Merge_Constraints('purge.Udc','UDC.Saddle','purge.Udc',mLambda,iDummy)
!                                                                      *
!***********************************************************************
!                                                                      *
! Get some basic information from the runfile.

call Get_Slapaf(iter,MaxItr,mTROld,lOld_Implicit,size(Coor,2),mLambda)
!                                                                      *
!***********************************************************************
!                                                                      *
! Save coordinates and gradients from this iteration onto the list.

Cx(:,:,iter) = Coor(:,:)
Gx(:,:,iter) = Grd(:,:)

if (iter > 1) then
  Temp = Zero
  do i=1,size(Coor,2)
    do j=1,3
      Temp = max(Temp,abs(Cx(j,i,iter)-Cx(j,i,iter-1)))
    end do
    if (Temp > Zero) exit
  end do
  if (Temp == Zero) then
    call WarningMessage(2,'Error in Init2')
    write(u6,*)
    write(u6,*) '****************** ERROR *********************'
    write(u6,*) 'Coordinates did not change!'
    write(u6,*) 'Maybe SEWARD is not inside the loop?'
    write(u6,*) '**********************************************'
    call Quit_OnUserError()
  end if
end if

! In case of a QM/MM geometry optimization, all the old MM gradients are
! replaced by the new one (both gradients are stored on the Runfile).

call Qpg_dArray('MM Grad',lMMGrd,nData)
lMMGrd = .false.
if (lMMGrd) then
  call mma_allocate(MMGrd,3*size(Coor,2),2,Label='MMGrd')
  call Get_dArray('MM Grad',MMGrd,3*size(Coor,2)*2)
  do i_N=1,iter-1
    write(u6,*) 'Grad at iteration :',i_N
    call RecPrt('Old:',' ',Gx(:,:,i_N),3,size(Coor,2))
    Gx(:,:,i_N) = Gx(:,:,i_N)+reshape(MMGrd(:,2)-MMGrd(:,1),[3,size(Coor,2)])
    call RecPrt('New:',' ',Gx(:,:,i_N),3,size(Coor,2))
  end do
  call mma_deallocate(MMGrd)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Pick up the reference structure from input or the run file.

if (Iter == 1) then

  ! Check if reference structure is defined by Gateway/Seward
  ! for the Saddle approach to find a TS.

  call qpg_dArray('Ref_Geom',Found,nData)
  if (Found) then
    if (.not. allocated(RefGeo)) call mma_allocate(RefGeo,3,size(Coor,2),Label='RefGeo')
    call Get_dArray('Ref_Geom',RefGeo,3*size(Coor,2))
  else

    ! Not defined: default reference structure to the starting structure.

    if (.not. allocated(RefGeo)) call mma_allocate(RefGeo,3,size(Coor,2),Label='RefGeo')
    RefGeo(:,:) = Cx(:,:,1)
    call Put_dArray('Ref_Geom',RefGeo,3*size(Coor,2))
  end if
else

  ! Pick up the reference structure.

  if (.not. allocated(RefGeo)) call mma_allocate(RefGeo,3,size(Coor,2),Label='RefGeo')
  call Get_dArray('Ref_Geom',RefGeo,3*size(Coor,2))
end if

! Align the reference structure to the current one, otherwise
! measuring distances does not make much sense

! (disabled for the moment, moving the reference affects the
!  computation of some vectors for MEP)
! If (iter > 1) Call Align(RefGeo,Coor,SIZE(Coor,2))
! Call RecPrt('Ref_Geom',' ',RefGeo,3,SIZE(Coor,2))
!                                                                      *
!***********************************************************************
!                                                                      *
! Get the energy of the last iteration

! Check if we are running in C&M mode.

call Get_iScalar('Columbus',columbus)

if (Columbus == 1) then
  call Get_dArray('MR-CISD energy',Columbus_Energy,2)
  Energy(iter) = Columbus_Energy(1)
else
  Is_Roots_Set = .false.
  call Qpg_iScalar('Number of roots',Is_Roots_Set)
  nRoots = 1
  if (Is_Roots_Set) call Get_iScalar('Number of roots',nRoots)
  !write(u6,*) 'Runfile'
  !write(u6,*) 'nRoots=',nRoots
  if (nRoots /= 1) then
    call Get_iScalar('NumGradRoot',iRoot)
    !write(u6,*) 'iRoot=',iRoot
    call mma_allocate(Tmp,nRoots,Label='Tmp')
    call Get_dArray('Last energies',Tmp,nRoots)
    !call RecPrt('Last Energies',' ',Tmp,1,nRoots)
    Energy(iter) = Tmp(iRoot)
    call mma_deallocate(Tmp)
  else
    call Get_dScalar('Last energy',Energy(iter))
  end if
end if
!                                                                      *
!*********** columbus interface ****************************************
!                                                                      *
! The dipole moment is required for numerical differentiation.
! Currently it is not written to the RUNFILE if we are in C&M mode.

if (columbus == 1) then
else
  Is_Roots_Set = .false.
  call Qpg_iScalar('Number of roots',Is_Roots_Set)
  nRoots = 1
  if (Is_Roots_Set) call Get_iScalar('Number of roots',nRoots)
  if (nRoots /= 1) then
    call Get_iScalar('NumGradRoot',iRoot)
    !write(u6,*) 'iRoot=',iRoot
    call mma_allocate(DMs,3,nRoots,Label='DMs')
    DMs(:,:) = Zero
    call Qpg_dArray('Last Dipole Moments',Found,nDip)
    if (Found .and. (nDip == 3*nRoots)) call Get_dArray('Last Dipole Moments',DMs,3*nRoots)
    DipM(:,iter) = DMs(:,iRoot)
    call mma_deallocate(DMs)
  else
    call Qpg_dArray('Dipole moment',Found,nDip)
    if (Found .and. (nDip == 3)) then
      call Get_dArray('Dipole moment',DipM(:,iter),3)
    else
      DipM(:,iter) = Zero
    end if
  end if
  !call RecPrt('Dipole Moment',' ',DipM(:,iter),1,3)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Test if this is a case of an intersection calculation. Depending
! on if we are running the calculation in M or C&M mode this is
! done in a bit different way.

if (Columbus == 1) then

  ! C&M mode

  call Get_iScalar('ColGradMode',iMode)
  if ((iMode == 2) .or. (iMode == 3)) then

    E0 = Columbus_Energy(2)
    Energy0(iter) = E0

    call qpg_dArray('Grad State2',Found,Length)
    if ((.not. Found) .or. (Length == 0)) call SysAbendmsg('Get_Molecule','Did not find:','Grad State2')
    call Get_dArray('Grad State2',Gx0(1,1,iter),Length)
    Gx0(:,:,iter) = -Gx0(:,:,iter)
    nSet = 2

  end if
  if (iMode == 3) call Get_dArray('NADC',NAC(:,:,iter),Length)

else

  ! M mode

  call f_Inquire('RUNFILE2',Exist_2)
  if (Exist_2) then
    call NameRun('RUNFILE2')

    Is_Roots_Set = .false.
    call Qpg_iScalar('Number of roots',Is_Roots_Set)
    nRoots = 1
    if (Is_Roots_Set) call Get_iScalar('Number of roots',nRoots)

    !write(u6,*) 'Runfile2'
    !write(u6,*) 'nRoots=',nRoots
    if (nRoots /= 1) then
      call Get_iScalar('NumGradRoot',iRoot)
      !write(u6,*) 'iRoot=',iRoot
      call mma_allocate(Tmp,nRoots,Label='Tmp')
      call Get_dArray('Last energies',Tmp,nRoots)
      !call RecPrt('Last Energies',' ',Tmp,1,nRoots)
      E0 = Tmp(iRoot)
      call mma_deallocate(Tmp)
    else
      call Get_dScalar('Last energy',E0)
    end if
    Energy0(iter) = E0

    nGrad = 3*size(Coor,2)
    call Get_dArray_chk('GRAD',Gx0(:,:,iter),nGrad)
    Gx0(:,:,iter) = -Gx0(:,:,iter)
    nSet = 2

    call NameRun('#Pop')
    TwoRunFiles = .true.
  end if
end if
if (NADC) nSet = 3
!                                                                      *
!***********************************************************************
!                                                                      *
! Pick up information from previous iterations
!                                                                      *
!***********************************************************************
!                                                                      *
if (iter /= 1) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call qpg_dArray('qInt',Found,nqInt)
  if (Found) then
    nQQ = nqInt/MaxItr
    call mma_allocate(qInt,nQQ,MaxItr,Label='qInt')
    call mma_allocate(dqInt,nQQ,MaxItr,Label='dqInt')
    call Get_dArray('qInt',qInt,nQQ*MaxItr)
    call Get_dArray('dqInt',dqInt,nQQ*MaxItr)
    if (nSet > 1) call mma_allocate(dqInt_Aux,nQQ,MaxItr,nSet-1,Label='dqInt_Aux')
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
Value_l = 20.0_wp
call Qpg_dScalar('Value_l',Found)
if (.not. Found) call Put_dScalar('Value_l',Value_l)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Init2
