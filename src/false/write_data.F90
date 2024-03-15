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
! Copyright (C) 2020, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine Write_Data()

use stdalloc, only: mma_Allocate, mma_Deallocate
use False_Global, only: Mode, Will_Print
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nAtoms, nCart, nHess, nRoots, nRelax, mRoots, i, j, k, l, LU, EOF
real(kind=wp), allocatable :: Energies(:), Gradient(:), NAC(:), Hessian(:), Dipoles(:,:), &
                              mEnergies(:), mGradient(:), mHessian(:), mNAC(:), mDipoles(:,:)
integer(kind=iwp), allocatable :: not_grad(:), not_nac(:,:)
character(len=16) :: line
logical(kind=iwp) :: Found, Too_Late
integer(kind=iwp), external :: IsFreeUnit, Read_Grad

call Get_nAtoms_All(nAtoms)
nCart = 3*nAtoms
nHess = nCart*(nCart+1)/2
if (Will_Print) write(u6,*)

! read data from interface output file
! write to RUNFILE or GRADS as we read it
LU = IsFreeUnit(11)
call Molcas_Open(LU,'OUTPUT')
Too_Late = .false.
nRoots = 0
do
  read(LU,100,iostat=EOF) line
  if (EOF < 0) exit
  call UpCase(line)
  select case (line)
    case ('[ROOTS]')
      read(LU,*) nRoots
      nRelax = nRoots
      call Check_nRoots(nRoots)
      if (Will_Print) write(u6,200) nRoots
      if (Mode == 'ADD') then
        call Get_iScalar('Number of roots',mRoots)
        call Check_nRoots(mRoots)
      else
        call Put_iScalar('Number of roots',nRoots)
        call Put_iScalar('Relax CASSCF root',nRelax)
      end if
      call mma_Allocate(not_grad,nRoots)
      call mma_Allocate(not_nac,nRoots,nRoots)
      not_grad(:) = 1
      not_nac(:,:) = 1
    case ('[RELAX ROOT]')
      call Check_Too_Late()
      call Check_nRoots()
      read(LU,*) nRelax
      call Check_Root_Number(nRelax)
      if (Will_Print) write(u6,201) nRelax
      call Put_iScalar('Relax CASSCF root',nRelax)
    case ('[ENERGIES]')
      call Check_nRoots()
      Too_Late = .true.
      call mma_Allocate(Energies,nRoots,label='Energies')
      read(LU,*) Energies(:)
      if (Will_Print) call RecPrt('Root energies','',Energies,nRoots,1)
      if (Mode == 'ADD') then
        call mma_Allocate(mEnergies,mRoots,label='mEnergies')
        call Get_dArray('Last energies',mEnergies,mRoots)
        if (nRoots > 1) then
          mEnergies(:) = mEnergies(:)+Energies(1)
        else
          mEnergies(:) = mEnergies(:)+Energies(:)
        end if
        call Store_Energies(mRoots,mEnergies,nRelax)
        call mma_Deallocate(mEnergies)
      else
        call Put_cArray('Relax Method','EXTERNAL',8)
        call Store_Energies(nRoots,Energies,nRelax)
      end if
      call mma_Deallocate(Energies)
    case ('[GRADIENT]')
      call Check_nRoots()
      read(LU,*) i
      call Check_Root_Number(i)
      call mma_Allocate(Gradient,nCart,label='Gradient')
      read(LU,*) Gradient(:)
      if (Will_Print) write(u6,202) i
      not_grad(i) = 0
      if (Mode == 'ADD') then
        call mma_Allocate(mGradient,nCart,label='mGradient')
        if (nRoots == 1) then
          ! modify all existing gradients
          do j=1,mRoots
            if (Read_Grad(mGradient,nCart,j,0,0) == 1) then
              mGradient(:) = mGradient(:)+Gradient(:)
              call Store_Grad(mGradient,nCart,j,0,0)
            end if
          end do
        else
          ! modify only this gradient
          if (Read_Grad(mGradient,nCart,i,0,0) == 1) then
            mGradient(:) = mGradient(:)+Gradient(:)
            call Store_Grad(mGradient,nCart,i,0,0)
          end if
        end if
        ! modify the runfile gradient
        call Qpg_dArray('GRAD',Found,i)
        if (Found .and. (i == nCart)) then
          call Get_dArray('GRAD',mGradient,nCart)
          mGradient(:) = mGradient(:)+Gradient(:)
          call Put_dArray('GRAD',mGradient,nCart)
        end if
        call mma_Deallocate(mGradient)
      else
        call Store_Grad(Gradient,nCart,i,0,0)
        call Qpg_dArray('GRAD',Found,i)
        if (Found .and. (i == nCart)) call Put_dArray('GRAD',Gradient,nCart)
      end if
      call mma_Deallocate(Gradient)
    case ('[NAC]')
      call Check_nRoots()
      read(LU,*) i,j
      call Check_Root_Pair(i,j)
      call mma_Allocate(NAC,nCart,label='NAC')
      read(LU,*) NAC(:)
      if (Will_Print) write(u6,203) i,j
      not_nac(i,j) = 0
      not_nac(j,i) = 0
      if (Mode == 'ADD') then
        call mma_Allocate(mNAC,nCart,label='mNAC')
        if (nRoots == 1) then
          ! modify all existing couplings
          do k=1,mRoots
            do l=k+1,mRoots
              if (Read_Grad(mNAC,nCart,0,k,l) == 1) then
                mNAC(:) = mNAC(:)+NAC(:)
                call Store_Grad(mNAC,nCart,0,k,l)
              end if
            end do
          end do
        else
          ! modify only this coupling
          if (Read_Grad(mNAC,nCart,0,i,j) == 1) then
            mNAC(:) = mNAC(:)+NAC(:)
            call Store_Grad(mNAC,nCart,0,i,j)
          end if
        end if
        call mma_Deallocate(mNAC)
      else
        call Store_Grad(NAC,nCart,0,i,j)
      end if
      call mma_Deallocate(NAC)
    case ('[HESSIAN]')
      call Check_nRoots()
      Too_Late = .true.
      read(LU,*) i
      call Check_Root_Number(i)
      if (i == nRelax) then
        call mma_Allocate(Hessian,nHess,label='Hessian')
        read(LU,*) Hessian(:)
        if (Mode == 'ADD') then
          call mma_Allocate(mHessian,nHess,label='mHessian')
          call Qpg_dArray('Analytic Hessian',Found,j)
          if (Found) then
            call Get_dArray_chk('Analytic Hessian',mHessian,nHess)
            mHessian(:) = mHessian(:)+Hessian(:)
            call Put_AnalHess(mHessian,nHess)
          end if
          call mma_Deallocate(mHessian)
        else
          call Put_AnalHess(Hessian,nHess)
        end if
        call mma_Deallocate(Hessian)
      end if
      if (Will_Print) write(u6,204) i
    case ('[DIPOLES]')
      call Check_nRoots()
      call mma_Allocate(Dipoles,3,nRoots,label='Dipoles')
      do i=1,nRoots
        read(LU,*) Dipoles(:,i)
      end do
      if (Will_Print) write(u6,*) 'Found dipole moments'
      if (Mode == 'ADD') then
        call mma_Allocate(mDipoles,3,mRoots,label='mDipoles')
        call Get_dArray('Last dipole momentes',mDipoles,3*mRoots)
        if (nRoots > 1) then
          do j=1,mRoots
            mDipoles(:,j) = mDipoles(:,j)+Dipoles(:,1)
          end do
        else
          mDipoles(:,:) = mDipoles(:,:)+Dipoles(:,:)
        end if
        call Put_dArray('Last dipole moments',mDipoles,3*nRoots)
        call mma_Deallocate(mDipoles)
      else
        call Put_dArray('Last dipole moments',Dipoles,3*nRoots)
      end if
      call mma_Deallocate(Dipoles)
    case default
  end select
end do
close(LU)

! mark non-computable gradients
if (allocated(not_grad)) then
  if (Mode == 'REPLACE') then
    do j=1,nRoots
      if (not_grad(j) /= 0) call Store_Not_Grad(j,0,0)
      do i=j+1,nRoots
        if (not_nac(i,j) /= 0) call Store_Not_Grad(0,i,j)
      end do
    end do
  end if
  call mma_Deallocate(not_grad)
  call mma_Deallocate(not_nac)
end if

return

100 format(a)
200 format('Found data for ',i3,' roots')
201 format('Relaxing on root ',i3)
202 format('Found gradient for root ',i3)
203 format('Found coupling vector for roots ',i3,' and ',i3)
204 format('Found Hessian for root ',i3)

contains

subroutine Check_nRoots(n)
  integer(kind=iwp), intent(in), optional :: n
  if (nRoots < 1) then
    if (present(n)) then
      call WarningMessage(2,'The number of roots must be positive.')
    else
      call WarningMessage(2,'[ROOTS] should be defined first.')
    end if
    call AbEnd()
  else if ((Mode == 'ADD') .and. present(n)) then
    ! if nRoots == 1 and mRoots /= 1, the same values will be added to all roots
    if ((nRoots /= 1) .and. (nRoots /= n)) then
      call WarningMessage(2,'The number of roots does not agree with the runfile')
      call AbEnd()
    end if
  end if
end subroutine Check_nRoots

subroutine Check_Too_Late()
  if (Too_Late) then
    call WarningMessage(2,'[RELAX] should have been given earlier.')
    call AbEnd()
  end if
end subroutine Check_Too_Late

subroutine Check_Root_Number(n)
  integer(kind=iwp), intent(in) :: n
  character(len=6) :: a, b
  call Check_nRoots()
  if ((n < 1) .or. (n > nRoots)) then
    write(a,'(I6)') n
    write(b,'(I6)') nRoots
    a = adjustl(a)
    b = adjustl(b)
    call WarningMessage(2,'Root number '//trim(a)//' must be between 1 and '//trim(b)//'.')
    call AbEnd()
  end if
end subroutine Check_Root_Number

subroutine Check_Root_Pair(n,m)
  integer(kind=iwp), intent(in) :: n, m
  character(len=6) :: a
  call Check_Root_Number(n)
  call Check_Root_Number(m)
  if (n == m) then
    write(a,'(I6)') n
    a = adjustl(a)
    call WarningMessage(2,'Roots in pair '//trim(a)//' '//trim(a)//' cannot be equal.')
    call AbEnd()
  end if
end subroutine Check_Root_Pair

end subroutine Write_Data
