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
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

subroutine Setup_Kriging(nRaw,nInter,qInt,Grad,Energy,Hessian_HMF,HDiag)

use kriging_mod, only: blavAI, nSet, layer_U, set_l
! This will be in the same module
!use kriging_procedures, only: set_l_Array
use Index_Functions, only: iTri, nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nRaw, nInter
real(kind=wp), intent(in) :: qInt(nInter,nRaw), Grad(nInter,nRaw,nSet), Energy(nRaw,nSet)
real(kind=wp), intent(inout), optional :: Hessian_HMF(nInter,nInter), HDiag(nInter)
integer(kind=iwp) :: i, iInter, ij, jInter
real(kind=wp) :: Value_l
real(kind=wp), allocatable :: Array_l(:), dqInt_s(:,:,:), Energy_s(:,:), Hessian(:,:), HTri(:), qInt_s(:,:)

!                                                                      *
!***********************************************************************
!                                                                      *
if (present(Hessian_HMF)) then

  call mma_allocate(layer_U,nInter,nInter,Label='layer_U')

  call unitmat(layer_U,nInter)

  call mma_allocate(Hessian,nInter,nInter,Label='Hessian')
  call mma_allocate(HTri,nTri_Elem(nInter),Label='HTri')

  do iInter=1,nInter
    do jInter=1,iInter
      ij = iTri(iInter,jInter)
      HTri(ij) = Hessian_HMF(iInter,jInter)
    end do
  end do
  call NIDiag_new(HTri,layer_U,nInter,nInter)
  Hessian(:,:) = Zero
  do i=1,nInter
    Hessian(i,i) = HTri(iTri(i,i))
  end do

  call mma_deallocate(HTri)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Select between setting all ls to a single value or go in
  ! multiple l-value mode in which the l-value is set such that
  ! the kriging hessian reproduce the diagonal value of the HMF
  ! Hessian of the current structure.

!# define _DEBUGPRINT_
# ifdef _DEBUGPRINT_
  call RecPrt('Setup_kriging: qInt',' ',qInt,nInter,nRaw)
  do i=1,nSet
    write(u6,*) 'iSet=',i
    call RecPrt('Setup_kriging: Energy',' ',Energy(:,i),1,nRaw)
    call RecPrt('Setup_kriging: Grad',' ',Grad(:,:,i),nInter,nRaw)
  end do
# endif
  call mma_allocate(Array_l,nInter,Label='Array_l')
  if (Set_l) then
    call Get_dScalar('Value_l',Value_l)
    Array_l(:) = Value_l
  else
    call Set_l_Array(Array_l,nInter,blavAI,Hessian=Hessian)
  end if
  call mma_deallocate(Hessian)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call mma_allocate(qInt_s,nInter,nRaw,Label="qInt_s")
  call mma_allocate(dqInt_s,nInter,nRaw,nSet,Label="dqInt_s")
  call mma_allocate(Energy_s,nRaw,nSet,Label="Energy_s")

  ! Transform to the basis which diagonalizes the HMF Hessian.

  call Trans_K(qInt,qInt_s,nInter,nRaw)
  do i=1,nSet
    call Trans_K(Grad(:,:,i),dqInt_s(:,:,i),nInter,nRaw)
  end do
  Energy_s(:,:) = Energy
  !                                                                    *
  !*********************************************************************
  !                                                                    *
# ifdef _DEBUGPRINT_
  call RecPrt('Setup_kriging: qInt_s',' ',qInt_s,nInter,nRaw)
  do i=1,nSet
    write(u6,*) 'iSet=',i
    call RecPrt('Setup_kriging: Energy_s',' ',Energy_s(:,i),1,nRaw)
    call RecPrt('Setup_kriging: dqInt_s',' ',dqInt_s(:,:,i),nInter,nRaw)
  end do
# endif
  call Start_Kriging(nRaw,nInter,qInt_s,dqInt_s,Energy_s)

  call mma_deallocate(Energy_s)
  call mma_deallocate(dqInt_s)
  call mma_deallocate(qInt_s)

else

  call mma_allocate(Array_l,nInter,Label='Array_l')
  if (Set_l) then
    call Get_dScalar('Value_l',Value_l)
    Array_l(:) = Value_l
  else
    call Set_l_Array(Array_l,nInter,blavAI,HDiag=HDiag)
  end if

  call Start_Kriging(nRaw,nInter,qInt,Grad,Energy)

end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Pass the l-values to the GEK routine. This will initiate the
! computation of the covariance matrix, and solve related GEK
! equations.

call Set_l_Kriging(Array_l,nInter)
call mma_deallocate(Array_l)

!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Setup_Kriging
