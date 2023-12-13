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
! Copyright (C) 2023, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine Prepare_Kriging(Model_E,Model_G,nData,nDim,iFirst)
! Prepare energy and gradients in the right format for Setup_Kriging

use Slapaf_Info, only: dqInt, dqInt_Aux, Energy, Energy0, NADC
use Kriging_mod, only: Model_Type, nSet
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Half, Pi
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nData, nDim, iFirst
real(kind=wp), intent(out) :: Model_E(nData,nSet), Model_G(nDim,nData,nSet)
integer(kind=iwp) :: i, iLast, iRef
real(kind=wp) :: c1, c2, c3, Diff, Omega, Phi
real(kind=wp), allocatable :: Aux(:), BP(:,:), M(:,:), Ref(:,:)
real(kind=wp), external :: DDot_

iLast = iFirst+nData-1

! Trivial case, single surface: just copy energy and gradient

Model_E(:,1) = Energy(iFirst:iLast)
Model_G(:,:,1) = -dqInt(:,iFirst:iLast)

! Get data for additional surfaces

do i=2,nSet
  if (i == 2) then
    Model_E(:,i) = Energy0(iFirst:iLast)
  else
    Model_E(:,i) = Zero
  end if
  Model_G(:,:,i) = -dqInt_Aux(:,iFirst:iLast,i-1)
end do

! For more than one surface, unfold the stored energies & gradients

if (nSet >= 2) then
  call mma_allocate(Aux,nDim,Label='Aux')

  ! Diabatize the surfaces
  ! (this will have to be undone in kriging_update)

  ! From the initially stored values (EA, EB: adiabatic energies)
  !   E1 = (EA+EB)/2   G1 = (gA+gB)/2
  !   E2 = (EB-EA)     G2 = (gB-gA)
  !   E3 = 0           G3 = h

  if ((nSet > 2) .and. NADC) then

    ! To the model values (alpha, beta: diabatic energies,
    !                      gamma: coupling)
    !   E1 = alpha       G1 = g_alpha
    !   E2 = beta        G2 = g_beta
    !   E3 = gamma       G3 = g_gamma

    ! Set model type

    call mma_allocate(Model_Type,nSet,Label='Model_Type')
    Model_Type(:) = 1
    Model_Type(3) = 2 ! Coupling surface tends to zero

    ! It's easier to with EDiff/2

    Model_E(:,2) = Half*Model_E(:,2)
    Model_G(:,:,2) = Half*Model_G(:,:,2)

    ! Reference branching plane, the (g h) matrix at the latest iter.

    call mma_allocate(Ref,nDim,2,Label='Ref')
    call mma_allocate(BP,nDim,2,Label='BP')
    iRef = nData
    Ref(:,1) = Model_G(:,iRef,2)
    Ref(:,2) = Model_G(:,iRef,3)

    ! Pseudoinverse of the reference (g h)
    ! (transposed, and ignoring a constant factor)

    call mma_allocate(M,nDim,2,Label='M')
    c1 = DDot_(nDim,Ref(:,2),1,Ref(:,2),1)
    c2 = DDot_(nDim,Ref(:,1),1,Ref(:,2),1)
    c3 = DDot_(nDim,Ref(:,1),1,Ref(:,1),1)
    M(:,1) = c1*Ref(:,1)-c2*Ref(:,2)
    M(:,2) = c3*Ref(:,2)-c2*Ref(:,1)

    ! Transform all iterations

    do i=1,nData

      ! Get the rotation angle,
      ! assuming it's zero at the latest iteration

      if (i == iRef) then
        Omega = Zero
      else

        ! Match the branching planes

        BP(:,1) = Model_G(:,i,2)
        BP(:,2) = Model_G(:,i,3)
        call Rotate_BP(BP,Ref,nDim,2,Phi)

        ! Angles from the g and h vectors
        !   R = (g_0 h_0)^+ (g h)
        !   c1 = omega_g = atan(R(2,1)/R(1,1))
        !   c2 = omega_h = atan(-R(1,2)/R(2,2))

        c1 = atan2(DDot_(nDim,M(:,2),1,BP(:,1),1),DDot_(nDim,M(:,1),1,BP(:,1),1))
        c2 = atan2(-DDot_(nDim,M(:,1),1,BP(:,2),1),DDot_(nDim,M(:,2),1,BP(:,2),1))

        ! Average the angles, taking the periodicity into account.
        ! Reverse the h vector if that gives a smaller difference

        Diff = modulo(c2-c1+Pi,Two*Pi)-Pi
        if (abs(Diff) > Half*Pi) then
          Model_G(:,i,3) = -Model_G(:,i,3)
          ! c2 -> c2+Pi
          Diff = modulo(c2-c1,Two*Pi)-Pi
        end if
        Omega = c1+Half*Diff
      end if

      ! Once the rotation angle is known, we can obtain the
      ! diabatic surfaces (two energies and coupling)

      c1 = sin(Omega)
      c2 = cos(Omega)
      Diff = Model_E(i,2)
      Model_E(i,3) = c1*Diff
      Model_E(i,2) = Model_E(i,1)+c2*Diff
      Model_E(i,1) = Model_E(i,1)-c2*Diff
      Aux(:) = Model_G(:,i,2)
      Model_G(:,i,2) = c2*Aux-c1*Model_G(:,i,3)
      Model_G(:,i,3) = c1*Aux+c2*Model_G(:,i,3)
      Aux(:) = Model_G(:,i,2)
      Model_G(:,i,2) = Model_G(:,i,1)+Aux
      Model_G(:,i,1) = Model_G(:,i,1)-Aux
    end do
    call mma_deallocate(M)
    call mma_deallocate(Ref)
    call mma_deallocate(BP)
  else

    ! Or just unfold to EA,EB and gA,gB

    do i=1,nData
      c1 = Model_E(i,1)
      Model_E(i,1) = c1+Half*Model_E(i,2)
      Model_E(i,2) = c1-Half*Model_E(i,2)
      Aux(:) = Model_G(:,i,1)
      Model_G(:,i,1) = Aux+Half*Model_G(:,i,2)
      Model_G(:,i,2) = Aux-Half*Model_G(:,i,2)
    end do
  end if
  call mma_deallocate(Aux)
end if

end subroutine Prepare_Kriging
