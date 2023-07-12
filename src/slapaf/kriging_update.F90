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
! Copyright (C) 2021, Roland Lindh                                     *
!***********************************************************************

subroutine Kriging_Update(nQQ,iter,qInt,E_Disp)

use Slapaf_Info, only: BMx_kriging, Curvilinear, Degen, dqInt, dqInt_Aux, Energy, Energy0, Gx, Gx0, NAC, NADC
use Kriging_Mod, only: nSet
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nQQ, iter
real(kind=wp), intent(in) :: qInt(nQQ)
real(kind=wp), intent(out) :: E_Disp
integer(kind=iwp) :: nAtoms
real(kind=wp) :: Diff, Omega
real(kind=wp), allocatable :: Aux(:,:), Demp(:), Temp(:), vAux(:)

call mma_allocate(Temp,nSet,Label='Temp')
call mma_allocate(Demp,nSet,Label='Demp')
call mma_allocate(Aux,nQQ,nSet,Label='Aux')

!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
call RecPrt('Kriging_Update: qInt',' ',qInt,nQQ,1)
#endif

call Energy_Kriging_layer(qInt,Temp,nQQ)

call Dispersion_Kriging_Layer(qInt,Demp,nQQ)

call Gradient_Kriging_layer(qInt,Aux,nQQ)

#ifdef _DEBUGPRINT_
call RecPrt('Kriging_Update: Temp',' ',Temp,1,nSet)
call RecPrt('Kriging_Update: Demp',' ',Demp,1,nSet)
call RecPrt('Kriging_Update: Aux',' ',Aux,nQQ,nSet)
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
! Refold or undo the diabatization

if (nSet >= 2) then
  call mma_allocate(vAux,nQQ,Label='vAux')
  if ((nSet > 2) .and. NADC) then
    Diff = Half*(Temp(2)-Temp(1))
    Omega = atan2(Temp(3),Diff)
    ! energies and dispersion
    Temp(1) = Half*(Temp(1)+Temp(2))
    Demp(1) = Half*(Demp(1)+Demp(2))
    Temp(2) = Two*sqrt(Diff**2+Temp(3)**2)
    Temp(3) = Zero
    ! gradients
    vAux(:) = Half*(Aux(:,2)-Aux(:,1))
    Aux(:,1) = Half*(Aux(:,1)+Aux(:,2))
    Aux(:,2) = Two*(cos(Omega)*vAux+sin(Omega)*Aux(:,3))
    Aux(:,3) = -sin(Omega)*vAux+cos(Omega)*Aux(:,3)
  else
    ! energies and dispersion
    Diff = Temp(1)-Temp(2)
    Temp(1) = Half*(Temp(1)+Temp(2))
    Temp(2) = Diff
    Demp(1) = Half*(Demp(1)+Demp(2))
    Demp(2) = Demp(1)
    ! gradients
    vAux(:) = Aux(:,1)-Aux(:,2)
    Aux(:,1) = Half*(Aux(:,1)+Aux(:,2))
    Aux(:,2) = vAux
  end if
  call mma_deallocate(vAux)
end if

!                                                                      *
!***********************************************************************
!                                                                      *

Energy(iter) = Temp(1)

E_Disp = Demp(1)

dqInt(:,iter) = -Aux(:,1)

!                                                                      *
!***********************************************************************
!                                                                      *
!  For the energy difference

if (nSet > 1) then

  Energy0(iter) = Temp(2)

  ! Right now we do not use the dispersion for the additional surfaces

  dqInt_Aux(:,iter,1) = -Aux(:,2)

  ! The computation of the gradients of the constraints is always done in
  ! Cartesian coordinates. The GEK, however, predicts them in internal coordinates.
  ! Thus, we need to transfrom the GEK predicted gradients to Cartesian and put
  ! them at the place where the code to compute the value and the Cartesian gradient
  ! can find them, that is, in Gx, Gx0

  ! dE/dx = dq/dx dE/dq

  nAtoms = size(Gx0,2)

  call DGEMM_('N','N',3*nAtoms,1,nQQ,One,BMx_kriging,3*nAtoms,dqInt(:,iter),nQQ,Zero,Gx(:,:,iter),3*nAtoms)
  call DGEMM_('N','N',3*nAtoms,1,nQQ,One,BMx_kriging,3*nAtoms,dqInt_Aux(:,iter,1),nQQ,Zero,Gx0(:,:,iter),3*nAtoms)

  ! Modify with degeneracy factors.

  if (Curvilinear) then
    Gx(:,:,iter) = Gx(:,:,iter)/Degen(:,:)
    Gx0(:,:,iter) = Gx0(:,:,iter)/Degen(:,:)
  end if

  if (nSet > 2) then
    dqInt_Aux(:,iter,2) = Aux(:,3)
    call DGEMM_('N','N',3*nAtoms,1,nQQ,One,BMx_kriging,3*nAtoms,dqInt_Aux(:,iter,2),nQQ,Zero,NAC(:,:,iter),3*nAtoms)
    if (Curvilinear) NAC(:,:,iter) = NAC(:,:,iter)/Degen(:,:)
  end if

end if

call mma_deallocate(Temp)
call mma_deallocate(Demp)
call mma_deallocate(Aux)

end subroutine Kriging_Update
