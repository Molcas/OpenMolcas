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
Subroutine Kriging_Update(nQQ,iter,qInt,E_Disp)
Use Slapaf_Info, only: Energy, dqInt, Energy0, dqInt_Aux, Gx, Gx0, NAC, BMx_kriging, Degen
Use Slapaf_Parameters, only: Curvilinear, NADC
Use Kriging_Mod, only: nSet
Implicit None
Integer nQQ, iter
Real*8  qInt(nQQ), E_Disp

#include "real.fh"
#include "stdalloc.fh"
Integer :: nAtoms
Real*8 :: Diff, Omega
Real*8, Allocatable :: Aux(:,:), Demp(:), Temp(:), vAux(:)

Call mma_allocate(Temp,nSet,Label='Temp')
Call mma_allocate(Demp,nSet,Label='Demp')
Call mma_allocate(Aux,nQQ,nSet,Label='Aux')

!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
Call RecPrt('Kriging_Update: qInt',' ',qInt,nQQ,1)
#endif

Call Energy_Kriging_layer(qInt,Temp,nQQ)

Call Dispersion_Kriging_Layer(qInt,Demp,nQQ)

Call Gradient_Kriging_layer(qInt,Aux,nQQ)

#ifdef _DEBUGPRINT_
Call RecPrt('Kriging_Update: Temp',' ',Temp,1,nSet)
Call RecPrt('Kriging_Update: Demp',' ',Demp,1,nSet)
Call RecPrt('Kriging_Update: Aux',' ',Aux,nQQ,nSet)
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
! Refold or undo the diabatization

If (nSet >= 2) Then
  Call mma_allocate(vAux,nQQ,Label='vAux')
  If ((nSet > 2) .And. NADC) Then
    Diff = Half*(Temp(2)-Temp(1))
    Omega = ATan2(Temp(3),Diff)
    ! energies and dispersion
    Temp(1) = Half*(Temp(1)+Temp(2))
    Demp(1) = Half*(Demp(1)+Demp(2))
    Temp(2) = Two*Sqrt(Diff**2+Temp(3)**2)
    Temp(3) = Zero
    ! gradients
    vAux(:) = Half*(Aux(:,2)-Aux(:,1))
    Aux(:,1) = Half*(Aux(:,1)+Aux(:,2))
    Aux(:,2) = Two*(Cos(Omega)*vAux+Sin(Omega)*Aux(:,3))
    Aux(:,3) = -Sin(Omega)*vAux+Cos(Omega)*Aux(:,3)
  Else
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
  End If
  Call mma_deallocate(vAux)
End If

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

  nAtoms = Size(Gx0,2)

  call DGEMM_('N','N',                                  &
              3*nAtoms,1,nQQ,                           &
              One,BMx_kriging,3*nAtoms,                 &
                  dqInt(:,iter),nQQ,                    &
             Zero,Gx(:,:,iter),3*nAtoms)
  call DGEMM_('N','N',                                  &
              3*nAtoms,1,nQQ,                           &
              One,BMx_kriging,3*nAtoms,                 &
                  dqInt_Aux(:,iter,1),nQQ,              &
             Zero,Gx0(:,:,iter),3*nAtoms)

  ! Modify with degeneracy factors.

  if (Curvilinear) Gx(:,:,iter) = Gx(:,:,iter)/Degen(:,:)
  if (Curvilinear) Gx0(:,:,iter) = Gx0(:,:,iter)/Degen(:,:)

  if (nSet > 2) then
    dqInt_Aux(:,iter,2) = Aux(:,3)
    call DGEMM_('N','N',                                &
                3*nAtoms,1,nQQ,                         &
                One,BMx_kriging,3*nAtoms,               &
                    dqInt_Aux(:,iter,2),nQQ,            &
               Zero,NAC(:,:,iter),3*nAtoms)
    if (Curvilinear) NAC(:,:,iter) = NAC(:,:,iter)/Degen(:,:)
  end if

end if

Call mma_deallocate(Temp)
Call mma_deallocate(Demp)
Call mma_deallocate(Aux)

End Subroutine Kriging_Update
