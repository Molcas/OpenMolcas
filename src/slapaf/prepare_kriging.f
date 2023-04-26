************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2023, Ignacio Fdez. Galvan                             *
************************************************************************

      Subroutine Prepare_Kriging(Model_E,Model_G,nData,nDim,iFirst)
*     Prepare energy and gradients in the right format for Setup_Kriging

      Use Kriging_mod, Only: Model_Type, nSet
      Use Slapaf_Info, Only: dqInt, dqInt_Aux, Energy, Energy0
      Use Slapaf_Parameters, Only: NADC
      Use stdalloc, Only: mma_allocate, mma_deallocate
      Use Constants, Only: Zero, Two, Half, Pi
      Use Definitions, Only: wp, iwp

      Implicit None
      Integer(kind=iwp), Intent(in) :: nData, nDim, iFirst
      Real(kind=wp), Intent(out) :: Model_E(nData,nSet),
     &                              Model_G(nDim,nData,nSet)
      Integer(kind=iwp) :: i, iLast, iRef
      Real(kind=wp) :: c1, c2, c3, Diff, Omega, Phi
      Real(kind=wp), Allocatable :: Aux(:), BP(:,:), M(:,:), Ref(:,:)
      Real(kind=wp), External :: DDot_

      iLast = iFirst+nData-1
*
*     Trivial case, single surface: just copy energy and gradient
*
      Model_E(:,1) = Energy(iFirst:iLast)
      Model_G(:,:,1) = -dqInt(:,iFirst:iLast)
*
*     Get data for additional surfaces
*
      Do i=2,nSet
        If (i == 2) Then
          Model_E(:,i) = Energy0(iFirst:iLast)
        Else
          Model_E(:,i) = Zero
        End If
        Model_G(:,:,i) = -dqInt_Aux(:,iFirst:iLast,i-1)
      End Do
*
*     For more than one surface, unfold the stored energies & gradients
*
      If (nSet >= 2) Then
        Call mma_allocate(Aux,nDim,Label='Aux')
*
*       Diabatize the surfaces
*       (this will have to be undone in kriging_update)
*
*       From the initially stored values (EA, EB: adiabatic energies)
*         E1 = (EA+EB)/2   G1 = (gA+gB)/2
*         E2 = (EB-EA)     G2 = (gB-gA)
*         E3 = 0           G3 = h
*
        If ((nSet > 2) .And. NADC) Then
*
*         To the model values (alpha, beta: diabatic energies,
*                              gamma: coupling)
*           E1 = alpha       G1 = g_alpha
*           E2 = beta        G2 = g_beta
*           E3 = gamma       G3 = g_gamma
*
*         Set model type
*
          Call mma_allocate(Model_Type,nSet,Label='Model_Type')
          Model_Type(:) = 1
          Model_Type(3) = 2 ! Coupling surface tends to zero
*
*         It's easier to with EDiff/2
*
          Model_E(:,2) = Half*Model_E(:,2)
          Model_G(:,:,2) = Half*Model_G(:,:,2)
*
*         Reference branching plane, the (g h) matrix at the latest iter.
*
          Call mma_allocate(Ref,nDim,2,Label='Ref')
          Call mma_allocate(BP,nDim,2,Label='BP')
          iRef = nData
          Ref(:,1) = Model_G(:,iRef,2)
          Ref(:,2) = Model_G(:,iRef,3)
*
*         Pseudoinverse of the reference (g h)
*         (transposed, and ignoring a constant factor)
*
          Call mma_allocate(M,nDim,2,Label='M')
          c1 = DDot_(nDim,Ref(:,2),1,Ref(:,2),1)
          c2 = DDot_(nDim,Ref(:,1),1,Ref(:,2),1)
          c3 = DDot_(nDim,Ref(:,1),1,Ref(:,1),1)
          M(:,1) = c1*Ref(:,1)-c2*Ref(:,2)
          M(:,2) = c3*Ref(:,2)-c2*Ref(:,1)
*
*         Transform all iterations
*
          Do i=1,nData
*
*           Get the rotation angle,
*           assuming it's zero at the latest iteration
*
            If (i == iRef) Then
              Omega = Zero
            Else
*
*             Match the branching planes
*
              BP(:,1) = Model_G(:,i,2)
              BP(:,2) = Model_G(:,i,3)
              Call Rotate_BP(BP,Ref,nDim,2,Phi)
*
*             Angles from the g and h vectors
*               R = (g_0 h_0)^+ (g h)
*               c1 = omega_g = atan(R(2,1)/R(1,1))
*               c2 = omega_h = atan(-R(1,2)/R(2,2))
*
              c1 = ATan2(DDot_(nDim,M(:,2),1,BP(:,1),1),
     &                   DDot_(nDim,M(:,1),1,BP(:,1),1))
              c2 = ATan2(-DDot_(nDim,M(:,1),1,BP(:,2),1),
     &                   DDot_(nDim,M(:,2),1,BP(:,2),1))
*
*             Average the angles, taking the periodicity into account.
*             Reverse the h vector if that gives a smaller difference
*
              Diff = Modulo(c2-c1+Pi,Two*Pi)-Pi
              If (Abs(Diff) > Half*Pi) Then
                Model_G(:,i,3) = -Model_G(:,i,3)
*               c2 -> c2+Pi
                Diff = Modulo(c2-c1,Two*Pi)-Pi
              End If
              Omega = c1+Half*Diff
            End If
*
*           Once the rotation angle is known, we can obtain the
*           diabatic surfaces (two energies and coupling)
*
            c1 = Sin(Omega)
            c2 = Cos(Omega)
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
          End Do
          Call mma_deallocate(M)
          Call mma_deallocate(Ref)
          Call mma_deallocate(BP)
        Else
*
*         Or just unfold to EA,EB and gA,gB
*
          Do i=1,nData
            c1 = Model_E(i,1)
            Model_E(i,1) = c1+Half*Model_E(i,2)
            Model_E(i,2) = c1-Half*Model_E(i,2)
            Aux(:) = Model_G(:,i,1)
            Model_G(:,i,1) = Aux+Half*Model_G(:,i,2)
            Model_G(:,i,2) = Aux-Half*Model_G(:,i,2)
          End Do
        End If
        Call mma_deallocate(Aux)
      End If

      End Subroutine Prepare_Kriging
