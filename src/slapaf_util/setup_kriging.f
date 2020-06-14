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
* Copyright (C) 2020, Roland Lindh                                     *
************************************************************************
      Subroutine SetUp_Kriging(nRaw,nInter,qInt,Grad,Energy,Hessian_HMF)
      Use kriging_mod, only: blavAI, set_l
      Use Limbo
      Implicit None
#include "stdalloc.fh"
#include "real.fh"
      Integer nRaw, nInter,i,iInter,jInter,ij
      Real*8 qInt(nInter,nRaw), Grad(nInter,nRaw), Energy(nRaw),
     &       Hessian_HMF(nInter,nInter)
      Real*8 Value_l
      Real*8, Allocatable:: Array_l(:), HTri(:), Hessian(:,:),
     &                      qInt_s(:,:), Grad_s(:,:)
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_allocate(U,nInter,nInter,Label='U')
*
      U(:,:)=Zero
      Forall (i=1:nInter) U(i,i)=One
*
      Call mma_allocate(Hessian,nInter,nInter,Label='Hessian')
      Call mma_allocate(HTri,nInter*(nInter+1)/2,Label='HTri')
*
      Do iInter = 1, nInter
         Do jInter = 1, iInter
            ij = iInter*(iInter-1)/2 + jInter
            HTri(ij)=Hessian_HMF(iInter,jInter)
         End Do
      End Do
*     Call TriPrt('HTri(raw)',' ',HTri,nInter)
      Call NIDiag_new(HTri,U,nInter,nInter,0)
*     U(:,:)=Zero
*     Forall (i=1:nInter) U(i,i)=One
*     Call TriPrt('HTri',' ',HTri,nInter)
      Hessian(:,:) = Zero
      Forall (i=1:nInter) Hessian(i,i)=HTri(i*(i+1)/2)
*
      Call mma_deallocate(HTri)
*                                                                      *
************************************************************************
*                                                                      *
*     Select between setting all ls to a single value or go in
*     multiple l-value mode in which the l-value is set such that
*     the kriging hessian reproduce the diagonal value of the HMF
*     Hessian of the current structure.
*
*#define _DEBUG_
#ifdef _DEBUG_
      Call RecPrt('Setup_kriging: Energy',' ',Energy,1,nRaw)
      Call RecPrt('Setup_kriging: qInt',' ',qInt,nInter,nRaw)
      Call RecPrt('Setup_kriging: Grad',' ',Grad,nInter,nRaw)
#endif
      Call mma_Allocate(Array_l,nInter,Label='Array_l')
      If (Set_l) Then
         Call Get_dScalar('Value_l',Value_l)
         Array_l(:)=Value_l
      Else
         Call Set_l_Array(Array_l,nInter,blavAI,Hessian)
      End If
      Call mma_deallocate(Hessian)
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_Allocate(qInt_s,nInter,nRaw,Label="qInt_s")
      Call mma_Allocate(Grad_s,nInter,nRaw,Label="Grad_s")
*
*     Transform to the basis which diagonalizes the HMF Hessian.
*
      Call Trans_K(U,qInt,qInt_s,nInter,nRaw)
      Call Trans_K(U,Grad,Grad_s,nInter,nRaw)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
      Call RecPrt('Setup_kriging: qInt_s',' ',qInt_s,nInter,nRaw)
      Call RecPrt('Setup_kriging: Grad_s',' ',Grad_s,nInter,nRaw)
#endif
      Call Start_Kriging(nRaw,nInter,qInt_s,Grad_s,Energy)
*
      Call mma_deAllocate(qInt_s)
      Call mma_deAllocate(Grad_s)
*                                                                      *
************************************************************************
*                                                                      *
*     Pass the l-values to the GEK routine. This will initiate the
*     computation of the covariance matrix, and solve related GEK
*     equations.
*
      Call Set_l_Kriging(Array_l,nInter)
      Call mma_deAllocate(Array_l)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End Subroutine Setup_Kriging
