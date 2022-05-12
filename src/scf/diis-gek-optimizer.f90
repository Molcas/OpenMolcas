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
! Copyright (C) 2022, Roland Lindh                                     *
!***********************************************************************
!#define _DEBUGPRINT_
Subroutine DIIS_GEK_Optimizer(dq,mOV)
!***********************************************************************
!                                                                      *
!     Object: Direct-inversion-in-the-iterative-subspace gradient-     *
!             enhanced kriging optimization.                           *
!                                                                      *
!     Author: Roland Lindh, Dept. of Chemistry -- BMC,                 *
!             University of Uppsala, SWEDEN                            *
!             May '22                                                  *
!***********************************************************************
use InfSO , only: iterso, Energy
use InfSCF, only: iter
use LnkLst, only: SCF_V, Init_LLs, LLx, LLGrad
use SCF_Arrays, only: HDiag
Implicit None
#include "real.fh"
#include "stdalloc.fh"
Integer, Intent(In):: mOV
Real*8,  Intent(Out):: dq(mOV)

Integer i, j, k, l, ipq, ipg, nDIIS, mDIIS, iFirst
Integer, External:: LstPtr
#ifdef _DEBUGPRINT_
Real*8, External::DDot_
#endif
Real*8, Allocatable:: q(:,:), g(:,:)
Real*8, Allocatable:: q_diis(:,:), g_diis(:,:), e_diis(:,:)
Real*8, Allocatable:: dq_diis(:)
Real*8, Allocatable:: H_Diis(:,:)
Real*8 :: gg
Character(Len=1) Step_Trunc, Step_Trunc_
Character(Len=6) UpMeth
Real*8 :: dqHdq, Disp, Fact
Real*8 :: StepMax=0.3D0
Real*8 :: Thr_RS=1.0D-7
Real*8 :: Beta_Disp=0.3D0
Integer, Parameter:: Max_Iter=50
Integer :: Iteration=0


#ifdef _DEBUGPRINT_
Write (6,*) 'Enter DIIS-GEK Optimizer'
#endif
If (.NOT.Init_LLs) Then
   Write (6,*) 'Link list not initiated'
   Call Abend()
End If

Call mma_allocate(q,mOV,iterso,Label='q')
Call mma_allocate(g,mOV,iterso,Label='g')

!Pick up coordinates and gradients in full space
nDIIS = 0
iFirst=iter-iterso+1
Do i = iFirst, iter
!  Write (6,*) 'i,iter=',i,iter
   nDIIS = nDIIS + 1

!  Coordinates
   ipq=LstPtr(i  ,LLx)
   q(:,nDIIS)=SCF_V(ipq)%A(:)

!  Gradients
   ipg=LstPtr(i  ,LLGrad)
   g(:,nDIIS)=SCF_V(ipg)%A(:)

End Do
#ifdef _DEBUGPRINT_
Call RecPrt('q',' ',q,mOV,iterso)
Call RecPrt('g',' ',g,mOV,iterso)
#endif


Call mma_allocate(e_diis,mOV,   iterso,Label='e_diis')
e_diis(:,:)=0.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find the coordinate, gradients and unit vectors of the reduced space
! Here we do this with the Gram-Schmidt procedure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set up the unit vectors from the gradient
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Do i = 1, nDIIS      ! normalize all the vectors
   gg = 0.0D0
   Do l = 1, mOV
      gg = gg + g(l,i)**2
   End Do
   e_diis(:,i) = g(:,i)/Sqrt(gg)
!   Write (6,*) i,i,DDot_(mOV,e_diis(:,i),1,e_diis(:,i),1)
End Do
!Write (6,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Gram-Schmidt ortho-normalization, eliminate linear dependence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
j = 1
Do i = 2, nDIIS
   Do k = 1, j
      gg = 0.0D0
      Do l = 1, mOV
         gg = gg + e_diis(l,i)*e_diis(l,k)
      End Do
      e_diis(:,i) = e_diis(:,i) - gg * e_diis(:,k)
   End Do
   gg = 0.0D0  ! renormalize  ! renormalize
   Do l = 1, mOV
      gg = gg + e_diis(l,i)**2
   End Do
!  Write (6,*) 'i,gg=',i,gg
   If (gg>1.0D-10) Then
      j = j + 1
      e_diis(:,j) = e_diis(:,i)/Sqrt(gg)
   End If
End Do
mDIIS=j
#ifdef _DEBUGPRINT_
Write (6,*) 'Check the ortonormality'
Do i = 1, mDIIS
   Do j = 1, i
      Write (6,*) i,j,DDot_(mOV,e_diis(:,i),1,e_diis(:,j),1)
   End Do
   Write (6,*)
End Do
Call RecPrt('e_diis',' ',e_diis,mOV,mDIIS)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computed the projected displacement coordinates. Note that the displacements are relative to the last coordinate, nDIIS.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Call mma_allocate(q_diis,mDIIS,nDIIS+Max_Iter,Label='q_diis')
q_diis(:,:)=0.0D0
Do i = 1, nDIIS
   Do k = 1, mDIIS
      gg = 0.0D0
      Do l = 1, mOV
         gg = gg + (q(l,i)-q(l,nDIIS))*e_diis(l,k)
      End Do
      q_diis(k,i)=gg
   End Do
End Do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computed the projected gradients
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Call mma_allocate(g_diis,mDIIS,nDIIS+Max_Iter,Label='g_diis')
g_diis(:,:)=0.0D0
Do i = 1, nDIIS
   Do k = 1, mDIIS
      gg = 0.0D0
      Do l = 1, mOV
         gg = gg + g(l,i)*e_diis(l,k)
      End Do
      g_diis(k,i)=gg
   End Do
End Do

#ifdef _DEBUGPRINT_
Call RecPrt('q_diis',' ',q_diis,mDIIS,nDIIS)
Call RecPrt('g_diis',' ',g_diis,mDIIS,nDIIS)
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Project the approximate Hessian to the subspace
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Call mma_allocate(H_diis,mDIIS,mDIIS,Label='H_diis')

Do i = 1, mDiis
   Do j = 1, mDiis
      gg = 0.0D0
      Do l = 1, mOV
         gg = gg + e_diis(l,i)*HDiag(l)*e_diis(l,j)
      End Do
      H_diis(i,j)=gg
   End Do
End Do
#ifdef _DEBUGPRINT_
Call RecPrt('H_diis((HDiag)',' ',H_diis,mDIIS,mDIIS)
#endif

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

Call mma_allocate(dq_diis,mDiis,Label='dq_Diis')

! We need to set the bias

Call Setup_Kriging(nDiis,mDiis,q_diis,g_diis,Energy(iFirst),H_diis)


!First implementation with a simple RS-RFO step
Iteration=nDiis
UpMeth=''
Step_Trunc=''
Fact=One
! Loop to enforce restricted variance
Do

   ! Compute the surrogate Hessian
   Call Hessian_Kriging_Layer(q(:,Iteration),H_diis,mDiis)
#ifdef _DEBUGPRINT_
   Call RecPrt('H_diis(updated)',' ',H_diis,mDIIS,mDIIS)
#endif

   Step_Trunc_=Step_Trunc
   dqHdq=Zero
   Call RS_RFO(H_diis,g_Diis(:,Iteration),mDiis,dq_diis,UpMeth,dqHdq,StepMax,Step_Trunc,Thr_RS)
   dq_diis(:)=-dq_diis(:)
   If (Step_Trunc//Step_Trunc_==' *') Step_Trunc='.'
#ifdef _DEBUGPRINT_
   Call RecPrt('dq_diis',' ',dq_diis,mDIIS,1)
#endif

   q_diis(:,Iteration+1) = q_diis(:,Iteration) + dq_diis(:)
   Call Dispersion_Kriging_Layer(q_diis(:,Iteration+1),Disp,mDIIS)
#ifdef _DEBUGPRINT_
   Write (6,*) 'Disp=',Disp
#endif

   Fact   =Half*Fact
   StepMax=Half*StepMax
   If (One-Disp/Beta_Disp>1.0D-3) Exit
   If ( (Fact<1.0D-5) .OR. (Disp<Beta_Disp) ) Exit
   Step_Trunc='*'

End Do

dq(:)=Zero
Do i = 1, nDIIS
   dq(:) = dq(:) + dq_diis(i)*e_diis(:,i)
End Do

#ifdef _DEBUGPRINT_
Call RecPrt('dq',' ',dq,mOV,1)
#endif

Call Finish_Kriging()
Call mma_deallocate(dq_diis)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

Call mma_deallocate(h_diis)

Call mma_deallocate(q_diis)
Call mma_deallocate(g_diis)
Call mma_deallocate(e_diis)

Call mma_deallocate(g)
Call mma_deallocate(q)

#ifdef _DEBUGPRINT_
Write (6,*) 'Exit DIIS-GEK Optimizer'
#endif
End Subroutine DIIS_GEK_Optimizer
