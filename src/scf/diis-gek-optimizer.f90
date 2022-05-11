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
Subroutine DIIS_GEK_Optimizer(Disp,mOV)
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
Real*8,  Intent(Out):: Disp(mOV)

Integer i, j, k, l, ipq, ipg, nDIIS, iFirst
Integer, External:: LstPtr
Real*8, External::DDot_
Real*8, Allocatable:: q(:,:), g(:,:)
Real*8, Allocatable:: q_diis(:,:), g_diis(:,:), e_diis(:,:)
Real*8, Allocatable:: dq_diis(:)
Real*8, Allocatable:: H_Diis(:,:)
Real*8 :: gg
Character(Len=1) Step_Trunc
Character(Len=6) UpMeth
Real*8 :: dqHdq, StepMax, Thr_RS


!Write (6,*) 'Enter DIIS-GEK Optimizer'
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
!Call RecPrt('q',' ',q,mOV,iterso)
!Call RecPrt('g',' ',g,mOV,iterso)


Call mma_allocate(q_diis,iterso,iterso,Label='q_diis')
q_diis(:,:)=0.0D0
Call mma_allocate(g_diis,iterso,iterso,Label='g_diis')
g_diis(:,:)=0.0D0
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
! Gram-Schmidt ortho-normalization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Do i = 2, nDIIS
   Do k = 1, i-1
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
   e_diis(:,i) = e_diis(:,i)/Sqrt(gg)
End Do
!Call RecPrt('e_diis',' ',e_diis,mOV,iterso)
!Do i = 1, nDIIS
!   Do j = 1, i
!      Write (6,*) i,j,DDot_(mOV,e_diis(:,i),1,e_diis(:,j),1)
!   End Do
!End Do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computed the projected coordinates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Do i = 1, nDIIS
   Do k = 1, nDIIS
      gg = 0.0D0
      Do l = 1, mOV
         gg = gg + (q(l,i)-q(l,1))*e_diis(l,k)
      End Do
      q_diis(k,i)=gg
   End Do
End Do
Call RecPrt('q_diis',' ',q_diis,iterso,iterso)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computed the projected gradients
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Do i = 1, nDIIS
   Do k = 1, i
      gg = 0.0D0
      Do l = 1, mOV
         gg = gg + g(l,i)*e_diis(l,k)
      End Do
      g_diis(k,i)=gg
   End Do
End Do
Call RecPrt('g_diis',' ',g_diis,iterso,iterso)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Project the approximate Hessian to the subspace
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Call mma_allocate(H_diis,nDIIS,nDIIS,Label='H_diis')

Do i = 1, nDiis
   Do j = 1, nDiis
      gg = 0.0D0
      Do l = 1, mOV
         gg = gg + e_diis(l,i)*HDiag(l)*e_diis(l,j)
      End Do
      H_diis(i,j)=gg
   End Do
End Do
Call RecPrt('H_diis',' ',H_diis,iterso,iterso)

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

Call mma_allocate(dq_diis,nDiis,Label='dq_Diis')

! We need to set the bias

Call Setup_Kriging(nDiis,nDiis,q_diis,g_diis,Energy(iFirst),H_diis)

! Compute the surrogate Hessian

Call Hessian_Kriging_Layer(q(:,nDiis),H_diis,nDiis)
Call RecPrt('H_diis',' ',H_diis,iterso,iterso)

!First implementation with a simple RS-RFO step

!#define _NEWCODE_
#ifdef _NEWCODE_
dqHdq=Zero
StepMax=0.3D0
Thr_RS=1.0D-7
UpMeth=''
Step_Trunc=''
Call RS_RFO(H_diis,g_Diis(:,nDiis),nDiis,dq_diis,UpMeth,dqHdq,StepMax,Step_Trunc,Thr_RS)
dq_diis(:)=-dq_diis(:)
Call RecPrt('dq_diis',' ',dq_diis,nDIIS,1)

Disp(:)=Zero
Do i = 1, nDIIS
   Disp(:) = Disp(:) + dq_diis(i)*e_diis(:,i)
End Do
#else
! This code gives right values on the first iterations
Block
Real*8, Allocatable:: Hessian(:,:)
Call mma_allocate(Hessian,mOV,mOV)
Hessian(:,:)=Zero
Do i = 1, mOV
  Hessian(i,i)=HDiag(i)
End Do
Call RS_RFO(Hessian,g(:,nDiis),mOV,Disp,UpMeth,dqHdq,StepMax,Step_Trunc,Thr_RS)
Disp(:)=-Disp(:)
Call mma_deallocate(Hessian)
End Block
#endif

Call RecPrt('Disp',' ',Disp,mOV,1)

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

!Write (6,*) 'Exit DIIS-GEK Optimizer'
End Subroutine DIIS_GEK_Optimizer
