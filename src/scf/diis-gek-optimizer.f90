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
Subroutine DIIS_GEK_Optimizer(dq,mOV,dqdq,UpMeth,Step_Trunc)
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
!use kriging_mod, only: blAI, mblAI, blaAI, blavAI
Implicit None
#include "real.fh"
#include "stdalloc.fh"
Integer, Intent(In):: mOV
Real*8,  Intent(Out):: dq(mOV)
Real*8,  Intent(Out):: dqdq
Character(Len=1), Intent(InOut):: Step_Trunc
Character(Len=6), Intent(InOut):: UpMeth

Integer i, j, k, l, ipq, ipg, nDIIS, mDIIS, iFirst
Integer, External:: LstPtr
Real*8, External::DDot_
Real*8, Allocatable:: q(:,:), g(:,:)
Real*8, Allocatable:: q_diis(:,:), g_diis(:,:), e_diis(:,:)
Real*8, Allocatable:: dq_diis(:)
Real*8, Allocatable:: H_Diis(:,:)
Real*8 :: gg
Character(Len=1) Step_Trunc_
Character(Len=6) UpMeth_
Real*8 :: dqHdq, Variance, Fact
Real*8 :: StepMax=0.D0
Real*8 :: StepMax_Seed=0.3D0
Real*8 :: Thr_RS=1.0D-7
Real*8 :: Beta_Disp=0.3D0
Real*8 :: FAbs, RMS, RMSMx, dEner
Real*8 :: ThrGrd=1.0D-6
Integer, Parameter:: Max_Iter=50
Integer :: Iteration=0
Integer :: Iteration_Micro=0
Integer :: Iteration_Total=0
Logical :: Converged=.FALSE.


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

UpMeth='RVO'

Step_Trunc=' '
Converged=.FALSE.

Iteration      =nDiis-1
Iteration_Micro=0
Iteration_Total=iter-1
Do While (.NOT.Converged) ! Micro iterate on the surrogate model
   Iteration_Micro = Iteration_Micro + 1
   Iteration_Total = Iteration_Total + 1
   Iteration       = Iteration       + 1
   If (Iteration_Micro==Max_Iter) Then
      Write (6,*)
      Write (6,*)  'DIIS-GEK-Optimizer: Iteration_Micro==Max_Iter'
      Write (6,*)  'Abend!'
      Write (6,*)
      Call Abend()
   End If

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _DEBUGPRINT_
   Write (6,*)
   Write (6,*) '================================'
   Write (6,*) 'Micro Iteration=',Iteration_Micro
   Write (6,*) '================================'
   Write (6,*)
   Write (6,*) '-----> Start RVO step'
#endif

   Fact=One
   StepMax=StepMax_Seed
   ! Loop to enforce restricted variance. Note, if the step restriction kicks no problem since we will still microiterate.
   Do

   ! Compute the surrogate Hessian
      Call Hessian_Kriging_Layer(q_diis(:,Iteration),H_diis,mDiis)
#ifdef _DEBUGPRINT_
      Call RecPrt('q_diis(:,Iteration)',' ',q_diis(:,Iteration),mDIIS,1)
      Call RecPrt('H_diis(updated)',' ',H_diis,mDIIS,mDIIS)
#endif

      Step_Trunc_=Step_Trunc
      Step_Trunc ='N'   ! set to not defined
      dqHdq=Zero
      Call RS_RFO(H_diis,g_Diis(:,Iteration),mDiis,dq_diis,UpMeth_,dqHdq,StepMax,Step_Trunc,Thr_RS)
      dq_diis(:)=-dq_diis(:)
      q_diis(:,Iteration+1) = q_diis(:,Iteration) + dq_diis(:)

#ifdef _DEBUGPRINT_
      Write (6,*)
      Write (6,*) 'Subiteration: Step_Trunc, StepMax:',Step_Trunc,StepMax
      Call RecPrt('dq_diis',' ',dq_diis,mDIIS,1)
      Call RecPrt('q_diis(:,Iteration+1)',' ',q_diis(:,Iteration+1),mDIIS,1)
#endif
      If (Step_Trunc.eq.'N') Step_Trunc=' '
      If (Step_Trunc//Step_Trunc_==' *') Step_Trunc='.'

      Call Dispersion_Kriging_Layer(q_diis(:,Iteration+1),Variance,mDIIS)
#ifdef _DEBUGPRINT_
      Write (6,*)
      Write (6,*) 'Beta_Disp=',Beta_Disp
      Write (6,*) 'Variance=',Variance
#endif

      Fact   =Half*Fact
      StepMax=Half*StepMax
      If (One-Variance/Beta_Disp>1.0D-3) Exit
      If ( (Fact<1.0D-5) .OR. (Variance<Beta_Disp) ) Exit
      Step_Trunc='*'

   End Do  ! Restricted variance
#ifdef _DEBUGPRINT_
   Write (6,*) '-----> Exit RVO step'
   Write (6,*)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Compute the energy and the gradient of the surrogate model
   Call Energy_Kriging_Layer(q_diis(:,Iteration+1),Energy(Iteration_Total+1),mDIIS)
   dEner = Energy(Iteration_Total+1) - Energy(Iteration_Total)
   Call Gradient_Kriging_Layer(q_diis(:,Iteration+1),g_diis(:,Iteration+1),mDIIS)

#ifdef _DEBUGPRINT_
   Write (6,*) 'Energy(Iteration_Total+1):',Energy(Iteration_Total+1)
   Call RecPrt('g_diis(:,Iteration+1)',' ',g_diis(:,Iteration+1),mDIIS,1)
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Check on convergence criteria.
!
   FAbs=Sqrt(DDot_(mDIIS,g_diis(:,Iteration+1),1,g_diis(:,Iteration+1),1)/DBLE(mDIIS))
   RMS =SQRT(DDot_(mDIIS,dq_diis(:),1,dq_diis(:),1)/DBLE(mDIIS))
   RMSMx=Zero
   Do i = 1, mDIIS
      RMSMx=Max(RMSMx,Abs(dq_diis(i)))
   End Do
#ifdef _DEBUGPRINT_
   Write (6,*)
   Write (6,*) 'FAbs=',FAbs
   Write (6,*) 'RMS=',RMS
   Write (6,*) 'RMSMx=',RMSMx
   Write (6,*)
#endif
   Converged = FAbs<ThrGrd
   Converged = Converged .AND. RMS<Four*ThrGrd
   Converged = Converged .AND. RMSMx<ThrGrd*Six
   Converged = Converged .AND. Step_Trunc==' '
   If (Step_Trunc.eq.'.') Step_Trunc=' '
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

End Do  ! While not converged
Write (UpMeth(5:6),'(I2)') Iteration_Micro


! Compute the displacement in the reduced space relative to the last structure of the full space
dq_diis(:)=q_diis(:,Iteration+1)-q_diis(:,nDIIS)
! Compute the displacement in the full space.
dq(:)=Zero
Do i = 1, SIZE(e_diis,2)
   dq(:) = dq(:) + dq_diis(i)*e_diis(:,i)
End Do
dqdq=Sqrt(DDot_(SIZE(dq),dq(:),1,dq(:),1))

#ifdef _DEBUGPRINT_
Call RecPrt('dq',' ',dq,SIZE(dq),1)
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
