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
Real*8, Allocatable:: Vec(:,:)
Real*8, Allocatable:: Val(:)
Real*8 :: gg
Character(Len=1) Step_Trunc_
Character(Len=6) UpMeth_
Real*8 :: dqHdq, Variance, Fact
Real*8 :: StepMax=0.D0
Real*8 :: StepMax_Seed=0.1D0
Real*8 :: Thr_RS=1.0D-7
Real*8 :: Beta_Disp_Seed=0.05D0
Real*8 :: Beta_Disp_Min=1.0D-3
Real*8 :: Beta_Disp
Real*8 :: FAbs, RMS, RMSMx
!Real*8 :: dEner
Real*8 :: ThrGrd=1.0D-6
Integer, Parameter:: Max_Iter=50
Integer :: Iteration=0
Integer :: Iteration_Micro=0
Integer :: Iteration_Total=0
Integer :: nWindow=10
Logical :: Converged=.FALSE., Terminate=.False.
Integer :: nExplicit
Real*8, Allocatable :: Probe(:)
Real*8 :: Test

Beta_Disp=Beta_Disp_Seed
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
iFirst=iter-Min(iterso,nWindow)+1
iFirst=iter-iterso+1    ! Disable the window
nDIIS = 0
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
Call RecPrt('q',' ',q,mOV,nDIIS)
Call RecPrt('g',' ',g,mOV,nDIIS)
#endif



!#define _REDUCED_SPACE_
#ifdef _REDUCED_SPACE_

Call mma_allocate(e_diis,mOV,   nDIIS,Label='e_diis')
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

#else

!#define _FULL_SPACE_
#ifdef _FULL_SPACE_

Call mma_allocate(e_diis,mOV,   mOV,Label='e_diis')
e_diis(:,:)=0.0D0
mDIIS=mOV
Do i = 1, mDIIS
   e_diis(i,i)=1.0D0
End Do

#else

Call mma_allocate(Probe,mOV,Label='Probe')

Probe(:) = Zero
#define _Case1_
#ifdef _Case1_
Do j = 1, nDIIS
!Do j = nDIIS, nDIIS
!Do j = Max(1,nDIIS-3), nDIIS
!  If (nDIIS==1) Then
      Do i = 1, mOV
!        Probe(i)=Probe(i)+g(i,j)**2
         Probe(i)=Probe(i)+(g(i,j)/HDiag(i))**2
      End Do
!  Else
!     Do i = 1, mOV
!        Probe(i)=Probe(i)+(g(i,j)-g(i,j-1))**2
!     End Do
!  End If
End Do
#else
Do j = 1, nDIIS
   Call SOrUpV(g(:,j),mOV,dq,'DISP','BFGS')
   Probe(:) = Probe(:) + dq(:)**2
End Do
#endif
!Write (6,*) 'Check HDiag'
Do i = 1, mOV
   Probe(i)=Sqrt(Probe(i))/DBLE(nDIIS)
!  If (HDiag(i)<Zero) Write(6,*) HDiag(i)
End Do
!Call RecPrt('Probe',' ',Probe(:),mOV,1)

!nExplicit=mOV    ! Full GEK
nExplicit=Min(mOV,10) ! Partial GEK
nExplicit=Min(mOV,20) ! Partial GEK
!nExplicit=0 ! DIIS-GEK
Call mma_allocate(e_diis,mOV,   nDIIS+nExplicit,Label='e_diis')
e_diis(:,:)=Zero

!Call RecPrt('g',' ',g(:,nDIIS),1, mOV)
! Pick up the explicitly uncontracted components
mDIIS = 0
Do i = 1, nExplicit
   j = 0
   Test = Zero
   Do k = 1, mOV
      If (Probe(k)>Test) Then
         Test=Probe(k)
         j=k
       Else If (Probe(k)<Zero) Then
         Write (6,*) 'Probe(k)<Zero'
         Call Abend()
       End If
   End Do
   If (j==0) Then
      Write(6,*) 'j==0'
      Call Abend()
   Else
!     Write (6,*) 'j,Probe(j)=',j, Probe(j)
      mDIIS = mDIIS + 1
      e_diis(j,mDIIS)=One
      Probe(j)=Zero
   End If
End Do

Call mma_deallocate(Probe)

Do i = 1, nDIIS      ! Pick up the gradient vectors as a seed for the reduced set of unit vectors.
   gg = 0.0D0
   Do l = 1, mOV
      gg = gg + g(l,i)**2
   End Do
!  If (gg<1.0D-10) Then
!     e_diis(:,i+nExplicit) = Zero
!  Else
      e_diis(:,i+nExplicit) = g(:,i)/Sqrt(gg)
!  End If
!  Write (6,*) i,i,DDot_(mOV,e_diis(:,i),1,e_diis(:,i),1)
End Do

! now orthogonalize all vectors
j = 1
Do i = 2, nDIIS + nExplicit
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
   If (gg>1.0D-10) Then   ! Skip vector if linear dependent.
      j = j + 1
      e_diis(:,j) = e_diis(:,i)/Sqrt(gg)
   End If
End Do
mDIIS=j

#endif

#endif

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
! Computed the projected displacement coordinates. Note that the displacements are relative to the last coordinate, q(:,nDIIS).
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
!Write (6,*) blAI, mblAI, blaAI, blavAI

UpMeth='RVO'
Terminate=.False.
Step_Trunc='N'
Converged=.FALSE.

Iteration      =nDiis-1
Iteration_Micro=0
Iteration_Total=iter-1
If (nDIIS>1) Beta_Disp=Min(Beta_Disp_Seed,Max(Beta_Disp_Min,Abs(Energy(iter)-Energy(iter-1))))
!Write (6,*) '->',Energy(iter)-Energy(iter-1),nDIIS,Beta_Disp
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
   Write (6,*) 'Step_Trunc:',Step_Trunc
   Write (6,*) '-----> Start RVO step'
#endif

   Fact=One
   StepMax=StepMax_Seed
   ! Loop to enforce restricted variance. Note, if the step restriction kicks no problem since we will still microiterate.
   Do

      ! Compute the surrogate Hessian
      Call Hessian_Kriging_Layer(q_diis(:,Iteration),H_diis,mDiis)

      Call mma_allocate(Val,mDIIS*(mDIIS+1)/2,Label='Val')
      Call mma_allocate(Vec,mDIIS,mDIIS,Label='Vec')

      Vec(:,:)=Zero
      Do i = 1, mDIIS
         Vec(:,:)=One
         Do j = 1, i
            Val((i-1)*i/2+j)=H_diis(i,j)
         End Do
      End Do

      Call NIDiag_new(Val,Vec,mDIIS,mDIIS)
      Call Jacord(Val,Vec,mDIIS,mDIIS)

!     If negative eigenvalues then correct and signal that the micro iterartions should be terminanted.
      Do i = 1, mDIIS
!        Write (*,*) Val(i*(i+1)/2)
         If (Val(i*(i+1)/2)<Zero) Then
            Terminate=.True.
            Do j = 1, mDIIS
               Do k = 1, mDIIS
                  H_Diis(j,k) = H_Diis(j,k) + Two*Abs(Val(i*(i+1)/2))*Vec(j,i)*Vec(k,i)
               End Do
            End Do
         End If
      End Do

      Call mma_deallocate(Vec)
      Call mma_deallocate(Val)

#ifdef _DEBUGPRINT_
      Call RecPrt('q_diis(:,Iteration)',' ',q_diis(:,Iteration),mDIIS,1)
      Call RecPrt('H_diis(updated)',' ',H_diis,mDIIS,mDIIS)
      Write (6,*) 'Step_Trunc:',Step_Trunc
#endif

      Step_Trunc_=Step_Trunc
      dqHdq=Zero
      Call RS_RFO(H_diis,g_Diis(:,Iteration),mDiis,dq_diis,UpMeth_,dqHdq,StepMax,Step_Trunc_,Thr_RS)
      dq_diis(:)=-dq_diis(:)
      q_diis(:,Iteration+1) = q_diis(:,Iteration) + dq_diis(:)
      dqdq=Sqrt(DDot_(SIZE(dq_diis),dq_diis(:),1,dq_diis(:),1))

#ifdef _DEBUGPRINT_
      Write (6,*)
      Write (6,*) 'Subiteration: Step_Trunc, StepMax:',Step_Trunc,StepMax, dqdq
      Write (6,*) 'Subiteration: Step_Trunc_        :',Step_Trunc_
      Call RecPrt('dq_diis',' ',dq_diis,mDIIS,1)
      Call RecPrt('q_diis(:,Iteration+1)',' ',q_diis(:,Iteration+1),mDIIS,1)
#endif
      If (Step_Trunc.eq.'N') Step_Trunc=' '   ! set to blank if not touched
      If (Step_Trunc//Step_Trunc_==' *') Step_Trunc='.' ! Mark that we have had a step Reduction

      Call Dispersion_Kriging_Layer(q_diis(:,Iteration+1),Variance,mDIIS)
#ifdef _DEBUGPRINT_
      Write (6,*)
      Write (6,*) 'Step_Trunc:',Step_Trunc
      Write (6,*) 'Beta_Disp=',Beta_Disp
      Write (6,*) 'Variance=',Variance
#endif

      Fact   =Half*Fact
      StepMax=Half*StepMax

      ! Note that we might have converged because the step restriction kicked in. However, we fill implicitly
      ! fix that during the second micro iteration.

      If (One-Variance/Beta_Disp>1.0D-3) Exit
      If ( (Fact<1.0D-5) .OR. (Variance<Beta_Disp) ) Exit
      Step_Trunc='*' ! This will only happen if variance restriction kicks in

   End Do  ! Restricted variance
#ifdef _DEBUGPRINT_
   Write (6,*)
   Write (6,*) 'Step_Trunc:',Step_Trunc
   Write (6,*) '-----> Exit RVO step'
   Write (6,*)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Compute the energy and the gradient of the surrogate model
   Call Energy_Kriging_Layer(q_diis(:,Iteration+1),Energy(Iteration_Total+1),mDIIS)
   Call Gradient_Kriging_Layer(q_diis(:,Iteration+1),g_diis(:,Iteration+1),mDIIS)

!  dEner = Energy(Iteration_Total+1) - Energy(Iteration_Total)

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
   If (Step_Trunc=='.') Step_Trunc=' '
   Converged = FAbs<ThrGrd
   Converged = Converged .AND. RMS<Four*ThrGrd
   Converged = Converged .AND. RMSMx<ThrGrd*Six
   Converged = Converged .AND. (Step_Trunc==' ' .or. Step_Trunc=='#')
!  Write (6,*) 'Step_Trunc:',Step_Trunc
   If (Step_Trunc=='*') Converged=.True.
   If (Terminate) Then
      Step_Trunc='#'
!     Exit
   End If
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

End Do  ! While not converged
#ifdef _DEBUGPRINT_
Write (6,*) 'Converged'
Write (6,*) 'Energy(Iteration_Total+1):',Energy(Iteration_Total+1)
#endif
Write (UpMeth(5:6),'(I2)') Iteration_Micro


! Compute the displacement in the reduced space relative to the last structure of the full space
dq_diis(:)=q_diis(:,Iteration+1)-q_diis(:,nDIIS)
! Compute the displacement in the full space.
dq(:)=Zero
Do i = 1, mDIIS
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
