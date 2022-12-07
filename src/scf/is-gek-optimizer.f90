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
Subroutine IS_GEK_Optimizer(dq,mOV,dqdq,UpMeth,Step_Trunc,TestThr)
!***********************************************************************
!                                                                      *
!     Object: iterative-subspace gradient-enhanced kriging             *
!             optimization.                                            *
!                                                                      *
!     Author: Roland Lindh, Dept. of Chemistry -- BMC,                 *
!             University of Uppsala, SWEDEN                            *
!             May '22, November '22                                    *
!***********************************************************************
use InfSO , only: iterso, Energy
use InfSCF, only: iter
use LnkLst, only: SCF_V, Init_LLs, LLx, LLGrad
use SCF_Arrays, only: HDiag
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half, One, Four, Six
use Kriging_mod, only: blAI, mblAI, blaAI, blavAI
#ifdef _DEBUGPRINT_
use Constants, only: Two
#endif
Implicit None
Integer, Intent(In):: mOV
Real*8,  Intent(Out):: dq(mOV)
Real*8,  Intent(Out):: dqdq
Character(Len=1), Intent(InOut):: Step_Trunc
Character(Len=6), Intent(InOut):: UpMeth
Real*8, Intent(In) :: TestThr

!local variable
Integer i, j, k, l, ipq, ipg, nDIIS, mDIIS, iFirst
Integer, External:: LstPtr
Real*8, External::DDot_
Real*8, Allocatable:: q(:,:), g(:,:)
Real*8, Allocatable:: q_diis(:,:), g_diis(:,:), e_diis(:,:)
Real*8, Allocatable:: dq_diis(:), dq_0(:)
Real*8, Allocatable:: H_Diis(:,:), HDiag_Diis(:)
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
Real*8 :: ThrGrd=1.0D-7
Integer, Parameter:: Max_Iter=50
Integer :: Iteration=0
Integer :: Iteration_Micro=0
Integer :: Iteration_Total=0
Integer :: nWindow=5
Logical :: Converged=.FALSE., Terminate=.False.
Real*8, Allocatable :: Probe(:)
Real*8 :: Test
Integer, Save :: nExplicit=-1
#ifdef _DEBUGPRINT_
Real*8 :: eg
Real*8, Allocatable:: Vec(:,:)
Real*8, Allocatable:: Val(:)
#endif

Interface

   Subroutine SetUp_Kriging(nRaw,nInter,qInt,Grad,Energy,Hessian_HMF,HDiag)
   Integer nRaw, nInter
   Real*8 qInt(nInter,nRaw), Grad(nInter,nRaw), Energy(nRaw)
   Real*8, Optional:: Hessian_HMF(nInter,nInter)
   Real*8, Optional:: HDiag(nInter)
   End Subroutine SetUp_Kriging

End Interface


Beta_Disp=Beta_Disp_Seed
#ifdef _DEBUGPRINT_
Write (6,*) 'Enter IS-GEK Optimizer'
#endif
If (.NOT.Init_LLs) Then
   Write (6,*) 'Link list not initiated'
   Call Abend()
End If

Call mma_allocate(dq_0,mOV,Label='dq_0')
dq_0(:)=dq(:)
Call mma_allocate(q,mOV,iterso,Label='q')
Call mma_allocate(g,mOV,iterso,Label='g')

!Pick up coordinates and gradients in full space
iFirst=iter-Min(iterso,nWindow)+1
j = 0
Do i = iFirst, iter
!  Write (6,*) 'i,iter=',i,iter
j = i - iFirst + 1

!  Coordinates
   ipq=LstPtr(i  ,LLx)
   q(:,j)=SCF_V(ipq)%A(:)

!  Gradients
   ipg=LstPtr(i  ,LLGrad)
   g(:,j)=SCF_V(ipg)%A(:)

End Do
nDIIS = iter-iFirst+1
#ifdef _DEBUGPRINT_
Write (6,*) 'nWindow=',nWindow
Write (6,*) 'nDIIS=',nDIIS
Call RecPrt('q',' ',q,mOV,nDIIS)
Call RecPrt('g',' ',g,mOV,nDIIS)
#endif



! Generate the probe array which will guide us when we select which degrees of freedom to uncontract from the iterative subspace.
! Note that we have two different ways to do this.
!==================================================================================================================================
#define _Case1_
#ifdef _Case1_
Call mma_allocate(Probe,mOV,Label='Probe')
!==================================================================================================================================
Probe(:) = Zero
Do j = 1, nDIIS
   Do i = 1, mOV
      Probe(i)=Probe(i)+(g(i,j)/HDiag(i))**2
   End Do
End Do
!==================================================================================================================================
#else
! this need to be debugged
!==================================================================================================================================
Do j = 1, nDIIS
   Call SOrUpV(g(:,j),mOV,dq,'DISP','BFGS')
   Probe(:) = Probe(:) + dq(:)**2
End Do
!==================================================================================================================================
#endif
!==================================================================================================================================
! Finalize the average probe by dividing with nDIIS.
Do i = 1, mOV
   Probe(i)=Sqrt(Probe(i))/DBLE(nDIIS)
!  If (HDiag(i)<Zero) Write(6,*) HDiag(i)
End Do
#ifdef _DEBUGPRINT_
Call RecPrt('Probe',' ',Probe(:),mOV,1)
#endif


If (nExplicit<0) Then
nExplicit = 0
Do i = 1, mOV
   j = 0
   Test = TestThr
   Do k = 1, mOV
      If (Probe(k)>Test) Then
         Test=Probe(k)
         j=k
       End If
   End Do
   If (j==0) Then
      Exit
   Else
      nExplicit = nExplicit + 1
      Probe(j)=-Probe(j)
      If (nExplicit==Int(mOV*0.50D0)) Exit
   End If
End Do
Probe(:)=-Probe(:)
End If
#ifdef _DEBUGPRINT_
   Write (6,*) 'TestThr  :',TestThr
   Write (6,*) 'nExplicit:',nExplicit
   Write (6,*) 'mOV      :',mOV
#endif

Call mma_allocate(e_diis,mOV,   nDIIS+nExplicit,Label='e_diis')
e_diis(:,:)=Zero

!Call RecPrt('g',' ',g(:,nDIIS),1, mOV)
! Pick up the explicitly uncontracted components by selecting the largest elements in e_diis
mDIIS = 0
Do i = 1, nExplicit
   j = 0
   Test = Zero
   Do k = 1, mOV
      If (Probe(k)>Test) Then
         Test=Probe(k)
         j=k
       End If
   End Do
   If (j==0) Then
      Write(6,*) 'j==0'
      Call Abend()
   Else
      mDIIS = mDIIS + 1
      e_diis(j,mDIIS)=One
      Probe(j)=-Probe(j)
   End If
End Do

Call mma_deallocate(Probe)

Do i = 1, nDIIS      ! Pick up the gradient vectors as a seed for the reduced set of unit vectors.
   gg = 0.0D0
   Do l = 1, mOV
      gg = gg + g(l,i)**2
   End Do
   e_diis(:,i+nExplicit) = g(:,i)/Sqrt(gg)
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
Write (6,*) '      mOV:',mOV
Write (6,*) 'nExplicit:',nExplicit
Write (6,*) '    mDIIS:',mDIIS

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
Call mma_allocate(HDiag_diis,mDIIS,Label='HDiag_diis')

Do i = 1, mDiis
   Do j = 1, mDiis
      gg = 0.0D0
      Do l = 1, mOV
         gg = gg + e_diis(l,i)*HDiag(l)*e_diis(l,j)
      End Do
      H_diis(i,j)=gg
   End Do
   HDiag_Diis(i)=H_Diis(i,i)
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

!We need to set the bias

blavAI=10.00D0
Call Setup_Kriging(nDiis,mDiis,q_diis,g_diis,Energy(iFirst),HDiag=HDiag_Diis)
If (.False.) Write (6,*) blAI, mblAI, blaAI, blavAI
Call mma_deallocate(HDiag_diis)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Here starts the code doing the actual optimization

UpMeth='RVO'
Terminate=.False.
Step_Trunc='N'
Converged=.FALSE.

Iteration      =nDiis-1
Iteration_Micro=0
Iteration_Total=iter-1
If (nDIIS>1) Beta_Disp=Min(Beta_Disp_Seed,Max(Beta_Disp_Min,Abs(Energy(iter)-Energy(iter-1))))
#ifdef _DEBUGPRINT_
Write (6,*) '->',Energy(iter)-Energy(iter-1),nDIIS,Beta_Disp
#endif

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

#ifdef _DEBUGPRINT_
   Write (6,*)
   Write (6,*) '================================'
   Write (6,*) 'Micro Iteration=',Iteration_Micro
   Write (6,*) '================================'
   Write (6,*)
   Write (6,*) 'Step_Trunc:',Step_Trunc
   Write (6,*) '-----> Start RVO step'
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   Fact=One
   StepMax=StepMax_Seed
   ! Loop to enforce restricted variance. Note, if the step restriction kicks in no problem since we will still microiterate.
   ! Normally a full step will be allowed -- no step restriction -- and the loop will be exited after the first iteration.
   Do

      ! Compute the surrogate Hessian
     !Call Hessian_Kriging_Layer(q_diis(:,Iteration),H_diis,mDiis)
      Call Hessian_Kriging(q_diis(:,Iteration),H_diis,mDiis)

#ifdef _FOR_DEBUGGING_
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
#ifdef _DEBUGPRINT_
         Write (6,*) 'Eigenvalue:',Val(i*(i+1)/2)
#endif
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
#endif

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

     !Call Dispersion_Kriging_Layer(q_diis(:,Iteration+1),Variance,mDIIS)
      Call Dispersion_Kriging(q_diis(:,Iteration+1),Variance,mDIIS)

      Fact   =Half*Fact
      StepMax=Half*StepMax

      ! Note that we might have converged because the step restriction kicked in. However, we fill implicitly
      ! fix that during the second micro iteration.

#ifdef _DEBUGPRINT_
      Write (6,*)
      Write (6,*) 'Step_Trunc:',Step_Trunc
      Write (6,*) 'Beta_Disp =',Beta_Disp
      Write (6,*) 'Variance  =',Variance
      Write (6,*) 'Fact      =',Fact
      Write (6,*) 'StepMax   =',StepMax
#endif
      If (One-Variance/Beta_Disp>1.0D-3) Exit
      If ( (Fact<1.0D-5) .OR. (Variance<Beta_Disp) ) Exit
      Step_Trunc='*' ! This will only happen if the variance restriction kicks in

   End Do  ! Restricted variance step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef _DEBUGPRINT_
   Write (6,*)
   Write (6,*) 'Step_Trunc:',Step_Trunc
   Write (6,*) '-----> Exit RVO step'
   Write (6,*)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Compute the energy and the gradient of the surrogate model
  !Call Energy_Kriging_Layer(q_diis(:,Iteration+1),Energy(Iteration_Total+1),mDIIS)
   Call Energy_Kriging(q_diis(:,Iteration+1),Energy(Iteration_Total+1),mDIIS)
  !Call Gradient_Kriging_Layer(q_diis(:,Iteration+1),g_diis(:,Iteration+1),mDIIS)
   Call Gradient_Kriging(q_diis(:,Iteration+1),g_diis(:,Iteration+1),mDIIS)

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
#ifdef _DEBUGPRINT_
   Write (6,*) 'Step_Trunc:',Step_Trunc
   Write (6,*) 'Converged:', Converged
#endif
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

nExplicit=-1     ! Comment this statement to make the size static

#ifdef _DEBUGPRINT_
Call RecPrt('dq_diis',' ',dq_diis(:),SIZE(dq_diis),1)
Call RecPrt('dq',' ',dq(:),SIZE(dq),1)
Call RecPrt('g_diis(:,Iteration+1)',' ',g_diis(:,Iteration+1),SIZE(g_diis,1),1)
#endif

Call Finish_Kriging()
Call mma_deallocate(dq_diis)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!


! Merge the precomputed displacement in the full space with the step in the reduced space.
! First, subtract from the full space step any displacement components in the reduced space.

Write (6,*)  '  ||dq||=',Sqrt(DDot_(mOV,dq(:),1,dq(:),1))
Write (6,*)  '||dq_0||=',Sqrt(DDot_(mOV,dq_0(:),1,dq_0(:),1))
Write (6,*)
Write (6,*)  '<dq_0|e(dq)>/||dq_0||=',DDot_(mOV,dq_0(:),1,dq(:),1)/Sqrt(DDot_(mOV,dq(:),1,dq(:),1)) &
                                     / Sqrt(DDot_(mOV,dq_0(:),1,dq_0(:),1))
Do i = 1, mDIIS
   dq_0(:) = dq_0(:) - DDot_(mOV,dq_0(:),1,e_diis(:,i),1) * e_diis(:,i)
End Do

!dq_0(:) = dq_0(:) - ( DDot_(mOV,dq_0(:),1,dq(:),1) / DDot_(mOV,dq(:),1,dq(:),1) ) * dq(:)

! Second, add the projected complementary full space step to the step in the reduced space.

Write (6,*)  '||dq_0||=',Sqrt(DDot_(mOV,dq_0(:),1,dq_0(:),1))
dq(:) = dq(:) + dq_0(:)
dqdq=Sqrt(DDot_(SIZE(dq),dq(:),1,dq(:),1))
Write (6,*)  '  ||dq||=',Sqrt(DDot_(mOV,dq(:),1,dq(:),1))


Call mma_deallocate(h_diis)

Call mma_deallocate(q_diis)
Call mma_deallocate(g_diis)
Call mma_deallocate(e_diis)

Call mma_deallocate(g)
Call mma_deallocate(q)
Call mma_deallocate(dq_0)

#ifdef _DEBUGPRINT_
Write (6,*) 'Exit IS-GEK Optimizer'
#endif
End Subroutine IS_GEK_Optimizer
