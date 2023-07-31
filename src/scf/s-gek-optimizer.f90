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
!#define _FULL_SPACE_ ! Debugging
!#define _KRYLOV_
!#define _HYBRID_
!#define _HYBRID2_
#define _HYBRID3_
Subroutine S_GEK_Optimizer(dq,mOV,dqdq,UpMeth,Step_Trunc)
!***********************************************************************
!                                                                      *
!     Object: subspace gradient-enhanced kriging optimization.         *
!                                                                      *
!     Author: Roland Lindh, Dept. of Chemistry -- BMC,                 *
!             University of Uppsala, SWEDEN                            *
!             May '22, November-December '22                           *
!***********************************************************************
use InfSO , only: iterso, Energy
use InfSCF, only: iter
use LnkLst, only: SCF_V, Init_LLs, LLx, LLGrad
use SCF_Arrays, only: HDiag
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half, One, Four, Six
use Kriging_mod, only: blAI, mblAI, blaAI, blavAI
use Kriging_procedures, only: Setup_Kriging
use Constants, only: Two
Implicit None
Integer, Intent(In):: mOV
Real*8,  Intent(InOut):: dq(mOV)
Real*8,  Intent(Out):: dqdq
Character(Len=1), Intent(InOut):: Step_Trunc
Character(Len=6), Intent(InOut):: UpMeth

!local variable
Integer i, j, k, l, ipq, ipg, nDIIS, mDIIS, iFirst
Integer, External:: LstPtr
Integer :: IterSO_Save, Iter_Save
Real*8, External::DDot_
Real*8, Allocatable:: q(:,:), g(:,:)
Real*8, Allocatable:: q_diis(:,:), g_diis(:,:), e_diis(:,:)
Real*8, Allocatable:: dq_diis(:), aux_a(:), aux_b(:)
Real*8, Allocatable:: H_Diis(:,:), HDiag_Diis(:)
Real*8 :: gg
Character(Len=1) Step_Trunc_
Character(Len=6) UpMeth_
Real*8 :: dqHdq, Variance(1), Fact
Real*8 :: StepMax=0.D0
Real*8 :: StepMax_Seed=0.1D0
Real*8 :: Thr_RS=1.0D-7
Real*8 :: Beta_Disp_Seed=0.05D0
Real*8 :: Beta_Disp_Min=5.0D-3
Real*8 :: Beta_Disp
Real*8 :: FAbs, RMS, RMSMx
!Real*8 :: dEner
Real*8 :: ThrGrd=1.0D-7
Integer, Parameter:: Max_Iter=50
Integer :: Iteration=0
Integer :: Iteration_Micro=0
Integer :: Iteration_Total=0
Integer :: nWindow=8
#if defined(_KRYLOV_)
Integer :: nKrylov=20
#endif
Logical :: Converged=.FALSE., Terminate=.False.
Integer :: nExplicit
Real*8, Allocatable:: Vec(:,:)
Real*8, Allocatable:: Val(:)

Beta_Disp=Beta_Disp_Seed
#ifdef _DEBUGPRINT_
Write (6,*) 'Enter S-GEK Optimizer'
#endif
If (.NOT.Init_LLs) Then
   Write (6,*) 'Link list not initiated'
   Call Abend()
End If

Call mma_allocate(q,mOV,iterso,Label='q')
Call mma_allocate(g,mOV,iterso,Label='g')

!Pick up coordinates and gradients in full space
iFirst=iter-Min(iterso,nWindow)+1
j = 0
Do i = iFirst, iter
j = i - iFirst + 1
 !Write (6,*) 'i,j,iter=',i,j,iter

!  Coordinates
   ipq=LstPtr(i  ,LLx)
   q(:,j)=SCF_V(ipq)%A(:)

!  Gradients
   ipg=LstPtr(i  ,LLGrad)
   g(:,j)=SCF_V(ipg)%A(:)

End Do

nDIIS = iter-iFirst+1

If (nDIIS==1) Then
#ifdef _DEBUGPRINT_
Write (6,*) 'Exit S-GEK Optimizer'
#endif
Call mma_deallocate(g)
Call mma_deallocate(q)
Return
End If

#ifdef _DEBUGPRINT_
Write (6,*) 'nWindow=',nWindow
Write (6,*) 'nDIIS=',nDIIS
Write (6,*) 'IterSO=',IterSO
Call RecPrt('q',' ',q,mOV,nDIIS)
Call RecPrt('g',' ',g,mOV,nDIIS)
Call RecPrt('g(:,nDIIS)',' ',g(:,nDIIS),mOV,1)
#endif

#if defined(_FULL_SPACE_)

! Set up the full space
nExplicit = mOV
Call mma_allocate(e_diis,mOV, nExplicit,Label='e_diis')
Do k = 1, nExplicit
   e_diis(:,k)=Zero
   e_diis(k,k)=One
End Do

#elif defined(_KRYLOV_)

! Set up unit vectors corresponding to a Krylov subspace for Adx=g
nExplicit = Min(nKrylov,mOV)
!nExplicit = mOV
Call mma_allocate(e_diis,mOV, nExplicit,Label='e_diis')
Call mma_allocate(Aux_a,mOV,Label='Aux_a')
Call mma_allocate(Aux_b,mOV,Label='Aux_b')
j=1
Aux_a(:) = dq(:)
e_diis(:,j) = Aux_a(:) / Sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))
j = j + 1
Aux_a(:) = g(:,nDIIS)
e_diis(:,j) = Aux_a(:) / Sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))
do j = 3, nExplicit
   Call SOrUpV(Aux_a(:),mOV,Aux_b(:),'GRAD','BFGS')
   e_diis(:,j) = Aux_b(:) / Sqrt(DDot_(mOV,Aux_b(:),1,Aux_b(:),1))
   Aux_a(:)=Aux_b(:)
end do
Call mma_deallocate(Aux_b)
Call mma_deallocate(Aux_a)

#elif defined(_HYBRID_)

! Set up unit vectors corresponding to the subspace which the BFGS update will span.
nExplicit = Min(2 * (nDIIS-1) + nKrylov,mOV)
Call mma_allocate(e_diis,mOV, nExplicit,Label='e_diis')

Call mma_allocate(Aux_a,mOV,Label='Aux_a')
Call mma_allocate(Aux_b,mOV,Label='Aux_b')
IterSO_save = IterSO
Iter_save = Iter
Iter = iFirst
IterSO = 1

j=0
Do k = 1, nDIIS-1
   j = j + 1
   Aux_a(:) = g(:,k+1) - g(:,k)
   e_diis(:,j) = Aux_a(:) / Sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))

   j = j + 1
   Aux_a(:) = q(:,k+1) - q(:,k)
   Call SOrUpV(Aux_a(:),mOV,Aux_b(:),'GRAD','BFGS')
   e_diis(:,j) = Aux_b(:) / Sqrt(DDot_(mOV,Aux_b(:),1,Aux_b(:),1))

   iter=iter+1
   iterSO=iterSO+1
End Do
IterSO = IterSO_save
Iter = Iter_save

!Add some unit vectors correponding to the Krylov subspace algorithm, g, Ag, A^2g, ....
j = j + 1
Aux_a(:) = dq(:)
e_diis(:,j) = Aux_a(:) / Sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))
j = j + 1
Aux_a(:) = g(:,nDIIS)
e_diis(:,j) = Aux_a(:) / Sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))
do k = j+1, nExplicit
   Call SOrUpV(Aux_a(:),mOV,Aux_b(:),'GRAD','BFGS')
   e_diis(:,k) = Aux_b(:) / Sqrt(DDot_(mOV,Aux_b(:),1,Aux_b(:),1))
   Aux_a(:)=Aux_b(:)
end do
Call mma_deallocate(Aux_b)
Call mma_deallocate(Aux_a)

#elif defined(_HYBRID2_)

! Set up unit vectors corresponding to the subspace which the BFGS update will span.
nExplicit = Min(2 * nDIIS - 1 + nKrylov,mOV)
Call mma_allocate(e_diis,mOV, nExplicit,Label='e_diis')

Call mma_allocate(Aux_a,mOV,Label='Aux_a')
Call mma_allocate(Aux_b,mOV,Label='Aux_b')
IterSO_save = IterSO
Iter_save = Iter
Iter = iFirst
IterSO = 1

j=0
Do k = 1, nDIIS
   j = j + 1
   Aux_a(:) = g(:,k)
   e_diis(:,j) = Aux_a(:) / Sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))

   if (j==2*nDIIS-1) Exit
   j = j + 1
   Aux_a(:) = q(:,k)
   Call SOrUpV(Aux_a(:),mOV,Aux_b(:),'GRAD','BFGS')
   Call RecPrt('Aux_b(:)',' ',Aux_b(:),1,mOV)
   e_diis(:,j) = Aux_b(:) / Sqrt(DDot_(mOV,Aux_b(:),1,Aux_b(:),1))

   Iter=Min(Iter+1,Iter_Save)
   IterSO=Min(IterSO+1,IterSO_Save)
End Do
IterSO = IterSO_save
Iter = Iter_save

!Add some unit vectors correponding to the Krylov subspace algorithm, g, Ag, A^2g, ....
j = j + 1
Aux_a(:) = dq(:)
e_diis(:,j) = Aux_a(:) / Sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))
j = j + 1
Aux_a(:) = g(:,nDIIS)
e_diis(:,j) = Aux_a(:) / Sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))
Aux_a(:)=e_diis(:,j)
do k = j+1, nExplicit
   Call SOrUpV(Aux_a(:),mOV,Aux_b(:),'GRAD','BFGS')
   e_diis(:,k) = Aux_b(:) / Sqrt(DDot_(mOV,Aux_b(:),1,Aux_b(:),1))
   Aux_a(:)=Aux_b(:)
   Aux_a(:)=e_diis(:,k)
end do
Call mma_deallocate(Aux_b)
Call mma_deallocate(Aux_a)

#elif defined(_HYBRID3_)

!nExplicit = 2 * (nDIIS - 1) + mOV + 2
nExplicit = 2 * (nDIIS - 1) + 2
Call mma_allocate(e_diis,mOV, nExplicit,Label='e_diis')

Call mma_allocate(Aux_a,mOV,Label='Aux_a')
Call mma_allocate(Aux_b,mOV,Label='Aux_b')
IterSO_save = IterSO
Iter_save = Iter
Iter = iFirst
IterSO = 1

j=0
Do k = 1, nDIIS-1
   j = j + 1
   Aux_a(:) = g(:,k+1) - g(:,k)
   e_diis(:,j) = Aux_a(:) / Sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))

   j = j + 1
   Aux_a(:) = q(:,k+1) - q(:,k)
!  Call SOrUpV(Aux_a(:),mOV,Aux_b(:),'GRAD','BFGS')
   Aux_b(:)=Aux_a(:)
   e_diis(:,j) = Aux_b(:) / Sqrt(DDot_(mOV,Aux_b(:),1,Aux_b(:),1))

   iter=iter+1
   iterSO=iterSO+1
End Do
IterSO = IterSO_save
Iter = Iter_save
Call mma_deallocate(Aux_b)

!Add some unit vectors correponding to the Krylov subspace algorithm, g, Ag, A^2g, ....
j = j + 1
Aux_a(:) = g(:,nDIIS)
e_diis(:,j) = Aux_a(:) / Sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))
j = j + 1
Aux_a(:) = dq(:)
e_diis(:,j) = Aux_a(:) / Sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))
Call mma_deallocate(Aux_a)

#endif

! Now orthogonalize all unit vectors
#ifdef _DEBUGPRINT_
If (Allocated(e_diis)) Call RecPrt('e_diis(unnorm)',' ',e_diis,mOV,nExplicit)
#endif
Do l = 1, 2
j = 1
Do i = 2, nExplicit
   Do k = 1, j
      gg = DDot_(mOV,e_diis(:,i),1,e_diis(:,k),1)
#ifdef _DEBUGPRINT_
!     Write (6,*) 'i,k,gg=',i,k,gg
#endif
      e_diis(:,i) = e_diis(:,i) - gg * e_diis(:,k)
   End Do
   gg = DDot_(mOV,e_diis(:,i),1,e_diis(:,i),1) ! renormalize
#ifdef _DEBUGPRINT_
   Write (6,*) 'j,i,gg=',j,i,gg
#endif
   If (gg>1.0D-17) Then   ! Skip vector if linear dependent.
      j = j + 1
      e_diis(:,j) = e_diis(:,i)/Sqrt(gg)
   End If
End Do
End Do
mDIIS=j
#ifdef _DEBUGPRINT_
Write (6,*) '      mOV:',mOV
#if defined(_HYBRID_)
Write (6,*) 'nExplicit:',nExplicit,'=', 2*(nDIIS-1), '+',nKrylov
#elif defined(_HYBRID2_)
Write (6,*) 'nExplicit:',nExplicit,'=', 2*nDIIS-1, '+',nKrylov
#else
Write (6,*) 'nExplicit:',nExplicit
#endif
Write (6,*) 'IterSO   :',IterSO
Write (6,*) '    nDIIS:',nDIIS
Write (6,*) '    mDIIS:',mDIIS
#endif

#ifdef _DEBUGPRINT_
Write (6,*) 'Check the ortonormality'
Do i = 1, mDIIS
   Do j = 1, i
      Write (6,*) i,j,DDot_(mOV,e_diis(:,i),1,e_diis(:,j),1)
   End Do
   Write (6,*)
End Do
If (Allocated(e_diis)) Call RecPrt('e_diis',' ',e_diis,mOV,mDIIS)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computed the projected displacement coordinates. Note that the displacements are relative to the last coordinate, q(:,nDIIS).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Call mma_allocate(q_diis,mDIIS,nDIIS+Max_Iter,Label='q_diis')
q_diis(:,:)=Zero
Do i = 1, nDIIS
   Do k = 1, mDIIS
      gg = Zero
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
g_diis(:,:)=Zero
Do i = 1, nDIIS
   Do k = 1, mDIIS
      gg = Zero
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
      gg = Zero
      Do l = 1, mOV
         gg = gg + e_diis(l,i)*HDiag(l)*e_diis(l,j)
      End Do
      H_diis(i,j)=gg
   End Do
   HDiag_Diis(i)=H_Diis(i,i)
End Do
#ifdef _DEBUGPRINT_
Call RecPrt('H_diis(HDiag)',' ',H_diis,mDIIS,mDIIS)
#endif

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

Call mma_allocate(dq_diis,mDiis,Label='dq_Diis')

!We need to set the bias

blavAI=10.00D0
 Call Setup_Kriging(nDiis,mDiis,q_diis,g_diis,Energy(iFirst),Hessian_HMF=H_diis)
!Call Setup_Kriging(nDiis,mDiis,q_diis,g_diis,Energy(iFirst),HDiag=HDiag_diis)
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
Write (6,*) 'Energy(iter)-Energy(iter-1)=', Energy(iter)-Energy(iter-1)
Write (6,*) 'nDIIS=', nDIIS
Write (6,*) 'Beta_Disp_Seed=', Beta_Disp_Seed
Write (6,*) 'Beta_Disp_Min=', Beta_Disp_Min
Write (6,*) 'Beta_Disp=', Beta_Disp
#endif

Do While (.NOT.Converged .and. nDIIS>1) ! Micro iterate on the surrogate model

   Iteration_Micro = Iteration_Micro + 1
   Iteration_Total = Iteration_Total + 1
   Iteration       = Iteration       + 1
   If (Iteration_Micro==Max_Iter) Then
      Write (6,*)
      Write (6,*)  'S-GEK-Optimizer: Iteration_Micro==Max_Iter'
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
   StepMax=StepMax_Seed*DBLE(Iteration_Micro)
   ! Loop to enforce restricted variance. Note, if the step restriction kicks in no problem since we will still microiterate.
   ! Normally a full step will be allowed -- no step restriction -- and the loop will be exited after the first iteration.
   Do

      ! Compute the surrogate Hessian
      Call Hessian_Kriging_Layer(q_diis(:,Iteration),H_diis,mDiis)
     !Call Hessian_Kriging(q_diis(:,Iteration),H_diis,mDiis)

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
      Write (6,*) 'Subiteration: Step_Trunc, StepMax,dqdq:',Step_Trunc,StepMax, dqdq
      Write (6,*) 'Subiteration: Step_Trunc_        :',Step_Trunc_
      Call RecPrt('dq_diis',' ',dq_diis,mDIIS,1)
      Call RecPrt('q_diis(:,Iteration+1)',' ',q_diis(:,Iteration+1),mDIIS,1)
#endif
      If (Step_Trunc.eq.'N') Step_Trunc=' '   ! set to blank if not touched
      If (Step_Trunc//Step_Trunc_==' *') Step_Trunc='.' ! Mark that we have had a step Reduction

      Call Dispersion_Kriging_Layer(q_diis(:,Iteration+1),Variance,mDIIS)
     !Call Dispersion_Kriging(q_diis(:,Iteration+1),Variance,mDIIS)

      ! Note that we might have converged because the step restriction kicked in. However, we fill implicitly
      ! fix that during the second micro iteration.

#ifdef _DEBUGPRINT_
      Write (6,*)
      Write (6,*) 'Step_Trunc:',Step_Trunc
      Write (6,*) 'Beta_Disp =',Beta_Disp
      Write (6,*) 'Variance  =',Variance(1)
      Write (6,*) 'Fact      =',Fact
      Write (6,*) 'StepMax   =',StepMax
#endif
      If (One-Variance(1)/Beta_Disp>1.0D-3) Exit
      If ( (Fact<1.0D-5) .OR. (Variance(1)<Beta_Disp) ) Exit
      Fact   =Half*Fact
      StepMax=Half*StepMax
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

  !Compute the energy and the gradient of the surrogate model

   Call Energy_Kriging_Layer(q_diis(:,Iteration+1),Energy(Iteration_Total+1),mDIIS)
  !Call Energy_Kriging(q_diis(:,Iteration+1),Energy(Iteration_Total+1),mDIIS)
   Call Gradient_Kriging_Layer(q_diis(:,Iteration+1),g_diis(:,Iteration+1),mDIIS)
  !Call Gradient_Kriging(q_diis(:,Iteration+1),g_diis(:,Iteration+1),mDIIS)

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
      Exit
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
Call RecPrt('dq_diis',' ',dq_diis(:),SIZE(dq_diis),1)
Write (6,*) '||dq||=', Sqrt(DDot_(SIZE(dq),dq(:),1,dq(:),1))
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

Call mma_deallocate(h_diis)
Call mma_deallocate(q_diis)
Call mma_deallocate(g_diis)
If (Allocated(e_diis)) Call mma_deallocate(e_diis)

Call mma_deallocate(g)
Call mma_deallocate(q)

#ifdef _DEBUGPRINT_
Write (6,*) 'Exit S-GEK Optimizer'
#endif

End Subroutine S_GEK_Optimizer
