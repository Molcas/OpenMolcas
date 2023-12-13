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
! Copyright (C) 1994,2004,2014,2017, Roland Lindh                      *
!               2014,2018, Ignacio Fdez. Galvan                        *
!***********************************************************************
!#define _DEBUGPRINT_
      Subroutine RS_RFO_SCF(g,nInter,dq,UpMeth,dqdq,dqHdq,StepMax_Seed,Step_Trunc)
!***********************************************************************
!                                                                      *
!     Object: Automatic restricted-step rational functional            *
!             optimization.                                            *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             December '94                                             *
!                                                                      *
!     Modified to the restricted-step RFO method of Besalu and Bofill. *
!     Ref: E. Besalu and J. M. Bofill, TCA, 100, 265-274 (1998), by    *
!     R. Lindh, Gyeongju, Korea.                                       *
!     Removed full diagonalizations, Ignacio Fdez. Galvan, Uppsala     *
!     Remove references to work, Roland Lindh, Harvard, Cambridge      *
!     Modified for SCF, Roland Lindh, Harvard, Cambridge               *
!***********************************************************************
      use stdalloc, only: mma_allocate, mma_deallocate
      use Constants, only: Zero, Half, One, Three, Pi
      use InfSCF, only: kOptim, Iter_x=>Iter, Iter_Start
      use InfSO, only: IterSO
      Implicit None
      Integer nInter
      Real*8 g(nInter), dq(nInter)
      Character UpMeth*6
      Real*8 dqdq, dqHdq, StepMax_Seed
      Character Step_Trunc*1

!     Local variables
      Integer :: Lu, IterMx, NumVal, iStatus, iRoot, I, Iter
      Real*8 :: GG, A_RFO, ZZ, Test, Fact, EigVal
      Real*8 :: A_RFO_Long, A_RFO_Short
      Real*8 :: DqDq_Long, DqDq_Short
      Real*8, External :: DDot_
      Real*8, Allocatable:: Tmp(:), Val(:), Vec(:,:)
      Logical Iterate, Restart
      Real*8, Save :: StepMax=One
      Real*8, Save :: Step_Lasttime=Pi
      Real*8, Parameter :: Thr=1.0D-4
      Real*8, Parameter :: StepMax_Min=1.0D-2
      Real*8, Parameter :: Step_Factor=Three
!
      UpMeth='RS-RFO'
      Step_Trunc=' '
      Lu=6
      gg=Sqrt(DDot_(nInter,g,1,g,1))
#ifdef _DEBUGPRINT_
      Write (6,*) 'StepMax_Seed=',StepMax_Seed
      Write (6,*) 'Sqrt(gg)=',gg
#endif


      If (Step_Lasttime==Pi) Step_Lasttime=Step_Lasttime/(gg*Step_Factor)

      StepMax=Min(Pi,StepMax_Seed*gg,Step_Lasttime*Step_Factor*gg)

!     Make sure that step restriction is not too tight.
      If (StepMax<StepMax_Min) StepMax=StepMax_Min
#ifdef _DEBUGPRINT_
      Write (6,*) 'StepMax=',StepMax
#endif

#ifdef _DEBUGPRINT_
      Write (Lu,*)
      Write (Lu,*) '***************************************************'
      Write (Lu,*) '********* S T A R T  O F  R S - R F O SCF *********'
      Write (Lu,*) '***************************************************'
      Call NrmClc(g,nInter,'RS-RFO','g(n)')
      Write (Lu,*) 'Trust radius=',StepMax
!
      Write (Lu,*)
      Write (Lu,*) 'RS-RF Optimization'
      Write (Lu,*) ' Iter    alpha    Sqrt(dqdq)  StepMax    EigVal'
#endif
!
      A_RFO=One   ! Initial seed of alpha
      IterMx=50
      Iter=0
      Iterate=.False.
      Restart=.False.
      NumVal=Min(3,nInter+1)
!     NumVal=Min(nInter+1,nInter+1)
      Call mma_allocate(Vec,(nInter+1),NumVal,Label='Vec')
      Call mma_allocate(Val,NumVal,Label='Val')
      Call mma_allocate(Tmp,nInter+1,Label='Tmp')
!
      Vec(:,:)=Zero
      Tmp(:)=Zero
 998  Continue
         Iter=Iter+1
!                                                                      *
!***********************************************************************
!                                                                      *
!        Execute step 1 of page 266                                    *
!                                                                      *
!***********************************************************************
!                                                                      *
!        Restore the vector from the previous iteration, if any
         call dcopy_(nInter+1,Tmp,1,Vec(:,1),1)
!
!        Call special Davidson routine which do not require the
!        augmented Hessian to be explicitly expressed but rather will
!        handle the gradient and Hessian part separated. The gradient
!        will be explicit, while the Hessian part will use an approach
!        which computes Hc, where c is a trial vector, from an initial
!        Hessian based on a diagonal approximation and a BFGS update.
!
         Call Davidson_SCF(g,nInter,NumVal,A_RFO,Val,Vec,iStatus)
         If (iStatus.gt.0) Then
            Call SysWarnMsg('RS_RFO SCF','Davidson procedure did not converge','')
         End If
!        Write (6,*) 'Val(:)=',Val(:)
!        Write (6,*) 'Vec(:,1)=',Vec(:,1)
!        Write (6,*) 'Vec(nInter+1,1)=',Vec(nInter+1,1)

!        Select a root with a negative value close to the current point

         iRoot=1
         dqdq=1.0D10
         Do i = 1, NumVal
            If (Val(i)<Zero) Then
               Tmp(:)=Vec(:,i)/Sqrt(A_RFO)
               ZZ=DDot_(nInter+1,Tmp(:),1,Tmp(:),1)
               Tmp(:)=Tmp(:)/Sqrt(ZZ)
               Tmp(:nInter)=Tmp(:nInter)/Tmp(nInter+1)
               Test=DDot_(nInter,Tmp,1,Tmp,1)
!              Write (6,*) 'Test=',Test
               If (Test<dqdq) Then
                  iRoot = i
                  dqdq = Test
               End If
            End If
         End Do
!        Write (6,*) 'iRoot,dqdq=',iRoot,dqdq
         If (iRoot/=1) Vec(:,1)=Vec(:,iRoot)
         call dcopy_(nInter+1,Vec(:,1),1,Tmp,1)
         Call DScal_(nInter,One/Sqrt(A_RFO),Vec(:,1),1)
!                                                                      *
!***********************************************************************
!                                                                      *
!        Execute step 2 on page 266                                    *
!                                                                      *
!***********************************************************************
!                                                                      *
!        Write (Lu,*) ' RF eigenvalue=',Val
         ZZ=DDot_(nInter+1,Vec(:,1),1,Vec(:,1),1)
         Call DScal_(nInter+1,One/Sqrt(ZZ),Vec(:,1),1)
!                                                                      *
!***********************************************************************
!                                                                      *
!       Execute step 3 of page 266                                     *
!                                                                      *
!***********************************************************************
!                                                                      *
!        Copy v^k_{n,i}
!
         call dcopy_(nInter,Vec(:,1),1,dq,1)
!
!        Pick v^k_{1,i}
!
         Fact=Vec(nInter+1,1)
!        Write (Lu,*) 'v^k_{1,i}=',Fact
!
!        Normalize according to Eq. (5)
!
         Call DScal_(nInter,One/Fact,dq,1)
!
!        Compute lambda_i according to Eq. (8a)
!
         EigVal=-DDot_(nInter,dq,1,g,1) ! note sign
!
!        Compute R^2 according to Eq. (8c)
!
         dqdq=DDot_(nInter,dq,1,dq,1)

         If (Sqrt(dqdq)>Pi) Then
!        If (Sqrt(dqdq)>Pi .or. Sqrt(dqdq)>StepMax.and.kOptim>1) Then
            If (kOptim/=1) Then
               Write (Lu,*) 'rs_rfo_SCF: Total displacement is too large.'
               Write (Lu,*) 'DD=',Sqrt(dqdq)
               Write (Lu,*)'Reset update depth in BFGS, redo the RS-RFO'
               Iter=Iter-1
               kOptim=1
               Iter_Start=Iter_x
               IterSO=1
               Go To 998
            End If
         End If
#ifdef _DEBUGPRINT_
         Write (Lu,'(I5,4ES11.3)') Iter,A_RFO,Sqrt(dqdq),StepMax,EigVal
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
!------- Initialize data for iterative scheme (only at first iteration)
!
         If (.Not.Iterate.Or.Restart) Then
            A_RFO_long=A_RFO
            dqdq_long=Sqrt(dqdq)
            A_RFO_short=Zero
            dqdq_short=dqdq_long+One
         End If
!        Write (6,*) 'dqdq_long=',dqdq_long
!        Write (6,*) 'dqdq_short=',dqdq_short
!                                                                      *
!***********************************************************************
!                                                                      *
!------- RF with constraints. Start iteration scheme if computed step
!        is too long.
!
         If ((Iter==1.or.Restart).and.dqdq>StepMax**2) Then
            Iterate=.True.
            Restart=.False.
         End If
!                                                                      *
!***********************************************************************
!                                                                      *
!        Procedure if the step length is not equal to the trust radius
!
         If (Iterate.and.Abs(StepMax-Sqrt(dqdq))>Thr) Then
            Step_Trunc='*'
!           Write (Lu,*) 'StepMax-Sqrt(dqdq)=',StepMax-Sqrt(dqdq)
!
!           Converge if small interval
!
            If ((dqdq.lt.StepMax**2).and.(Abs(A_RFO_long-A_RFO_short).lt.Thr)) Go To 997
            Call Find_RFO_Root(A_RFO_long,dqdq_long,       &
                               A_RFO_short,dqdq_short,     &
                               A_RFO,Sqrt(dqdq),StepMax)
!           Write (6,*) 'A_RFO_Short=',A_RFO_Short
!           Write (6,*) 'A_RFO_Long=',A_RFO_Long
!           Write (6,*) 'dqdq_long=',dqdq_long
!           Write (6,*) 'dqdq_short=',dqdq_short
            If (A_RFO.eq.-One) Then
               A_RFO=One
               Step_Trunc=' '
               Restart=.True.
               Iterate=.False.
            End If
            If (Iter.gt.IterMx) Then
               Write (Lu,*) ' Too many iterations in RF'
               Go To 997
            End If
            Go To 998
         End If
!
 997  Continue
      Call mma_deallocate(Tmp)
      dqHdq=dqHdq+EigVal*Half
      Step_Lasttime=Sqrt(dqdq)/gg
#ifdef _DEBUGPRINT_
      Write (Lu,*)
      Write (Lu,*) 'Rational Function Optimization, Lambda=',EigVal
      Write (Lu,*)
      Write (Lu,*) 'EigVal,dqHdq=',EigVal,dqHdq
      Call NrmClc(g,nInter,'RS-RFO','g(n)')
      Call NrmClc(dq,nInter,'RS-RFO','dX(n)')
      Write (Lu,*) '***************************************************'
      Write (Lu,*) '************* E N D  O F  R S - R F O SCF *********'
      Write (Lu,*) '***************************************************'
      Write (Lu,*)
#endif
      Call mma_deallocate(Vec)
      Call mma_deallocate(Val)
      End Subroutine RS_RFO_SCF
