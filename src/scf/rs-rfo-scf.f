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
* Copyright (C) 1994,2004,2014,2017, Roland Lindh                      *
*               2014,2018, Ignacio Fdez. Galvan                        *
************************************************************************
      Subroutine RS_RFO_SCF(HDiag,g,nInter,dq,UpMeth,dqdq,dqHdq,StepMax,
     &                      Step_Trunc)
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
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
      Integer nInter
      Real*8 HDiag(nInter), g(nInter), dq(nInter)
      Character UpMeth*6, Step_Trunc*1
      Real*8 dqdq, dqHdq, StepMax
*     Local variables
      Real*8, Allocatable:: Tmp(:), Val(:), Vec(:,:)
      Logical Iterate, Restart
*
      UpMeth='RS-RFO'
      Step_Trunc=' '
      Lu=6
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
*     Call RecPrt('rs-rfo: HDiag',' ',HDiag,1,nInter)
*     Call RecPrt('rs-rfo: g',' ',g,1,nInter)
      Write (Lu,*)
      Write (Lu,*) '***************************************************'
      Write (Lu,*) '********* S T A R T  O F  R S - R F O *************'
      Write (Lu,*) '***************************************************'
      Call NrmClc(g,nInter,'RS-RFO','g(n)')
      Write (Lu,*) 'Trust radius=',StepMax
*
      Write (Lu,*)
      Write (Lu,*) 'RS-RF Optimization'
      Write (Lu,*) ' Iter    alpha    Sqrt(dqdq)  StepMax    EigVal'
#endif
*
      A_RFO=One   ! Initial seed of alpha
      IterMx=25
      Iter=0
      Iterate=.False.
      Restart=.False.
      Thr=1.0D-3
      NumVal=Min(1,nInter+1)
!     NumVal=Min(nInter+1,nInter+1)
      Call mma_allocate(Vec,(nInter+1),NumVal,Label='Vec')
      Call mma_allocate(Val,NumVal,Label='Val')
      Call mma_allocate(Tmp,nInter+1,Label='Tmp')
*
      Vec(:,:)=Zero
      Tmp(:)=Zero
 998  Continue
         Iter=Iter+1
#ifdef _DEBUGPRINT_
         Write (Lu,*) 'Iter=',Iter
         Write (Lu,*) 'A_RFO=',A_RFO
#endif
*                                                                      *
************************************************************************
*                                                                      *
*        Execute step 1 of page 266                                    *
*                                                                      *
************************************************************************
*                                                                      *
*        Restore the vector from the previous iteration, if any
         call dcopy_(nInter+1,Tmp,1,Vec(:,1),1)
*
*        Call special Davidson routine which do not require the
*        augmented Hessian to be explicitly expressed by rather will
*        handle the gradient and Hessian part separated. The gradient
*        will be explicit, while the Hessian part will use an approach
*        which computes Hc, where c is a trial vector, from an initial
*        Hessian based on a diagonal approximation and a BFGS update.
*
!#define _DEBUGCode_
#ifdef _DEBUGCode_
         Call Plain_rs_rfo()
#else
         Call Davidson_SCF(HDiag,g,nInter,NumVal,A_RFO,Val,Vec,iStatus)
         If (iStatus.gt.0) Then
            Call SysWarnMsg('RS_RFO',
     &       'Davidson procedure did not converge','')
         End If
#endif
*        Write (6,*) 'Val(:)=',Val(:)
*        Write (6,*) 'Vec(:,1)=',Vec(:,1)
         call dcopy_(nInter+1,Vec(:,1),1,Tmp,1)
         Call DScal_(nInter,One/Sqrt(A_RFO),Vec(:,1),1)
*                                                                      *
************************************************************************
*                                                                      *
*        Execute step 2 on page 266                                    *
*                                                                      *
************************************************************************
*                                                                      *
*        Write (Lu,*) ' RF eigenvalue=',Val
         ZZ=DDot_(nInter+1,Vec(:,1),1,Vec(:,1),1)
         Call DScal_(nInter+1,One/Sqrt(ZZ),Vec(:,1),1)
*                                                                      *
************************************************************************
*                                                                      *
*       Execute step 3 of page 266                                     *
*                                                                      *
************************************************************************
*                                                                      *
*        Copy v^k_{n,i}
*
         call dcopy_(nInter,Vec(:,1),1,dq,1)
*
*        Pick v^k_{1,i}
*
         Fact=Vec(nInter+1,1)
*        Write (Lu,*) 'v^k_{1,i}=',Fact
*
*        Normalize according to Eq. (5)
*
         Call DScal_(nInter,One/Fact,dq,1)
*
*        Compute lambda_i according to Eq. (8a)
*
         EigVal=-DDot_(nInter,dq,1,g,1) ! note sign
*
*        Compute R^2 according to Eq. (8c)
*
         dqdq=DDot_(nInter,dq,1,dq,1)
#ifdef _DEBUGPRINT_
         Write (Lu,'(I5,4E11.3)') Iter,A_RFO,Sqrt(dqdq),StepMax,EigVal
#endif
*                                                                      *
************************************************************************
*                                                                      *
*------- Initialize data for iterative scheme (only at first iteration)
*
         If (.Not.Iterate.Or.Restart) Then
            A_RFO_long=A_RFO
            dqdq_long=Sqrt(dqdq)
            A_RFO_short=Zero
            dqdq_short=dqdq_long+One
         End If
*                                                                      *
************************************************************************
*                                                                      *
*------- RF with constraints. Start iteration scheme if computed step
*        is too long.
*
         If ((Iter.eq.1.or.Restart).and.dqdq.gt.StepMax**2) Then
            Iterate=.True.
            Restart=.False.
         End If
*                                                                      *
************************************************************************
*                                                                      *
*        Procedure if the step length is not equal to the trust radius
*
         If (Iterate.and.Abs(StepMax-Sqrt(dqdq)).gt.Thr) Then
            Step_Trunc='*'
*           Write (Lu,*) 'StepMax-Sqrt(dqdq)=',StepMax-Sqrt(dqdq)
*
*           Converge if small interval
*
            If ((dqdq.lt.StepMax**2).and.
     &          (Abs(A_RFO_long-A_RFO_short).lt.Thr)) Go To 997
            Call Find_RFO_Root(A_RFO_long,dqdq_long,
     &                         A_RFO_short,dqdq_short,
     &                         A_RFO,Sqrt(dqdq),StepMax)
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
*
 997  Continue
      Call mma_deallocate(Tmp)
      dqHdq=dqHdq+EigVal*Half
#ifdef _DEBUGPRINT_
      Write (Lu,*)
      Write (Lu,*) 'Rational Function Optimization, Lambda=',EigVal
      Write (Lu,*)
      Write (Lu,*) 'EigVal,dqHdq=',EigVal,dqHdq
      Call NrmClc(g,nInter,'RS-RFO','g(n)')
      Call NrmClc(dq,nInter,'RS-RFO','dX(n)')
      Write (Lu,*) '***************************************************'
      Write (Lu,*) '************* E N D  O F  R S - R F O *************'
      Write (Lu,*) '***************************************************'
      Write (Lu,*)
#endif
*
      Call mma_deallocate(Vec)
      Call mma_deallocate(Val)
*
#ifdef _DEBUGCode_
      Contains

      Subroutine Plain_rs_rfo()
      use LnkLst, only: SCF_V
#include "mxdm.fh"
#include "infscf.fh"
#include "llists.fh"
      Real*8, Allocatable :: H(:,:) , H_Aug(:,:)
      Real*8, Allocatable :: dq(:), dg(:), Hdq(:)
      Real*8, Allocatable :: EVec(:,:), EVal(:)
      Integer i, j

      Call mma_allocate(H,nInter,nInter,Label='H')
      Call mma_allocate(H_Aug,nInter+1,nInter+1,Label='H')

      H(:,:)=Zero
      Do i = 1, nInter
         H(i,i)=HDiag(i)
      End Do
#ifdef _DEBUGPRINT_
      Call RecPrt('H',' ',H,nInter,nInter)
#endif

!     Update matrix with BFGS

      Call mma_allocate( dq,nInter,Label=' dq')  ! q(:,i)-q(:,i-1)
      Call mma_allocate( dg,nInter,Label=' dg')  ! g(:,i)-g(:,i-1)
      Call mma_allocate(Hdq,nInter,Label='Hdq')  ! Scratch
      Do i = 3, iter
         Write (6,*) 'i,iter=',i,iter

!        delta step
         ipdq=LstPtr(i  ,LLx)
         dq(:)=SCF_V(ipdq)%A(:)
         ipdq=LstPtr(i-1,LLx)
         dq(:)=dq(:)-SCF_V(ipdq)%A(:)

!        delta gradient
         ipdg=LstPtr(i  ,LLGrad)
         dg(:)=SCF_V(ipdg)%A(:)
         ipdg=LstPtr(i-1,LLGrad)
         dg(:)=dg(:)-SCF_V(ipdg)%A(:)

         Call DFP(H,nInter,Hdq,dq,dg)

      End Do
      Call mma_deallocate(Hdq)
      Call mma_deallocate( dg)
      Call mma_deallocate( dq)
#ifdef _DEBUGPRINT_
      Call RecPrt('H(updated)',' ',H,nInter,nInter)
      Call RecPrt('g         ',' ',g,1,nInter)
#endif

!     Set up the augmented Hessian
      H_Aug(:,:)=Zero
      Do i = 1, nInter
         H_Aug(nInter+1,i)=g(i)/Sqrt(A_RFO)
         H_Aug(i,nInter+1)=g(i)/Sqrt(A_RFO)
         Do j = 1, nInter
            H_Aug(i,j)=H(i,j)/A_RFO
         End Do
      End Do
#ifdef _DEBUGPRINT_
      Call RecPrt('H_Aug',' ',H_Aug,nInter+1,nInter+1)
#endif

!     Diagonalize

      nH = nInter+1
      Call mma_allocate(EVec,nH,nH,Label='EVec')
      EVec(:,:)=Zero
      Call mma_allocate(EVal,nH*(nH+1)/2,Label='EVal')

      Do i = 1, nH
         Do j = 1, i
            ij = i*(i-1)/2 + j
            EVal(ij)=H_Aug(i,j)
         End Do
         EVec(i,i)=One
      End Do

      Call NIDiag_new(EVal,EVec,nH,nH)
      Call JacOrd    (EVal,EVec,nH,nH)
#ifdef _DEBUGPRINT_
      Call RecPrt('EVec',' ',EVec,nH,nH)
      Call TriPrt('EVal',' ',EVal,nH)
#endif

      Do i = 1, Size(Val)
         Val(i)=EVal(i*(i+1)/2)
         Vec(:,i)=EVec(:,i)
      End Do

      Call mma_deallocate(EVal)
      Call mma_deallocate(EVec)
      Call mma_deallocate(H_Aug)
      Call mma_deallocate(H)
      End Subroutine Plain_rs_rfo

      Subroutine DFP(B,nDim,Bd,Delta,Gamma)
      Implicit Real*8 (a-h,o-z)
      Real*8 B(nDim,nDim), Bd(nDim),Gamma(nDim),Delta(nDim)
      Real*8, Parameter :: Thr=1.0D-8
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Call RecPrt('DFP: B',' ',B,nDim,nDim)
*     Call RecPrt('DFP: Bd',' ',Bd,1,nDim)
      Call RecPrt('DFP: Gamma',' ',Gamma,1,nDim)
      Call RecPrt('DFP: Delta',' ',Delta,1,nDim)
#endif
      Call DGEMM_('N','N',
     &            nDim,1,nDim,
     &            1.0d0,B,nDim,
     &                  Delta,nDim,
     &            0.0d0,Bd,nDim)
      gd=DDot_(nDim,Gamma,1,Delta,1)
      dBd=DDot_(nDim,Delta,1,Bd,1)
#ifdef _DEBUGPRINT_
      Call RecPrt('DFP: Bd',' ',Bd,1,nDim)
      Write (6,*) 'gd=',gd
      Write (6,*) 'dBd=',dBd
      Write (6,*) 'Thr=',Thr
#endif
      If (gd<0.0D0) Then
         Call MSP(B,Gamma,Delta,nDim)
      Else
*
         Do i = 1, nDim
            Do j = 1, nDim
               B(i,j) = B(i,j) + (Gamma(i)*Gamma(j))/gd
     &                         - (Bd(i)*Bd(j))/dBd
            End Do
         End Do
      End If
*
#ifdef _DEBUGPRINT_
      Call RecPrt('DFP: B',' ',B,nDim,nDim)
#endif
      End Subroutine DFP

      Subroutine MSP(B,Gamma,Delta,nDim)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 B(nDim,nDim),Gamma(nDim),Delta(nDim)
*
*                              T       T            ( T)
*                    |(1-phi)/d g phi/d d|        | (g )
*                    |     T       T    T         | ( T)
*     B = B + (g  d )|phi/d d    -Phi*d g/(d d)**2| (d )
*
*
*
      gd= DDot_(nDim,Gamma,1,Delta,1)
      dd= DDot_(nDim,Delta,1,Delta,1)
      gg= DDot_(nDim,Gamma,1,Gamma,1)
      phi=(One-((gd**2)/(dd*gg)))
      e_msp=(gd/dd)**2*((Two/(One-Phi*Sqrt(Phi)))-One)
#ifdef _DEBUGPRINT_
      Call RecPrt(' MSP: Hessian',' ',B,nDim,nDim)
      Call RecPrt(' MSP: Delta',' ',Delta,nDim,1)
      Call RecPrt(' MSP: Gamma',' ',Gamma,nDim,1)
      Write (6,*) 'MSP: Phi=',Phi
      Write (6,*) 'gd,dd,gg=', gd,dd,gg
      Write (6,*) 'MSP: a=',Sqrt(Phi)
      Write (6,*) 'MSP: E_msp=',E_msp
#endif
      Do i = 1, nDim
         Do j = 1, nDim
            B(i,j) = B(i,j)
     &             + ((One-phi)/gd)*Gamma(i)*Gamma(j)
     &             + phi*( (Gamma(i)*Delta(j)+Delta(i)*Gamma(j))/dd
     &                   - gd*Delta(i)*Delta(j)/dd**2 )
         End Do
      End Do
#ifdef _DEBUGPRINT_
      Call RecPrt(' MSP: Updated Hessian',' ',B,nDim,nDim)
#endif
      Return
      End Subroutine MSP

#endif

      End Subroutine RS_RFO_SCF
