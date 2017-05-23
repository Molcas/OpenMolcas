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
*               2014, Ignacio Fdez. Galvan                             *
************************************************************************
*define _DEBUG_SORUPV_
#ifdef _DEBUG_SORUPV_
      Subroutine RS_RFO_SCF(HDiag,g,nInter,dq,UpMeth,dqdq,dqHdq,StepMax,
     &                      Step_Trunc,MemRsv,iter_SCF)
#else
      Subroutine RS_RFO_SCF(HDiag,g,nInter,dq,UpMeth,dqdq,dqHdq,StepMax,
     &                      Step_Trunc,MemRsv)
#endif
************************************************************************
*                                                                      *
*     Object: Automatic restricted-step rational functional            *
*             optimization.                                            *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             December '94                                             *
*                                                                      *
*     Modified to the restricted-step RFO method of Besalu and Bofill. *
*     Ref: E. Besalu and J. M. Bofill, TCA, 100, 265-274 (1998), by    *
*     R. Lindh, Gyeongju, Korea.                                       *
*     Removed full diagonalizations, Ignacio Fdez. Galvan, Uppsala     *
*     Remove references to work, Roland Lindh, Harvard, Cambridge      *
*     Modified for SCF, Roland Lindh, Harvard, Cambridge               *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
      Integer nInter
      Real*8 HDiag(nInter), g(nInter), dq(nInter)
      Character UpMeth*6, Step_Trunc*1
      Real*8 dqdq, dqHdq, StepMax
#ifdef _DEBUG_SORUPV_
      Logical Update_H
#endif
*     Local variables
      Real*8, Dimension(:), Allocatable:: Tmp, Val, Vec
      Logical Iterate
      Real*8 Lambda
*
#ifdef _DEBUG_SORUPV_
      Update_H = .True.
#endif
      UpMeth='RS-RFO'
      Step_Trunc=' '
      Lu=6
*define _DEBUG_
#ifdef _DEBUG_
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
      Thr=1.0D-3
      NumVal=1
      Call mma_allocate(Vec,(nInter+1)*NumVal,Label='Vec')
      Call mma_allocate(Val,NumVal,Label='Val')
      Call mma_allocate(Tmp,nInter+1,Label='Tmp')
*
      Call DZero(Vec,(nInter+1)*NumVal)
      Call DZero(Tmp,nInter+1)
 998  Continue
         Iter=Iter+1
#ifdef _DEBUG_
*        Write (Lu,*) 'Iter=',Iter
*        Write (Lu,*) 'A_RFO=',A_RFO
#endif
*                                                                      *
************************************************************************
*                                                                      *
*        Execute step 1 of page 266                                    *
*                                                                      *
************************************************************************
*                                                                      *
*        Restore the vector from the previous iteration, if any
         call dcopy_(nInter+1,Tmp,1,Vec,1)
*
*        Call special Davidson routine which do not require the
*        augmented Hessian to be explicitly expressed by rather will
*        handle the gradient and Hessian part separated. The gradient
*        will be explicit, while the Hessian part will use an approach
*        which computes Hc, where c is a trial vector, from the use
*        an Hessian based on a diagonal approximation and a BFGS update.
*
*        UNDER DEVELOPMENT
*
#ifdef _DEBUG_SORUPV_
         Call Davidson_SCF(HDiag,g,nInter,NumVal,A_RFO,Val,Vec,MemRsv,
     &                     iStatus,iter_SCF,Update_H)
#else
         Call Davidson_SCF(HDiag,g,nInter,NumVal,A_RFO,Val,Vec,MemRsv,
     &                     iStatus)
#endif
*
*        END DEVELOPMENT
*
         If (iStatus.gt.0) Then
            Call SysWarnMsg('RS_RFO',
     &       'Davidson procedure did not converge','')
         End If
         call dcopy_(nInter+1,Vec,1,Tmp,1)
         Call DScal_(nInter,One/Sqrt(A_RFO),Vec,1)
*                                                                      *
************************************************************************
*                                                                      *
*        Execute step 2 on page 266                                    *
*                                                                      *
************************************************************************
*                                                                      *
*        Write (Lu,*) ' RF eigenvalue=',Val
         ZZ=DDot_(nInter+1,Vec,1,Vec,1)
         Call DScal_(nInter+1,One/Sqrt(ZZ),Vec,1)
*                                                                      *
************************************************************************
*                                                                      *
*       Execute step 3 of page 266                                     *
*                                                                      *
************************************************************************
*                                                                      *
*        Copy v^k_{n,i}
*
         call dcopy_(nInter,Vec,1,dq,1)
*
*        Pick v^k_{1,i}
*
         Fact=Vec(nInter+1)
*        Write (Lu,*) 'v^k_{1,i}=',Fact
*
*        Normalize according to Eq. (5)
*
         Call DScal_(nInter,One/Fact,dq,1)
*
*        Compute lambda_i according to Eq. (8a)
*
         EigVal=-DDot_(nInter,dq,1,g,1) ! note sign
         Lambda = EigVal
*
*        Compute R^2 according to Eq. (8c)
*
         dqdq=DDot_(nInter,dq,1,dq,1)
#ifdef _DEBUG_
         Write (Lu,'(I5,4E11.3)') Iter,A_RFO,Sqrt(dqdq),StepMax,EigVal
#endif
*                                                                      *
************************************************************************
*                                                                      *
*------- Initialize data for iterative scheme (only at first iteration)
*
         If (.Not.Iterate) Then
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
         If (Iter.eq.1.and.dqdq.gt.StepMax**2) Iterate=.True.
*                                                                      *
************************************************************************
*                                                                      *
*        Procedure if the step length is not equal to the trust radius
*
         If (Iterate.and.Abs(StepMax-Sqrt(dqdq)).gt.Thr) Then
            Step_Trunc='*'
*           Write (Lu,*) 'StepMax-Sqrt(dqdq)=',StepMax-Sqrt(dqdq)
            Call Find_RFO_Root(A_RFO_long,dqdq_long,
     &                         A_RFO_short,dqdq_short,
     &                         A_RFO,Sqrt(dqdq),StepMax)
            If (Iter.gt.IterMx) Then
               Write (Lu,*) ' Too many iterations in RF'
               Call Abend()
*              Go To 997
            End If
            Go To 998
         End If
*
*997  Continue
      Call mma_deallocate(Tmp)
      dqHdq=dqHdq+EigVal*Half
#ifdef _DEBUG_
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
      Return
      End
