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
* Copyright (C) 1994,2004,2014,2017,2019, Roland Lindh                 *
*               2014,2018, Ignacio Fdez. Galvan                        *
************************************************************************
      Subroutine RS_RFO(H,g,nInter,dq,UpMeth,dqHdq,StepMax,Step_Trunc)
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
*     Remove references to work, Roland Lindh                          *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
      Integer nInter
      Real*8 H(nInter,nInter), g(nInter), dq(nInter)
      Character UpMeth*6, Step_Trunc*1
      Real*8 StepMax
*
*     Local variables
      Real*8, Dimension(:), Allocatable:: Tmp, Val, Matrix
      Real*8, Dimension(:,:), Allocatable:: Vec
      Logical Iterate, Restart
      Real*8 Lambda
*
      UpMeth='RS-RFO'
      Lu=6
*#define _DEBUG_
*define _DEBUG2_
#ifdef _DEBUG_
      Call RecPrt(' In RS_RFO: H',' ',H,nInter,nInter)
      Call RecPrt(' In RS_RFO: g',' ', g,nInter,1)
      Call RecPrt(' In RS_RFO: q',' ', q,nInter,1)
      Write (Lu,*) 'Trust radius=',StepMax
*
      Write (Lu,*)
      Write (Lu,*) 'RS-RF Optimization'
      Write (Lu,*) ' Iter   alpha          dqdq    StepMax     EigVal'
#endif
*
      A_RFO=One   ! Initial seed of alpha
      IterMx=25
      Iter=0
      Iterate=.False.
      Restart=.False.
      Thr=1.0D-7
*ifdef _DEBUG_
*     NumVal=nInter+1
*else
      NumVal=Min(6,nInter)+1
*endif
      Call mma_allocate(Vec,(nInter+1),NumVal,Label='Vec')
      Call mma_allocate(Val,NumVal,Label='Val')
      Call mma_allocate(Matrix,(nInter+1)*(nInter+2)/2,Label='Matrix')
      Call mma_allocate(Tmp,nInter+1,Label='Tmp')
*
      Vec(:,:)=0.0D0
      Tmp(:)=0.0
 998  Continue
         Iter=Iter+1
#ifdef _DEBUG2_
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
*        Set up the augmented Hessian of Eq. (2)
*
*        Assume that the S-matrix is trivially diagonal
*
         Do i = 1, nInter
            Do j = 1, i
               ij = i*(i-1)/2 + j
               Matrix(ij)=Half*(H(i,j)+H(j,i))/A_RFO
            End Do
         End Do
         j=  nInter+1
         Do i = 1, nInter
            ij = j*(j-1)/2 + i
            Matrix(ij)=-g(i)/Sqrt(A_RFO) ! note sign
         End Do
         jj = j*(j+1)/2
         Matrix(jj)=Zero
#ifdef _DEBUG2_
         Call TriPrt('R_Tri',' ',Matrix,nInter+1)
#endif
*
*        Restore the vector from the previous iteration, if any
         call dcopy_(nInter+1,Tmp(1),1,Vec(1,1),1)
         Call Davidson(Matrix,nInter+1,NumVal,Val,Vec,iStatus)
         If (iStatus.gt.0) Then
           Call SysWarnMsg('RS_RFO',
     &       'Davidson procedure did not converge','')
         End If
*
*        Pick up the root which represents the shortest displacement.
*
#ifdef _DEBUG_
*        Call RecPrt('Val',' ',Val,1,NumVal)
*        Call RecPrt('Vec',' ',Vec,nInter+1,NumVal)
#endif
         iRoot=-1
         Dist=1.0D99
         Do iVal = 1, NumVal
            If (Vec(nInter+1,iVal).eq.0.0D0) Cycle
            VV=DDot_(nInter,Vec(1,iVal),1,Vec(1,iVal),1)
            ZZ = VV/A_RFO + Vec(nInter+1,iVal)**2
            Fact=Vec(nInter+1,iVal)/Sqrt(ZZ)
            dqdq=VV/(A_RFO*Fact**2*ZZ)
#ifdef _DEBUG_
*           Write (6,*)
*           Write (6,*) 'iVal,A_RFO=',iVal,A_RFO
*           Write (6,*) 'ZZ=',ZZ
*           Write (6,*) 'Fact=',Fact
*           Write (6,*) 'dqdq=',dqdq
#endif
            If (dqdq.lt.Dist) Then
               iRoot=iVal
               Dist=dqdq
            End If
         End Do
         If (iRoot.eq.-1) Then
            Write (6,*)
            Write (6,*) 'RS-RFO: Illegal iroot value!'
            Call Abend()
         End If
         call dcopy_(nInter+1,Vec(1,iRoot),1,Tmp,1)
         Call DScal_(nInter,One/Sqrt(A_RFO),Vec(1,iRoot),1)
*                                                                      *
************************************************************************
*                                                                      *
*        Execute step 2 on page 266                                    *
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG2_
         Write (Lu,*) ' RF eigenvalue=',Val
#endif
         ZZ=DDot_(nInter+1,Vec(1,iRoot),1,Vec(1,iRoot),1)
         Call DScal_(nInter+1,One/Sqrt(ZZ),Vec(1,iRoot),1)
*                                                                      *
************************************************************************
*                                                                      *
*       Execute step 3 of page 266                                     *
*                                                                      *
************************************************************************
*                                                                      *
*        Copy v^k_{n,i}
*
         call dcopy_(nInter,Vec(1,iRoot),1,dq,1)
*
*        Pick v^k_{1,i}
*
         Fact=Vec(nInter+1,iRoot)
#ifdef _DEBUG2_
         Write (Lu,*) 'v^k_{1,i}=',Fact
#endif
*
*        Normalize according to Eq. (5)
*
         Call DScal_(nInter,One/Fact,dq,1)
#ifdef _DEBUG_
*           Write (6,*)
*           Write (6,*) 'iRoot=',iRoot
*           Write (6,*) 'ZZ=',ZZ
*           Write (6,*) 'Fact=',Fact
*           Write (6,*) 'dqdq=',Sqrt(DDot_(nInter,dq,1,dq,1))
#endif
*
*        Compute lambda_i according to Eq. (8a)
*
         EigVal=-DDot_(nInter,dq,1,g,1) ! note sign
         Lambda = EigVal
*
*        Compute R^2 according to Eq. (8c)
*
         dqdq=Sqrt(DDot_(nInter,dq,1,dq,1))
#ifdef _DEBUG_
         Write (Lu,'(I5,5(E12.5,1x))') Iter,A_RFO,dqdq,StepMax,EigVal
#endif
*                                                                      *
************************************************************************
*                                                                      *
*------- Initialize data for iterative scheme (only at first iteration)
*
         If (.Not.Iterate.Or.Restart) Then
            A_RFO_long=A_RFO
            dqdq_long=dqdq
            A_RFO_short=Zero
            dqdq_short=dqdq_long+One
         End If
*                                                                      *
************************************************************************
*                                                                      *
*------- RF with constraints. Start iteration scheme if computed step
*        is too long.
*
         If ((Iter.eq.1.or.Restart).and.dqdq.gt.StepMax) Then
            Iterate=.True.
            Restart=.False.
         End If
*                                                                      *
************************************************************************
*                                                                      *
*        Procedure if the step length is not equal to the trust radius
*
         If (Iterate.and.Abs(StepMax-dqdq).gt.Thr) Then
            Step_Trunc='*'
#ifdef _DEBUG2_
            Write (Lu,*) 'StepMax-dqdq=',StepMax-dqdq
#endif
*
*           Converge if small interval
*
            If ((dqdq.lt.StepMax).and.
     &          (Abs(A_RFO_long-A_RFO_short).lt.Thr)) Go To 997
            Call Find_RFO_Root(A_RFO_long,dqdq_long,
     &                         A_RFO_short,dqdq_short,
     &                         A_RFO,dqdq,StepMax)
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
#ifdef _DEBUG_
      Write (Lu,*)
      Write (Lu,*) 'Rational Function Optimization, Lambda=',EigVal
      Write (Lu,*)
      Write (Lu,*) 'EigVal,dqHdq=',EigVal,dqHdq
      Call RecPrt(' In RS_RFO: g',' ', g,nInter,1)
      Call RecPrt(' In RS_RFO:dq',' ',dq,nInter,1)
#endif
#define _CHECK_UPDATE_
#ifdef _CHECK_UPDATE_
      Thr_Check=1.0D2
      Do i = 1, nInter
         If (ABS(dq(i)).gt.Thr_Check) Then
            Write (6,*) 'RS_RFO: ABS(dq(i)).gt.Thr_Check'
            Write (6,*) '        Probably an error.'
            Call Abend()
         End If
      End Do
#endif
*
      Call mma_deallocate(Vec)
      Call mma_deallocate(Val)
      Call mma_deallocate(Matrix)
*
      Return
      End
