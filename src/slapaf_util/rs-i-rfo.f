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
* Copyright (C) 1994,1997, Roland Lindh                                *
************************************************************************
      Subroutine RS_I_RFO(H,g,nInter,dq,UpMeth,dqHdq,StepMax,
     &                    Step_Trunc,Thr_RS)
************************************************************************
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             December '94                                             *
*                                                                      *
*                                                                      *
*             Solve      | H    g | | d |     | 1 0 | | d |            *
*                        |  T     | |   | = e |     | |   |            *
*                        | g    0 | | 1 |     | 0 1 | | 1 |            *
*                                                                      *
*             this corresponds to                                      *
*                                                                      *
*             H d + g = e d                                            *
*                                                                      *
*             and                                                      *
*                                                                      *
*               T                                                      *
*             g  d = e                                                 *
*                                                                      *
*             Modified from single negative eigenvalue to an arbitrary *
*             number, June '97, R. Lindh                               *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
      Real*8 H(nInter,nInter), g(nInter), dq(nInter)
      Real*8, Allocatable:: Val(:), Tmp(:,:), Vec(:,:), Mat(:)
*
      Character*6 UpMeth
      Character*1 Step_Trunc
      Logical Found
*
      Lu =6
*#define _DEBUG_
#ifdef _DEBUG_
      Call RecPrt(' In RS_I_RFO: H','(10f10.6)',H,nInter,nInter)
      Call RecPrt(' In RS_I_RFO: g','(10f10.6)', g,nInter,1)
      Call RecPrt(' In RS_I_RFO:dq','(10f10.6)',dq,nInter,1)
#endif
*
      NumVal=Min(2,nInter)
      nVStep=2
      Found=.False.
      Thr=1.0D-6
      Call mma_allocate(Vec,nInter,NumVal,Label='Vec')
      Vec(:,:)=Zero
      Call mma_allocate(Val,NumVal,Label='Val')
      Val(:)=Zero
      Call mma_allocate(Mat,nInter*(nInter+1)/2,Label='Mat')
      Do i = 1, nInter
         Do j = 1, i
            ij = i*(i-1)/2+j
            Mat(ij)=H(i,j)
         End Do
      End Do
*
*---- Find the negative eigenvalue(s)
*     Stop when the highest eigenvalue found is larger than Thr
      Do While (.Not.Found)
        Call Davidson(Mat,nInter,NumVal,Val,Vec,iStatus)
        If (iStatus.gt.0) Then
          Call SysWarnMsg('RS_I_RFO',
     &      'Davidson procedure did not converge','')
        End If
        nNeg=0
        Do i = 1, NumVal
           If (Val(i).lt.Zero) nNeg=nNeg+1
        End Do
        If ( (Val(NumVal).gt.Thr .and. nNeg.gt.0) .or.
     &      (NumVal.ge.nInter)
     &     ) Then
          Found=.True.
        Else
*----     Increase the number of eigenpairs to compute
          Call mma_allocate(Tmp,nInter,NumVal,Label='Tmp')
          call dcopy_(NumVal*nInter,Vec,1,Tmp,1)
          Call mma_deallocate(Vec)
          Call mma_deallocate(Val)
*
          NumVal=Min(NumVal+nVStep,nInter)
*
          Call mma_allocate(Vec,nInter,NumVal,Label='Vec')
          Call mma_allocate(Val,NumVal,Label='Val')
          Vec(:,:)=Zero
          Vec(:,1:NumVal-nVStep) = Tmp(:,:)
          Val(:)=Zero
          Call mma_deallocate(Tmp)
        End If
      End Do
      Call mma_deallocate(Mat)
*
      nNeg=0
      i=NumVal
      Do While ((i.ge.0).and.(nNeg.eq.0))
         If (Val(i).lt.Zero) nNeg=i
         i=i-1
      End Do
      If (nNeg.eq.0) Then
         Write (Lu,*) 'Warning RS-I-RFO: Neq.eq.0'
         Call RecPrt(' In RS_I_RFO: Eigenvalues',' ',Val,1,NumVal)
*        Call Abend()
      End If
#ifdef _DEBUG_
      Call RecPrt(' In RS_I_RFO: Eigenvalues',' ',Val,1,NumVal)
      Call RecPrt(' In RS_I_RFO: Eigenvectors',' ',Vec,nInter,NumVal)
      Write (Lu,*) ' nNeg=',nNeg
#endif
*
*     Transform the gradient and Hessian to generate the
*     corresponding entities for the image function. This
*     corresponds to an elementary Householder orthogonal
*     transformation.
*
      If (nNeg.gt.0) Then
         Call mma_allocate(Tmp,nInter,1,Label='Tmp')
         call dcopy_(nInter,g,1,Tmp(1,1),1)
*
         Do iNeg=1,nNeg
           gi = DDot_(nInter,g,1,Vec(1,iNeg),1)
           Call DaXpY_(nInter,-Two*gi,Vec(1,iNeg),1,g,1)
           Fact = Two * Val(iNeg)
           Do j = 1, nInter
              Do i = 1, nInter
                 H(i,j) = H(i,j) - Fact * Vec(i,iNeg) * Vec(j,iNeg)
              End Do
           End Do
         End Do
      End If
*
      Call mma_deallocate(Vec)
      Call mma_deallocate(Val)
*
      Call RS_RFO(H,g,nInter,dq,UpMeth,dqHdq,StepMax,Step_Trunc,
     &            Thr_RS)
*
*     Restore the original gradient
*
      If (nNeg.gt.0) Then
         call dcopy_(nInter,Tmp(1,1),1,g,1)
         Call mma_deallocate(Tmp)
      End If
*
      UpMeth='RSIRFO'
*
#ifdef _DEBUG_
      Call RecPrt(' In RS_I_RFO: g','(10f10.6)', g,nInter,1)
      Call RecPrt(' In RS_I_RFO:dq','(10f10.6)',dq,nInter,1)
#endif
*
      Return
      End
