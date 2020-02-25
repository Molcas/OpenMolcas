************************************************************************
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
*               2014, Ignacio Fdez. Galvan                             *
************************************************************************
      Subroutine RS_P_RFO(H,g,nInter,dq,UpMeth,dqHdq,StepMax,Step_Trunc)
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
*                                                                      *
*             Removed full diagonalizations, April '14, I. Fdez. Galvan*
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 H(nInter,nInter), g(nInter), dq(nInter), Lambda
*
      Character*6 UpMeth
      Character*1 Step_Trunc
      Logical Found,Iterate
*
*     Call QEnter('RS_P_RFO')
      iRout = 215
      Lu =6
      iPrint = nPrint(iRout)
      If (iPrint.ge.99) Then
         Call RecPrt(' In RS_P_RFO: H','(10f10.6)',H,nInter,nInter)
         Call RecPrt(' In RS_P_RFO: g','(10f10.6)', g,nInter,1)
         Call RecPrt(' In RS_P_RFO:dq','(10f10.6)',dq,nInter,1)
      End If
*
      UpMeth='RSPRFO'
*
      NumVal=Min(2,nInter)
      nVStep=2
      Found=.False.
      Thr=1.0D-6
      Call GetMem('Vector','Allo','Real',ipVec,nInter*NumVal)
      Call GetMem('Values','Allo','Real',ipVal,NumVal)
      Call GetMem('Matrix','Allo','Real',ipMat,nInter*(nInter+1)/2)
      Call DZero(Work(ipVec),NumVal*nInter)
      Do i = 1, nInter
         Do j = 1, i
            ij = i*(i-1)/2+j-1
            Work(ipMat+ij)=H(i,j)
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*---- Find the negative eigenvalue(s)
*     Stop when the highest eigenvalue found is larger than Thr
      Do While (.Not.Found)
        Call Davidson(Work(ipMat),nInter,NumVal,
     &                Work(ipVal),Work(ipVec),iStatus)
        If (iStatus.gt.0) Then
          Call SysWarnMsg('RS_P_RFO',
     &      'Davidson procedure did not converge','')
        End If
        If ((Work(ipVal+NumVal-1).gt.Thr).or.(NumVal.ge.nInter)) Then
          Found=.True.
        Else
*----     Increase the number of eigenpairs to compute
          Call Allocate_Work(ipTmp,NumVal*nInter)
          call dcopy_(NumVal*nInter,Work(ipVec),1,Work(ipTmp),1)
          Call GetMem('Vector','Free','Real',ipVec,nInter*NumVal)
          Call GetMem('Values','Free','Real',ipVal,NumVal)
          NumVal=Min(NumVal+nVStep,nInter)
          Call GetMem('Vector','Allo','Real',ipVec,nInter*NumVal)
          Call GetMem('Values','Allo','Real',ipVal,NumVal)
         call dcopy_((NumVal-nVStep)*nInter,Work(ipTmp),1,Work(ipVec),1)
          Call DZero(Work(ipVec+(NumVal-nVStep)*nInter),nVStep*nInter)
          Call DZero(Work(ipVal),NumVal)
          Call Free_Work(ipTmp)
        End If
      End Do
      Call GetMem('Matrix','Free','Real',ipMat,nInter*(nInter+1)/2)
*
      nNeg=0
      i=NumVal
      Do While ((i.ge.0).and.(nNeg.eq.0))
         i=i-1
         If (Work(ipVal+i).lt.Zero) nNeg=i+1
      End Do
      If (iPrint.ge.99) Then
         Call RecPrt(' In RS_P_RFO: Eigenvalues',' ',Work(ipVal),
     &               1,NumVal)
         Call RecPrt(' In RS_P_RFO: Eigenvectors',' ',Work(ipVec),
     &               nInter,NumVal)
         Write (Lu,*) ' nNeg=',nNeg
      End If
*
      If (iPrint.ge.6) Then
         Write (Lu,*)
         Write (Lu,*) 'RS-P-RF Optimization'
         Write (Lu,*) ' Iter   alpha        dqdq  StepMax'//
     &                '   EigVal_r  EigVal_t'
      End If
*
*     Write (Lu,*) 'Trust radius=',StepMax
      A_RFO=One   ! Initial seed of alpha
      IterMx=25
      Iter=0
      Iterate=.False.
      Thr=1.0D-7
      If (nNeg.gt.0) Then
         mInter=nNeg+1
         Call GetMem('StepN','Allo','Real',ipNStep,nInter)
         Call GetMem('GradN','Allo','Real',ipNGrad,nInter)
         Call GetMem('VectorN','Allo','Real',ipNVec,mInter)
         Call GetMem('ValuesN','Allo','Real',ipNVal,1)
         Call GetMem('MatrixN','Allo','Real',ipNMat,mInter*(mInter+1)/2)
         Call Allocate_Work(ipNTmp,mInter)
         Call DZero(Work(ipNTmp),mInter)
      End If
      mInter=nInter+1
      Call GetMem('StepP','Allo','Real',ipPStep,nInter)
      Call GetMem('GradP','Allo','Real',ipPGrad,nInter)
      Call GetMem('VectorP','Allo','Real',ipPVec,mInter)
      Call GetMem('ValuesP','Allo','Real',ipPVal,1)
      Call GetMem('MatrixP','Allo','Real',ipPMat,mInter*(mInter+1)/2)
      Call Allocate_Work(ipPTmp,mInter)
      Call DZero(Work(ipPTmp),mInter)
 998  Continue
         Iter=Iter+1
*        Write (Lu,*) 'Iter=',Iter
*        Write (Lu,*) 'A_RFO=',A_RFO
         Call DZero(dq,nInter)
         If (nNeg.gt.0) Then
*           write(Lu,*)
*           write(Lu,*) 'Process negative eigenvalues.'
*           write(Lu,*)
            mInter=nNeg+1
*                                                                      *
************************************************************************
*                                                                      *
*--------   Build the augmented matrix of the "negative" subspace
*            diagonal: negative eigenvalues (divided by alpha)
*            last row/column: components of the gradient along the
*                     negative eigenvectors (divided by sqrt(alpha))
*           Project the gradient in the negative subspace (but expressed
*           in the full space)
*           (Note that, since we are interested in the largest eigenvalue,
*            the whole matrix is multiplied by -1, so we actually find the
*            smallest eigenvalue and change its sign)
            Call DZero(Work(ipNGrad),nInter)
            Call DZero(Work(ipNMat),mInter*(mInter+1)/2)
            j=mInter*(mInter-1)/2
            Do i=1,nNeg
               Work(ipNMat+i*(i+1)/2-1)=-Work(ipVal+i-1)/A_RFO
               gv=DDot_(nInter,g,1,Work(ipVec+(i-1)*nInter),1)
               Work(ipNMat+j+i-1)=gv/Sqrt(A_RFO)
               Call DaXpY_(nInter,gv,Work(ipVec+(i-1)*nInter),1,
     &                              Work(ipNGrad),1)
            End Do
*
*--------   Solve the partial RFO system for the negative subspace
            call dcopy_(mInter,Work(ipNTmp),1,Work(ipNVec),1)
            Call Davidson(Work(ipNMat),mInter,1,Work(ipNVal),
     &                                          Work(ipNVec),iStatus)
            call dcopy_(mInter,Work(ipNVec),1,Work(ipNTmp),1)
            If (iStatus.gt.0) Then
               Call SysWarnMsg('RS_P_RFO',
     &              'Davidson procedure did not converge','')
            End If
            Call DScal_(1,-One,Work(ipNVal),1)
*
*--------   Scale the eigenvector (combines eqs. (5) and (23))
*           Convert to full space and add to complete step
            Call DScal_(nNeg,One/(Sqrt(A_RFO)*Work(ipNVec+nNeg)),
     &                      Work(ipNVec),1)
            Call dGeMV_('N',nInter,nNeg,One,Work(ipVec),nInter,
     &                                     Work(ipNVec),1,
     &                                 Zero,Work(ipNStep),1)
            Call DaXpY_(nInter,One,Work(ipNStep),1,dq,1)
            dqdq_max=Sqrt(DDot_(nInter,Work(ipNStep),1,Work(ipNStep),1))
*           write (Lu,*) 'dqdq_max=',dqdq_max
!           Sign
            EigVal_r=-DDot_(nInter,Work(ipNStep),1,Work(ipNGrad),1)
            If (iPrint.ge.99) Then
               Call RecPrt('dq_r',' ',Work(ipNStep),1,nInter)
               Call RecPrt(' g_r',' ',Work(ipNGrad),1,nInter)
               Write (Lu,*) 'Lambda=',EigVal_r
            End If
            If (EigVal_r.lt.-Thr) Then
               Write (Lu,*)
               Write (Lu,*) 'W A R N I N G !'
               Write (Lu,*) 'EigVal_r.lt.Zero',EigVal_r
               Write (Lu,*)
            End If
         Else
            EigVal_r=Zero
            dqdq_max=Zero
         End If
*
*        write(Lu,*)
*        write(Lu,*) 'Process positive eigenvalues.'
*        write(Lu,*)
         mInter=nInter+1
*                                                                      *
************************************************************************
*                                                                      *
*----    Build the augmented matrix of the "positive" subspace
*        Instead of reducing the dimensions, the negative eigenvectors
*        are simply projected out from the gradient and the eigenvalues
*        are shifted to a large positive arbitrary value (10), to avoid
*        interferences
         call dcopy_(nInter,g,1,Work(ipPGrad),1)
         Do i=1,nNeg
           gv=DDot_(nInter,Work(ipPGrad),1,Work(ipVec+(i-1)*nInter),1)
           Call DaXpY_(nInter,-gv,Work(ipVec+(i-1)*nInter),1,
     &                           Work(ipPGrad),1)
         End Do
         Do j=1,nInter
           call dcopy_(j,H(1,j),1,Work(ipPMat+j*(j-1)/2),1)
           Do i=1,nNeg
             ii=(i-1)*nInter-1
             Do k=1,j
               jk=j*(j-1)/2+k-1
               Work(ipPMat+jk)=Work(ipPMat+jk)-(Work(ipVal+i-1)-Ten)*
     &                           Work(ipVec+ii+j)*Work(ipVec+ii+k)
             End Do
           End Do
           Call DScal_(j,One/A_RFO,Work(ipPMat+j*(j-1)/2),1)
         End Do
         Call DZero(Work(ipPMat+mInter*(mInter-1)/2),mInter)
         Call DaXpY_(nInter,-One/Sqrt(A_RFO),Work(ipPGrad),1,
     &                     Work(ipPMat+mInter*(mInter-1)/2),1)
*
*----    Solve the partial RFO system for the positive subspace
         call dcopy_(mInter,Work(ipPTmp),1,Work(ipPVec),1)
         Call Davidson(Work(ipPMat),mInter,1,Work(ipPVal),Work(ipPVec),
     &                 iStatus)
         If (iStatus.gt.0) Then
           Call SysWarnMsg('RS_P_RFO',
     &          'Davidson procedure did not converge','')
         End If
         call dcopy_(mInter,Work(ipPVec),1,Work(ipPTmp),1)
         call dcopy_(nInter,Work(ipPVec),1,Work(ipPStep),1)
*
*----    Scale the eigenvector (combines eqs. (5) and (23))
*        Add to complete step
         Call DScal_(nInter,One/(Sqrt(A_RFO)*Work(ipPVec+nInter)),
     &                     Work(ipPStep),1)
         Call DaXpY_(nInter,One,Work(ipPStep),1,dq,1)
         dqdq_min=Sqrt(DDot_(nInter,Work(ipPStep),1,Work(ipPStep),1))
*        write (Lu,*) 'dqdq_min=',dqdq_min
         EigVal_t=-DDot_(nInter,Work(ipPStep),1,Work(ipPGrad),1) ! Sign
         If (iPrint.ge.99) Then
           Call RecPrt('dq_t',' ',Work(ipPStep),1,nInter)
           Call RecPrt(' g_t',' ',Work(ipPGrad),1,nInter)
           Write (Lu,*) 'Lambda=',EigVal_t
         End If
         If (EigVal_t.gt.Thr) Then
           Write (Lu,*)
           Write (Lu,*) 'W A R N I N G !'
           Write (Lu,*) 'EigVal_t.gt.Zero',EigVal_t
           Write (Lu,*)
         End If
*
      Lambda = EigVal_t + EigVal_r
      dqdq=Sqrt(DDot_(nInter,dq,1,dq,1))
*
      If (iPrint.ge.6)
     &  Write (Lu,'(I5,5F10.5)') Iter,A_RFO,dqdq,StepMax,
     &                           EigVal_r,EigVal_t
*                                                                      *
************************************************************************
*                                                                      *
*------- Initialize data for iterative scheme (only at first iteration)
*
         If (.Not.Iterate) Then
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
         If (Iter.eq.1.and.dqdq.gt.StepMax) Iterate=.True.
*                                                                      *
************************************************************************
*                                                                      *
*        Procedure if the step length is not equal to the trust radius
*
         If (Iterate.and.Abs(StepMax-dqdq).gt.Thr) Then
            Step_Trunc='*'
*           Write (Lu,*) 'StepMax-dqdq=',StepMax-dqdq
            Call Find_RFO_Root(A_RFO_long,dqdq_long,
     &                         A_RFO_short,dqdq_short,
     &                         A_RFO,dqdq,StepMax)
            If (Iter.gt.IterMx) Then
               Write (Lu,*) ' Too many iterations in RF'
               Go To 997
            End If
            Go To 998
         End If
*
 997  Continue
      Call GetMem('Vector','Free','Real',ipVec,nInter*NumVal)
      Call GetMem('Values','Free','Real',ipVal,NumVal)
      If (nNeg.gt.0) Then
         mInter=nNeg+1
         Call GetMem('StepN','Free','Real',ipNStep,nInter)
         Call GetMem('GradN','Free','Real',ipNGrad,nInter)
         Call GetMem('VectorN','Free','Real',ipNVec,mInter)
         Call GetMem('ValuesN','Free','Real',ipNVal,1)
         Call GetMem('MatrixN','Free','Real',ipNMat,mInter*(mInter+1)/2)
         Call Free_Work(ipNTmp)
      End If
      mInter=nInter+1
      Call GetMem('StepP','Free','Real',ipPStep,nInter)
      Call GetMem('GradP','Free','Real',ipPGrad,nInter)
      Call GetMem('VectorP','Free','Real',ipPVec,mInter)
      Call GetMem('ValuesP','Free','Real',ipPVal,1)
      Call GetMem('MatrixP','Free','Real',ipPMat,mInter*(mInter+1)/2)
      Call Free_Work(ipPTmp)
*
      If (iPrint.ge.6) Then
         Write (Lu,*)
         Write (Lu,*)
         Write (Lu,*) 'Rational Function Optimization: Lambda=',Lambda
         Write (Lu,*)
      End If
      dqHdq=dqHdq+Lambda*Half
*
      If (iPrint.ge.99) Then
         Write (Lu,*) 'EigVal,dqHdq=',Lambda,dqHdq
         Call RecPrt(' In RS_P_RFO: g','(10f10.6)', g,nInter,1)
         Call RecPrt(' In RS_P_RFO:dq','(10f10.6)',dq,nInter,1)
      End If
*
*     Call QExit('RS_P_RFO')
      Return
      End
