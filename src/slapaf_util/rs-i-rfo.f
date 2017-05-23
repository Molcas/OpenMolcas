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
      Subroutine RS_I_RFO(H,g,nInter,dq,UpMeth,dqHdq,StepMax,Step_Trunc)
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
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 H(nInter,nInter), g(nInter), dq(nInter)
*
      Character*6 UpMeth
      Character*1 Step_Trunc
      Logical Found
*
*     Call QEnter('RS_I_RFO')
      iRout = 215
      Lu =6
      iPrint = nPrint(iRout)
      If (iPrint.ge.99) Then
         Call RecPrt(' In RS_I_RFO: H','(10f10.6)',H,nInter,nInter)
         Call RecPrt(' In RS_I_RFO: g','(10f10.6)', g,nInter,1)
         Call RecPrt(' In RS_I_RFO:dq','(10f10.6)',dq,nInter,1)
      End If
*
      NumVal=2
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
*
*---- Find the negative eigenvalue(s)
*     Stop when the highest eigenvalue found is larger than Thr
      Do While (.Not.Found)
        Call Davidson(Work(ipMat),nInter,NumVal,
     &                Work(ipVal),Work(ipVec),iStatus)
        If (iStatus.gt.0) Then
          Call SysWarnMsg('RS_I_RFO',
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
          NumVal=NumVal+nVStep
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
         Call RecPrt(' In RS_I_RFO: Eigenvalues',' ',Work(ipVal),
     &               1,NumVal)
         Call RecPrt(' In RS_I_RFO: Eigenvectors',' ',Work(ipVec),
     &               nInter,NumVal)
         Write (Lu,*) ' nNeg=',nNeg
      End If
*
*     Transform the gradient and Hessian to generate the
*     corresponding entities for the image function. This
*     corresponds to an elementary Householder orthogonal
*     transformation.
*
      Call Allocate_Work(ipTmp,nInter)
      call dcopy_(nInter,g,1,Work(ipTmp),1)
*
      Do iNeg=1,nNeg
        iVec = ipVec+(iNeg-1)*nInter
        gi = DDot_(nInter,g,1,Work(iVec),1)
        Call DaXpY_(nInter,-Two*gi,Work(iVec),1,g,1)
        Fact = Two * Work(ipVal+iNeg-1)
        Do j = 1, nInter
           Do i = 1, nInter
              H(i,j) = H(i,j) - Fact * Work(iVec+i-1) * Work(iVec+j-1)
           End Do
        End Do
      End Do
*
      Call GetMem('Vector','Free','Real',ipVec,nInter*NumVal)
      Call GetMem('Values','Free','Real',ipVal,NumVal)
*
      Call RS_RFO(H,g,nInter,dq,UpMeth,dqHdq,StepMax,Step_Trunc)
*
*     Restore the original gradient
*
      call dcopy_(nInter,Work(ipTmp),1,g,1)
      Call Free_Work(ipTmp)
*
      UpMeth='RSIRFO'
*
      If (iPrint.ge.99) Then
         Call RecPrt(' In RS_I_RFO: g','(10f10.6)', g,nInter,1)
         Call RecPrt(' In RS_I_RFO:dq','(10f10.6)',dq,nInter,1)
      End If
*
*     Call QExit('RS_I_RFO')
      Return
      End
