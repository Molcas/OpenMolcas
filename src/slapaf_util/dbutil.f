************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
*                                                                      *
      Subroutine dBMult(dCdQ,QC,nQQ,nDim,nLambda)
*                                                                      *
************************************************************************
*                                                                      *
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "db.fh"
      Real*8 dCdQ(nQQ,nLambda), QC(nDim**2,nLambda)
*
      Call FZero(QC,nDim**2*nLambda)
      If (ip_dB.eq.ip_Dummy) Then
         Return
      End If
      Call Allocate_Work(ipX,mq*nLambda)
      Call FZero(Work(ipX),mq*nLambda)
      Call Allocate_Work(ipK,mq*nQQ)
      Call Get_dArray('K',Work(ipK),mq*nQQ)
*
      Call DGEMM_('N','N',mq,nLambda,nQQ,
     &            1.0D0,Work(ipK),mq,
     &                  dCdQ,nQQ,
     &            0.0D0,Work(ipX),mq)
      Call Free_Work(ipK)
*
      idB = 0
      Do iq = 0, mq-1
         nElem = iWork(ip_nqB + iq)
         Do iElem = idB, idB + (nElem**2)-1
            dBqR=Work(ip_dB+iElem)
            iDim=iWork(iElem*2 + ip_idB  )
            jDim=iWork(iElem*2 + ip_idB+1)
            ijDim = (jDim-1)*nDim + iDim
            Do iLambda = 1, nLambda
               X = Work((iLambda-1)*mq + iq + ipX)
               QC(ijDim,iLambda) = QC(ijDim,iLambda)
     &                           + X * dBqR
            End Do
         End Do
         idB = idB + nElem**2
      End Do
      Call Free_Work(ipX)
*
      Return
      End
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine dBPrint(nQQ,nDim)
*                                                                      *
************************************************************************
*                                                                      *
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "db.fh"
      If (ip_dB.eq.ip_Dummy) Return
      Call GetMem('dBQQ','Allo','Real',ipBuf,nDim**2)
      Call Allocate_Work(ipK,mq*nQQ)
      Call Get_dArray('K',Work(ipK),mq*nQQ)
      Do iQQ = 1, nQQ
         Call FZero(Work(ipBuf),nDim**2)
         idB = 0
         Do iq = 0, mq-1
            nElem = iWork(ip_nqB + iq)
            rK = Work((iQQ-1)*mq + iq + ipK)
            Do iElem = idB, idB + (nElem**2)-1
               dBqR=Work(ip_dB+iElem)
               iDim=iWork(iElem*2 + ip_idB  )
               jDim=iWork(iElem*2 + ip_idB+1)
               ijDim = (jDim-1)*nDim + iDim + ipBuf-1
               Work(ijDim) = Work(ijDim) + rK * dBqR
            End Do
            idB = idB + nElem**2
         End Do
         Call RecPrt('dBQQ',' ',Work(ipBuf),nDim,nDim)
      End Do
      Call Free_Work(ipBuf)
*
      Return
      End
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine dBuu(uM12,nQQ,nDim,g,Hss,Inv)
*                                                                      *
************************************************************************
*                                                                      *
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"
#include "real.fh"
#include "db.fh"
      Real*8 uM12(nDim), g(nQQ), Hss(nDim,nDim)
      Logical Inv
*
      If (ip_dB.eq.ip_Dummy) Then
         Call FZero(Hss,nDim**2)
         Return
      End If
      Call Allocate_Work(ipY,mq)
      Call Allocate_Work(ipK,mq*nQQ)
      Call Get_dArray('K',Work(ipK),mq*nQQ)
*
*     Compute Y(qR) = Sum_(Q) g(Q) rK(qR,Q)
*
      Call FZero(Work(ipY),mq)
      Do iQQ = 1, nQQ
         ipKk = (iQQ-1)*mq + ipK
         Call DaXpY_(mq,g(iQQ),Work(ipKk),1,Work(ipY),1)
      End Do
      Call Free_Work(ipK)
*
*     Compute Temp = Sum_(qR) Y(qR) * dB(qR)
*
      Call GetMem('temp','Allo','Real',ipTemp,nDim**2)
      Call FZero(Work(ipTemp),nDim**2)
      idB = 0
      Do iq = 0, mq-1
         YqR = Work(ipY+iq)
         nElem = iWork(ip_nqB + iq)
         Do iElem = idB, idB + (nElem**2)-1
            dBqR=Work(ip_dB+iElem)
            iDim=iWork(iElem*2 + ip_idB  )
            jDim=iWork(iElem*2 + ip_idB+1)
            ijDim = (jDim-1)*nDim + iDim + ipTemp-1
            Work(ijDim) = Work(ijDim) + YqR * dBqR
         End Do
         idB = idB + nElem**2
      End Do
      Call Free_Work(ipY)
*
      If (Inv) Call DScal_(nDim**2,-One,Work(ipTemp),1)
*
      Do i = 1, nDim
         Do j = 1, nDim
            xx = Sqrt(uM12(i)*uM12(j))
            ij = ipTemp + (j-1)*nDim + i - 1
            Hss(i,j) = Hss(i,j) +  Work(ij) / xx
         End Do
      End Do
*
      Call Free_Work(ipTemp)
*
      Return
      End
