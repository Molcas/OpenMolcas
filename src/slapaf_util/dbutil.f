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
#include "stdalloc.fh"
#include "db.fh"
      Real*8 dCdQ(nQQ,nLambda), QC(nDim**2,nLambda)
      Real*8, Allocatable:: X(:,:), K(:,:)
*
      QC(:,:)=0.0D0

      If (ip_dB.eq.ip_Dummy) Then
C        Write (6,*) 'FAST out'
         Return
      End If

      Call mma_allocate(X,mq,nLambda,Label='X')
      X(:,:)=0.0D0
      Call mma_allocate(K,mq,nQQ,Label='K')
      Call Get_dArray('K',K,mq*nQQ)
*
      Call DGEMM_('N','N',mq,nLambda,nQQ,
     &            1.0D0,K,mq,
     &                  dCdQ,nQQ,
     &            0.0D0,X,mq)
      Call mma_deallocate(K)
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
               QC(ijDim,iLambda) = QC(ijDim,iLambda)
     &                           + X(1+iq,iLambda) * dBqR
            End Do
         End Do
         idB = idB + nElem**2
      End Do
      Call mma_deallocate(X)
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
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "db.fh"
      Real*8, Allocatable:: dBQQ(:,:), K(:,:)

      If (ip_dB.eq.ip_Dummy) Return
      Call mma_allocate(dBQQ,nDim, nDim,Label='dBQQ')
      Call mma_allocate(K,mq,nQQ,Label='K')
      Call Get_dArray('K',K,mq*nQQ)
      Do iQQ = 1, nQQ
         dBQQ(:,:) = Zero
         idB = 0
         Do iq = 0, mq-1
            nElem = iWork(ip_nqB + iq)
            rK = K(iq+1,iQQ)
            Do iElem = idB, idB + (nElem**2)-1
               dBqR=Work(ip_dB+iElem)
               iDim=iWork(iElem*2 + ip_idB  )
               jDim=iWork(iElem*2 + ip_idB+1)
               dBQQ(iDim,jDim) = dBQQ(iDim,jDim) + rK * dBqR
            End Do
            idB = idB + nElem**2
         End Do
         Call RecPrt('dBQQ',' ',dBQQ,nDim,nDim)
      End Do
      Call mma_deallocate(dBQQ)
      Call mma_deallocate(K)
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
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "db.fh"
      Real*8 uM12(nDim), g(nQQ), Hss(nDim,nDim)
      Logical Inv
      Real*8, Allocatable:: Y(:), K(:,:), Temp(:,:)
*
      If (ip_dB.eq.ip_Dummy) Then
         Hss(:,:)=Zero
         Return
      End If

      Call mma_allocate(Y,mq,Label='Y')
      Call mma_allocate(K,mq,nQQ,Label='K')
      Call Get_dArray('K',K,mq*nQQ)
*
*     Compute Y(qR) = Sum_(Q) g(Q) rK(qR,Q)
*
      Y(:) = Zero
      Do iQQ = 1, nQQ
         Call DaXpY_(mq,g(iQQ),K(:,iQQ),1,Y,1)
      End Do
      Call mma_deallocate(K)
*
*     Compute Temp = Sum_(qR) Y(qR) * dB(qR)
*
      Call mma_allocate(Temp,nDim,nDim,Label='Temp')
      Temp(:,:)=Zero

      idB = 0
      Do iq = 0, mq-1
         YqR = Y(1+iq)
         nElem = iWork(ip_nqB + iq)
         Do iElem = idB, idB + (nElem**2)-1
            dBqR=Work(ip_dB+iElem)
            iDim=iWork(iElem*2 + ip_idB  )
            jDim=iWork(iElem*2 + ip_idB+1)
            Temp(iDim,jDim) = Temp(iDim,jDim) + YqR * dBqR
         End Do
         idB = idB + nElem**2
      End Do
      Call mma_deallocate(Y)
*
      If (Inv) Call DScal_(nDim**2,-One,Temp,1)
*
      Do i = 1, nDim
         Do j = 1, nDim
            xx = Sqrt(uM12(i)*uM12(j))
            Hss(i,j) = Hss(i,j) + Temp(i,j) / xx
         End Do
      End Do
*
      Call mma_deallocate(Temp)
*
      Return
      End
