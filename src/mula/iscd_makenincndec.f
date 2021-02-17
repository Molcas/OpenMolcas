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
      Subroutine ISCD_MakenIncDec(n_max,nOrd,nOsc,lNMAT,lNINC,lNDEC,
     &                            lBatch,nBatch,leftBatch,nIndex,
     &                            Graph2, nMat,nInc,nDec)
C!
      Implicit Real*8 ( a-h,o-z )
      Integer nMat(nOsc,lBatch), nInc(nOsc,lBatch), nDec(nOsc,lBatch)
      Integer Graph2(n_max+1,n_max+1,nOsc)
      Integer n_max, nOrd, lBatch,nBatch,leftBatch
      Integer lNMAT, lNINC, lNDEC
#include "WrkSpc.fh"
#include "io_mula.fh"
      Integer nIndex(3,0:maxMax_n)

CGGt -------------------------------------------------------------------
c      Write(6,*)
c      Write(6,*)'CGGt[ISCD_Mk_nIncDec] Infos:                   '
c      Write(6,*)'     nMat(',nOsc,',',lBatch,')'
c      Write(6,*)'     n_max,nOrd,nOsc==',n_max,nOrd,nOsc
c      Write(6,*)'     lBatch,nBatch,leftBatch==',lBatch,nBatch,leftBatch
c      Write(6,*)'----------------------------------------------'
c      Write(6,*)'  The nIndex file:    '                          ! CGGt
c      Do i=1,nBatch+1                                             ! CGGt
c      Write(6,*) i,': ',nIndex(1,i)                               ! CGGt
c      EndDo                                                       ! CGGt
c      Write(6,*)'----------------------------------------------'
c      Call XFlush(6)                                              ! CGGt
CGGt -------------------------------------------------------------------

C!
C!---- Initialize
C!
      Call GetMem('iVecI','Allo','INTE',ipiVecI,nOsc)
      Call GetMem('iVecD','Allo','INTE',ipiVecD,nOsc)
C!
C!---- Macrocycle iBatch
C!     Reading nMat
C!
      iIndex = 0
      jIndex = 0
      DO iBatch = 1, nBatch
        do ii=1,lBatch
          do iv=1,nOsc
            nInc(iv,ii)=-1
            nDec(iv,ii)=-1
          enddo
        enddo
        kIndex = nIndex(1,iBatch)
c      Write(6,*)'          iBatch=',iBatch,'  kIndex=',kIndex     ! CGGt
        Call iDaFile(lNMAT,2,nMat,nOsc*lBatch,kIndex)
CGGt -------------------------------------------------------------------
c      Do i=0,lBatch-1                                             ! CGGt
c      Write(6,*) i+(iBatch-1)*lBatch,': ',(nMat(k,i+1),k=1,nOsc)  ! CGGt
c      EndDo                                                       ! CGGt
c      Write(6,*)'----------------------------------------------'  ! CGGt
c      Call XFlush(6)                                              ! CGGt
CGGt -------------------------------------------------------------------
        Do ii = 1, lBatch
C!
C!---- Create nInc.
C!
          do iv=1, nOsc
            iWork(ipiVecI+iv-1) = nMat(iv,ii)
          enddo
          Do j = 1,nOsc
            iWork(ipiVecI+j-1) = iWork(ipiVecI+j-1)+1
            nInc(j,ii)=iDetnr(iWork(ipiVecI),Graph2,nOsc,n_max)
            iWork(ipiVecI+j-1) = iWork(ipiVecI+j-1)-1
          EndDo
C!
C!---- Create nDec.
C!
          Do j = 1,nOsc
            If (nMat(j,ii).ne.0)Then
              do iv=1,nOsc
                iWork(ipiVecD+iv-1) = nMat(iv,ii)
              enddo
              iWork(ipiVecD+j-1)=iWork(ipiVecD+j-1)-1
              nDec(j,ii)=iDetnr(iWork(ipiVecD),Graph2,nosc,n_max)
              do iv=1,nOsc
                iWork(ipiVecD+iv-1) = iWork(ipiVecD+j-1)+1
              enddo
            Else
              nDec(j,ii)=-1
            End IF
          End Do
C!
        EndDo
        nIndex(2,iBatch) = iIndex
        Call iDaFile(lNINC,1,nInc,nOsc*lBatch,iIndex)
        nIndex(3,iBatch) = jIndex
        Call iDaFile(lNDEC,1,nDec,nOsc*lBatch,jIndex)
CGGt -------------------------------------------------------------------
c      Do i=0,lBatch-1                                             ! CGGt
c      Write(6,*) i+(iBatch-1)*lBatch,':I',(nInc(k,i+1),k=1,nOsc)  ! CGGt
c      Write(6,*) i+(iBatch-1)*lBatch,':D',(nDec(k,i+1),k=1,nOsc)  ! CGGt
c      EndDo                                                       ! CGGt
c      Call XFlush(6)                                              ! CGGt
CGGt -------------------------------------------------------------------
c      Write(6,*)'            nInc Written at ',nIndex(2,iBatch)   ! CGGt
c      Write(6,*)'            nDec Written at ',nIndex(3,iBatch)   ! CGGt
c      Write(6,*)'----------------------------------------------'  ! CGGt
c      Call XFlush(6)                                              ! CGGt
CGGt -------------------------------------------------------------------
C!
      END DO
C!
      IF (leftBatch.GT.0) then
        kIndex = nIndex(1,nBatch+1)
c      Write(6,*)'          nBatch+1',nBatch+1,'  kIndex=',kIndex  ! CGGt
        Call iDaFile(lNMAT,2,nMat,nOsc*lBatch,kIndex)
CGGt -------------------------------------------------------------------
c      Write(6,*)'  --------- last Batch'                          ! CGGt
c      Do i=0,lBatch-1                                             ! CGGt
c      Write(6,*) i+(iBatch-1)*lBatch,': ',(nMat(k,i+1),k=1,nOsc)  ! CGGt
c      EndDo                                                       ! CGGt
c      Call XFlush(6)                                              ! CGGt
CGGt -------------------------------------------------------------------
C!
c      Write(6,*)'CGGt nBatch*lBatch,nOrd==',nBatch*lBatch,nOrd    ! CGGt
c      Call XFlush(6)                                              ! CGGt
        do ii=1,lBatch
          do iv=1,nOsc
            nInc(iv,ii)=-1
            nDec(iv,ii)=-1
          enddo
        enddo
        Do iOrd = nBatch*lBatch , nOrd
          ii = 1+ iOrd - nBatch*lBatch
C!
C!---- Create nInc.
C!
          do iv=1, nOsc
            iWork(ipiVecI+iv-1) = nMat(iv,ii)
          enddo
          Do j = 1,nOsc
            iWork(ipiVecI+j-1) = iWork(ipiVecI+j-1)+1
            nInc(j,ii)=iDetnr(iWork(ipiVecI),Graph2,nOsc,n_max)
            iWork(ipiVecI+j-1) = iWork(ipiVecI+j-1)-1
          EndDo
C!
C!---- Create nDec.
C!
          Do j = 1,nOsc
            If (nMat(j,ii).ne.0)Then
              do iv=1,nOsc
                iWork(ipiVecD+iv-1) = nMat(iv,ii)
              enddo
              iWork(ipiVecD+j-1)=iWork(ipiVecD+j-1)-1
              nDec(j,ii)=iDetnr(iWork(ipiVecD),Graph2,nosc,n_max)
              do iv=1,nOsc
                iWork(ipiVecD+iv-1) = iWork(ipiVecD+j-1)+1
              enddo
            Else
              nDec(j,ii)=-1
            End IF
          End Do
C!
        EndDo
C!
        nIndex(2,nBatch+1) = iIndex
        Call iDaFile(lNINC,1,nInc,nOsc*lBatch,iIndex)
        nIndex(3,nBatch+1) = jIndex
        Call iDaFile(lNDEC,1,nDec,nOsc*lBatch,jIndex)
CGGt -------------------------------------------------------------------
c      Do i=0,lBatch-1                                             ! CGGt
c      Write(6,*) i+(iBatch-1)*lBatch,':I',(nInc(k,i+1),k=1,nOsc)  ! CGGt
c      Write(6,*) i+(iBatch-1)*lBatch,':D',(nDec(k,i+1),k=1,nOsc)  ! CGGt
c      EndDo                                                       ! CGGt
c      Call XFlush(6)                                              ! CGGt
CGGt -------------------------------------------------------------------
c      Write(6,*)'            nInc Written at ',nIndex(2,nBatch+1) ! CGGt
c      Write(6,*)'            nDec Written at ',nIndex(3,nBatch+1) ! CGGt
c      Write(6,*)'----------------------------------------------'  ! CGGt
c      Call XFlush(6)                                              ! CGGt
CGGt -------------------------------------------------------------------
      ENDIF
C!
c      Write(6,*)'----------------------------------------------'  ! CGGt
c      Call XFlush(6)                                              ! CGGt
C!
      Call GetMem('iVecD','Free','INTE',ipiVecD,nOsc)
      Call GetMem('iVecI','Free','INTE',ipiVecI,nOsc)
      Return
      End


      Subroutine ISCD_ReloadNMAT(lnTabDim,nOrd,nOsc,lNMAT0,lNMAT,
     &               lBatch,nBatch,leftBatch, nIndex,nTabDim,nMat0,nMat)
C!
      Implicit Real*8 ( a-h,o-z )
      Integer nMat(nOsc,lBatch)
      Integer nMat0(nOsc), nTabDim(0:lnTabDim)
      Integer lnTabDim,nOrd,nOsc,lNMAT0,lNMAT
      Integer lBatch,nBatch,leftBatch
#include "WrkSpc.fh"
#include "io_mula.fh"
      Integer nIndex(3,0:maxMax_n)
C!
C!---- Initialize
C!
CGGt -------------------------------------------------------------------
c      Write(6,*)
c      Write(6,*)'CGGt[ISCD_ReloadNMAT] Infos:                   '
c      Write(6,*)'     nMat(',nOsc,',',lBatch,')'
c      Write(6,*)'     lnTabDim,nOrd,nOsc==',lnTabDim,nOrd,nOsc
c      Write(6,*)'     lBatch,nBatch,leftBatch==',lBatch,nBatch,leftBatch
c      Write(6,*)'     lnTabDim+1=',lnTabDim+1,':'
c      Do i=0, lnTabDim
c      Write(6,*) i,' raed at ',nTabDim(i)
c      iIndex0 = nTabDim(i)                                        ! CGGt
c      Call iDaFile(lNMAT0,2,nMat0,nOsc,iIndex0)
c      Write(6,*) i,' read at',nTabDim(i),'  M:',(nMat0(j),j=1,nOsc)
c      EndDo
c      Write(6,*)'----------------------------------------------'
c      Call XFlush(6)
CGGt -------------------------------------------------------------------
      jIndex = 0
      Rewind(lNMAT0)
      Do iBatch = 1, nBatch
        Do ii = 0, lBatch-1
          iOrd = ii + (iBatch-1)*lBatch
          iIndex0 = nTabDim(iOrd)
          Call iDaFile(lNMAT0,2,nMat0,nOsc,iIndex0)
          Do iOsc = 1, nOsc
            nMat(iOsc,ii+1) = nMat0(iOsc)
          EndDo
        EndDo
        nIndex(1,iBatch) = jIndex
        Call iDaFile(lNMAT,1,nMat,nOsc*lBatch,jIndex)
CGGt -------------------------------------------------------------------
c      Write(6,*)'  --------- iBatch =',iBatch                     ! CGGt
c      Do i=0,lBatch-1                                             ! CGGt
c      Write(6,*) i+(iBatch-1)*lBatch,':M',(nMat(k,i+1),k=1,nOsc)  ! CGGt
c      EndDo                                                       ! CGGt
c      Call XFlush(6)                                              ! CGGt
CGGt -------------------------------------------------------------------
      EndDo
      If (leftBatch.GT.0) then
        Do iOrd = nBatch*lBatch , nOrd
          ii = iOrd - nBatch*lBatch
          iIndex0 = nTabDim(iOrd)
          Call iDaFile(lNMAT0,2,nMat0,nOsc,iIndex0)
          Do iOsc = 1, nOsc
            nMat(iOsc,ii+1) = nMat0(iOsc)
          EndDo
        EndDo
        nIndex(1,nBatch+1) = jIndex
        Call iDaFile(lNMAT,1,nMat,nOsc*lBatch,jIndex)
CGGt -------------------------------------------------------------------
c      Write(6,*)'  --------- last Batch'                          ! CGGt
c      Do i=0,lBatch-1                                             ! CGGt
c      Write(6,*) i+(iBatch-1)*lBatch,':M',(nMat(k,i+1),k=1,nOsc)  ! CGGt
c      EndDo                                                       ! CGGt
c      Call XFlush(6)                                              ! CGGt
CGGt -------------------------------------------------------------------
      EndIf
CGGt -------------------------------------------------------------------
c      Write(6,*)'----------------------------------------------'
c      Write(6,*)'  The nIndex file:    '                          ! CGGt
c      Do i=1,nBatch+1                                             ! CGGt
c      Write(6,*) i,': ',nIndex(1,i)                               ! CGGt
c      EndDo                                                       ! CGGt
c      Write(6,*)'----------------------------------------------'
c      Call XFlush(6)                                              ! CGGt
CGGt -------------------------------------------------------------------
C!
      Return
      End
