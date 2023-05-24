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
      SUBROUTINE CHO_FINAL(WriteBookmarks)
      use ChoArr, only: iSOShl
C
C     Purpose: Cholesky finalizations.
C
      use ChoBkm, only: BkmVec, BkmThr, nRow_BkmVec, nCol_BkmVec,
     &                   nRow_BkmThr, nCol_BkmThr
      use stdalloc, only: mma_allocate, mma_deallocate
      use ChoIni
      Implicit None
      Logical WriteBookmarks
#include "cholesky.fh"
#include "choorb.fh"

      Integer, Allocatable:: BkmDim(:), iScratch(:)
      Real*8, Allocatable:: Scratch(:)
      INTEGER CHOISINI, IREO, l
      INTEGER NUMV(8)
#if defined (_DEBUGPRINT_)
      Integer is1CCD
#endif

C     Write NUMCHO array, shell indices, and threshold to runfile.
C     ------------------------------------------------------------

      CALL CHO_P_GETGV(NUMV,NSYM)
      CALL PUT_IARRAY('NUMCHO',NUMV,NSYM)
      CALL PUT_IARRAY('iSOShl',ISOSHL,NBAST)
      CALL PUT_DSCALAR('Cholesky Threshold',THRCOM)
#if defined (_DEBUGPRINT_)
      ! This is needed in order for bookmark tests in cho_x_init to work
      If (WriteBookmarks) Then
         If (Cho_1Center) Then
            is1CCD=1
         Else
            is1CCD=0
         End If
         Call Put_iScalar('1C-CD',is1CCD)
      End If
#endif

C     Write bookmarks to runfile.
C     First, transpose array.
C     ---------------------------

      If (WriteBookmarks) Then
         Call mma_allocate(BkmDim,4,Label='BkmDim')
         BkmDim(1)=nCol_BkmVec
         BkmDim(2)=nRow_BkmVec
         BkmDim(3)=nCol_BkmThr
         BkmDim(4)=nRow_BkmThr
         Call Put_iArray('Cholesky BkmDim',BkmDim,SIZE(BkmDim))
         Call mma_deallocate(BkmDim)
         If (nRow_BkmVec.gt.0 .and. nCol_BkmVec.gt.0 .and.
     &       nRow_BkmThr.gt.0 .and. nCol_BkmThr.gt.0) Then
            l=nRow_BkmVec*nCol_BkmVec
            Call mma_allocate(iScratch,l,Label='iScratch')
            Call iTrnsps(nSym,nCol_BkmVec,BkmVec,iScratch)
            Call Put_iArray('Cholesky BkmVec',iScratch,l)
            Call mma_deallocate(iScratch)
            Call mma_deallocate(BkmVec)
            nRow_BkmVec=0
            nCol_BkmVec=0
            l=nRow_BkmThr*nCol_BkmThr
            Call mma_allocate(Scratch,l,Label='Scratch')
            Call Trnsps(nSym,nCol_BkmThr,BkmThr,Scratch)
            Call Put_dArray('Cholesky BkmThr',Scratch,l)
            Call mma_deallocate(Scratch)
            Call mma_deallocate(BkmThr)
            nRow_BkmThr=0
            nCol_BkmThr=0
         End If
      End If
      If (Allocated(BkmVec)) Then
         Call mma_deallocate(BkmVec)
         nRow_BkmVec=0
         nCol_BkmVec=0
      End If
      If (Allocated(BkmThr)) Then
         Call mma_deallocate(BkmThr)
         nRow_BkmThr=0
         nCol_BkmThr=0
      End If

C     Write vector file address mode to runfile.
C     ------------------------------------------

      CALL PUT_ISCALAR('ChoVec Address',CHO_ADRVEC)

C     Write reorder mark to runfile.
C     ------------------------------

      IF (CHO_REORD) THEN
         IREO = 1
      ELSE
         IREO = 0
      END IF
      CALL PUT_ISCALAR('Cholesky Reorder',IREO)

C     Set initialization integer flag to "not set".
C     ---------------------------------------------

      CHOISINI = CHOINICHECK + 1
      CALL PUT_ISCALAR('ChoIni',CHOISINI)
      END
