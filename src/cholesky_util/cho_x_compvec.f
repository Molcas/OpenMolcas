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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine Cho_X_CompVec(irc,
     &                         NVT,l_NVT,
     &                         nBlock,l_nBlock,
     &                         nV,l_nV1,l_nV2,
     &                         iV1,l_iV11,l_iV12,
     &                         ip_Z,l_Z1,l_Z2,
     &                         Free_Z)
C
C     Thomas Bondo Pedersen, April 2010.
C
C     Purpose: Compute Cholesky vectors from Z vectors in core
C              and integrals computed on the fly:
C
C              L(uv,J) = 1/Z(J,J)
C                      * [ (uv|J) - sum[K=1,J-1] L(uv,K)*Z(J,K) ]
C
C     NVT(i): total number of Cholesky vectors (all nodes), sym. i
C     nBlock(i): number of vector blocks, sym. i
C     nV(i,j): number of vectors in block i of sym. j
C     iV1(i,j): index of first vector in block i of sym. j
C     ip_Z(i,j): pointer to Z block i of sym. j
C
C     If (Free_Z): Z array will be de-allocated here using ip_Z(1,1) as
C     start point of the array - i.e. assuming that the Z blocks are
C     stored as one coherent array. If this is not the case, simply put
C     Free_Z=.False. Letting this routine deallocate Z maximizes the
C     memory available for vector distribution.
C
C     Vectors are distributed across nodes and stored according to
C     reduced set 1 (all of them!).
C
#if defined (_DEBUGPRINT_)
      use ChoArr, only: iSP2F
#endif
      use ChoSwp, only: iQuAB, pTemp, iQuAB_here, nnBstRSh
      Implicit None
      Integer irc
      Integer l_NVT
      Integer l_nBlock
      Integer l_nV1, l_nV2
      Integer l_iV11, l_iV12
      Integer l_Z1, l_Z2
      Integer NVT(l_NVT)
      Integer nBlock(l_nBlock)
      Integer nV(l_NV1,l_NV2)
      Integer iV1(l_IV11,l_iV12)
      Integer ip_Z(l_Z1,l_Z2)
      Logical Free_Z
#include "cholesky.fh"
#include "choptr.fh"
#include "chosew.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "choprint.fh"

      Character*2  Unt
      Character*13 SecNam
      Parameter (SecNam='Cho_X_CompVec')

      Integer N2
      Parameter (N2=InfVec_N2)

      Integer  Cho_F2SP
      External Cho_F2SP

#if defined (_DEBUGPRINT_)
      Integer nBlock_Max, nnB
      Integer nTot, nTot2
      Integer myDebugInfo
      Parameter (myDebugInfo=100)
      Real*8  Tol
      Parameter (Tol=1.0d-14)
#endif

      Integer iAdr(8)
      Integer iSym, n
      Integer iSP_, iSP_1, iSP_2, iSp, nSP, nSP_Max, nSP_this_batch
      Integer ip_Int, l_Int
      Integer kOffI, kOffZ
      Integer jBlock, kBlock
      Integer J_inBlock, K_inBlock
      Integer kI, kL, kZ
      Integer ldL, ldZ
      Integer ip_Wrk, l_Wrk
      Integer MaxQual_SAVE
      Integer ip_InfVec_T
      Integer ip_Zd, l_Zd, incZd
      Integer kZd
      Integer ip_Tmp, l_Tmp
      Integer ip_ListSP, l_ListSP
      Integer iCountSP
      Integer lTot, Left
      Integer ip_BatchDim, l_BatchDim, nBatch

      Real*8 C0, C1, W0, W1
      Real*8 X0, X1, Y0, Y1
      Real*8 Byte, PMem, PDone
      Real*8 TotMem, TotCPU, TotWall

      Integer i, j, k
      Integer IndRSh, iTri, InfVcT
      IndRSh(i)=iWork(ip_IndRSh-1+i)
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j
      InfVcT(i,j,k)=iWork(ip_InfVec_T-1+MaxVec*N2*(k-1)+MaxVec*(j-1)+i)

      ! Init return code
      irc=0

#if defined (_DEBUGPRINT_)
      ! Check input
      If (l_NVT.lt.nSym .or. l_nBlock.lt.nSym .or.
     &    l_nV2.lt.nSym .or. l_iV12.lt.nSym .or.
     &    l_Z2.lt.nSym) Then
         irc=-2
         Return
      End If
      nBlock_Max=nBlock(1)
      Do iSym=2,nSym
         nBlock_Max=max(nBlock_Max,nBlock(iSym))
      End Do
      nnB=nBlock_Max*(nBlock_Max+1)/2
      If (l_nV1.lt.nBlock_Max .or. l_iV11.lt.nBlock_Max .or.
     &    l_Z1.lt.nnB) Then
         irc=-1
         Return
      End If
#endif

      ! Parallel runs: open tmp vector files
      Call Cho_XCV_TmpFiles(irc,1)
      If (irc .ne. 0) Then
         Write(LuPri,'(A,A,I8,A)')
     &   SecNam,': [1] Error in Cho_XCV_TmpFiles! (Return code:',irc,')'
         irc=1
         Return
      End If

      ! Get pointer to InfVec array for all vectors
      Call Cho_X_GetIP_InfVec(ip_InfVec_T)

      ! Copy reduced set 1 to location 2.
      ! I.e. make rs1 the "current" reduced set.
      Call Cho_X_RSCopy(irc,1,2)
      If (irc .ne. 0) Then
         Write(LuPri,'(A,A,I8,A)')
     &   SecNam,': Error in Cho_X_RSCopy! (Return code:',irc,')'
         irc=1
         Return
      End If

      ! Allocate tmp array to keep track of which shell pairs are
      ! computed on this node.
      l_Tmp=nnShl
      Call GetMem('XCVTMP','Allo','Inte',ip_Tmp,l_Tmp)

      ! Allocate and set list of parent shell pairs to compute
      ! NOTE: the tmp array allocated above is used here!!!
      Call iZero(iWork(ip_Tmp),nnShl)
      nSP=0
      Do iSym=1,nSym
         Do J=1,NVT(iSym)
            iSP=Cho_F2SP(IndRSh(InfVcT(J,1,iSym))) ! Parent SP
            If (iSP.gt.0 .and. iSP.le.nnShl) Then
               If (iWork(ip_Tmp-1+iSP) .eq. 0) Then
                  nSP=nSP+1
               End If
               iWork(ip_Tmp-1+iSP)=iWork(ip_Tmp-1+iSP)+1
            Else
               Call Cho_Quit(SecNam//': Parent SP error!!',103)
            End If
         End Do
      End Do
      l_ListSP=nSP
      Call GetMem('XCVLSP','Allo','Inte',ip_ListSP,l_ListSP)
      nSP=0
      Do iSP=0,nnShl-1
         If (iWork(ip_Tmp+iSP).gt.0) Then
            iWork(ip_ListSP+nSP)=iSP+1
            nSP=nSP+1
         End If
      End Do
#if defined (_DEBUGPRINT_)
      If (nSP.ne.l_ListSP) Then
         Call Cho_Quit('SP counting error [1] in '//SecNam,103)
      End If
      nTot=NVT(1)
      Do iSym=2,nSym
         nTot=nTot+NVT(iSym)
      End Do
      nTot2=0
      Do i=0,nSP-1
         nTot2=nTot2+iWork(ip_Tmp-1+iWork(ip_ListSP+i))
      End Do
      If (iPrint.ge.myDebugInfo) Then
         Call Cho_Head('Shell pair info from '//SecNam,'*',80,LuPri)
         Write(LuPri,'(/,A,I8,/)')
     &   'Number of shell pairs giving rise to vectors:',nSP
         Do i=0,nSP-1
            iSP=iWork(ip_ListSP+i)
            Call Cho_InvPck(iSP2F(iSP),j,k,.True.)
            Write(Lupri,'(A,I8,A,I4,A,I4,A,I8,A)')
     &      'Shell pair',iSP,' (shells',j,',',k,
     &      ') gives rise to',iWork(ip_Tmp-1+iSP),
     &      ' vectors'
         End Do
         Write(LuPri,'(A,I8,/,A,I8,/,A,I8)')
     &   'Total number of vectors (from SP).................',nTot2,
     &   'Total number of vectors (from NVT)................',nTot,
     &   'Difference........................................',nTot2-nTot
      End If
      If (nTot.ne.nTot2) Then
         Call Cho_Quit('SP counting error [2] in '//SecNam,103)
      End If
#endif

      ! Re-allocate and set qualification arrays
      pTemp => iQuAB
      MaxQual_SAVE=MaxQual
      MaxQual=NVT(1)
      Do iSym=2,nSym
         MaxQual=max(MaxQual,NVT(iSym))
      End Do
      Call mma_allocate(iQuAB_here,MaxQual,nSym,Label='iQuAB_here')
      iQuAB => iQuAB_here
      Do iSym=1,nSym
         nQual(iSym)=NVT(iSym)
         Do J=1,nQual(iSym)
            iQuAB(J,iSym)=InfVcT(J,1,iSym) ! parent product for vector J
         End Do
      End Do

      ! Modify diagonal of Z matrix
      ! Z(J,J) <- 1/Z(J,J)
      ! If Z memory should not be de-allocated, save a backup of the
      ! diagonal.
      If (Free_Z) Then
         l_Zd=1
         incZd=0
      Else
         l_Zd=NVT(1)
         Do iSym=2,nSym
            l_Zd=l_Zd+NVT(iSym)
         End Do
         incZd=1
      End If
      Call GetMem('XCVZd','Allo','Real',ip_Zd,l_Zd)
      kZd=ip_Zd
      Do iSym=1,nSym
         Do jBlock=1,nBlock(iSym)
            kOffZ=ip_Z(iTri(jBlock,jBlock),iSym)-1
            Do J_inBlock=1,nV(jBlock,iSym)
#if defined (_DEBUGPRINT_)
               ! Check for division by zero or negative diagonal
               ! This would be a bug....
               If (abs(Work(kOffZ+iTri(J_inBlock,J_inBlock))).lt.Tol
     &             .or. Work(kOffZ+iTri(J_inBlock,J_inBlock)).lt.-Tol)
     &         Then
                  Write(LuPri,'(//,A,A)') SecNam,': Ooooops!'
                  Write(LuPri,'(A)')
     &            '....division by small or negative number....'
                  Write(LuPri,'(A,I8)') 'iSym=',iSym
                  Write(LuPri,'(A,I8)') 'jBlock=',jBlock
                  Write(LuPri,'(A,I8)') 'J=',
     &            iV1(jBlock,iSym)+J_inBlock-1
                  Write(LuPri,'(A,1P,D15.6)')
     &            'Z(J,J)=',Work(kOffZ+iTri(J_inBlock,J_inBlock))
                  Write(LuPri,'(A,1P,D15.6)')
     &            'Tolerance=',Tol
                  Call Cho_Quit(
     &            'Division by small or negative number in '//SecNam,
     &            103)
               End If
#endif
               Work(kZd)=Work(kOffZ+iTri(J_inBlock,J_inBlock))
               Work(kOffZ+iTri(J_inBlock,J_inBlock))=1.0d0/Work(kZd)
               kZd=kZd+incZd
            End Do
         End Do
      End Do

      ! Distribute shell pairs across nodes according to dimension
      ! (helps load balance the linear algebra part below)
      Call Cho_XCV_Distrib_SP(iWork(ip_Tmp),l_Tmp,iCountSP)

      ! Allocate batch dimension array
      l_BatchDim=max(iCountSP,1)
      Call GetMem('XCVnBt','Allo','Inte',ip_BatchDim,l_BatchDim)

      ! Allocate offset array for batched SP loop
      l_iOff_Batch=nSym*nnShl
      Call GetMem('XCVOFB','Allo','Inte',ip_iOff_Batch,l_iOff_Batch)

      ! Split remaining memory in two parts.
      ! One for integrals/vectors, one for Seward.
      Call GetMem('XCVMx1','Max ','Real',ip_Wrk,l_Wrk)
      If (l_Wrk.lt.3) Then
         Call Cho_Quit('Insufficient memory in '//SecNam,101)
      End If
      l_Int=INT(DBLE(l_Wrk)*2.0d0/3.0d0)
      Call GetMem('XCVInt','Allo','Real',ip_Int,l_Int)
      Call GetMem('XCVMx2','Max ','Real',ip_Wrk,l_Wrk)
      Call xSetMem_Ints(l_Wrk)

      ! Print
      If (iPrint.ge.Inf_Pass) Then
         Call Cho_Head('Generation of Cholesky vectors','=',80,LuPri)
         Call Cho_Word2Byte(l_Int,8,Byte,Unt)
         Write(LuPri,'(/,A,I12,A,F10.3,1X,A,A)')
     &   'Memory available for integrals/vectors..',l_Int,' (',Byte,Unt,
     &   ')'
         Call Cho_Word2Byte(l_Wrk,8,Byte,Unt)
         Write(LuPri,'(A,I12,A,F10.3,1X,A,A)')
     &   'Memory available for Seward.............',l_Wrk,' (',Byte,Unt,
     &   ')'
         Write(LuPri,'(A,I12)')
     &   'Number of shell pairs, total (nnShl)....',nnShl
         Write(LuPri,'(A,I12)')
     &   'Number of shell pairs, this node........',iCountSP
         Write(LuPri,'(//,65X,A)') 'Time/min'
         Write(LuPri,'(1X,A,5X,A,5X,A,2X,A,10X,A,16X,A,6X,A)')
     &   'Batch','iSP1','iSP2','%Done','Memory','CPU','Wall'
         Write(LuPri,'(A,A)')
     &   '------------------------------------------------------------',
     &   '----------------'
         Call Cho_Flush(LuPri)
         TotMem=0.0d0
         TotCPU=0.0d0
         TotWall=0.0d0
      End If

      ! Compute Cholesky vectors in batched loop over shell pairs
      Call iZero(iAdr,nSym) ! disk addresses
      iSP_1=1
      nBatch=0
      Do While (iSP_1.le.iCountSP)
         If (iPrint.ge.Inf_Pass) Then
            Call Cho_Timer(X0,Y0)
         End If
         ! Set batch info
         Left=l_Int
         nSP_Max=iCountSP-iSP_1+1
         nSP_this_batch=0
         Do While (Left.gt.0 .and. nSP_this_batch.lt.nSP_Max)
            iSP=iSP_1+nSP_this_batch
            n=nnBstRSh(1,iWork(ip_Tmp-1+iSP),2)*NVT(1)
            Do iSym=2,nSym
               n=n+nnBstRSh(iSym,iWork(ip_Tmp-1+iSP),2)*NVT(iSym)
            End Do
            If (n.le.Left) Then
               Left=Left-n
               nSP_this_batch=nSP_this_batch+1
            Else
               Left=0 ! break while loop
            End If
         End Do
         If (nSP_this_batch.lt.1) Then
            Call Cho_Quit(
     &         SecNam//': Insufficient memory for shell pair batch',101)
         End If
         iSP_2=iSP_1+nSP_this_batch-1
         Do iSym=1,nSym
            nDim_Batch(iSym)=0
            Do iSP_=iSP_1,iSP_2
               iSP=iWork(ip_Tmp-1+iSP_)
               iWork(ip_iOff_Batch-1+nSym*(iSP-1)+iSym)=nDim_Batch(iSym)
               nDim_Batch(iSym)=nDim_Batch(iSym)+nnBstRSh(iSym,iSP,2)
            End Do
         End Do
         ! Calculate integrals (uv|J) for uv in shell pair batch and
         ! for all J (shell pairs that give rise to vectors are
         ! listed in ListSP).
         Call Cho_XCV_GetInt(irc,iWork(ip_Tmp-1+iSP_1),nSP_this_batch,
     &                       iWork(ip_ListSP),l_ListSP,
     &                       NVT,l_NVT,Work(ip_Int),l_Int)
         If (irc .ne. 0) Then
            Write(LuPri,'(A,A,I8)')
     &      SecNam,': Cho_XCV_GetInt returned code',irc
            Call Cho_Quit(SecNam//': Error in Cho_XCV_GetInt',104)
         End If
         ! Convert integrals into Cholesky vectors in each symmetry
         Call Cho_Timer(C0,W0)
         kOffI=ip_Int
         Do iSym=1,nSym
            ldL=max(nDim_Batch(iSym),1)
            Do jBlock=1,nBlock(iSym)
               kOffZ=ip_Z(iTri(jBlock,jBlock),iSym)-1
               Do J_inBlock=1,nV(jBlock,iSym)
                  ! Convert integral column into Cholesky vector
                  ! L(uv,J)=(uv|J)/Z(J,J)
                  ! Note that the inverse was taken before the integral
                  ! loop - i.e. Z(J,J) <- 1/Z(J,J)
                  J=iV1(jBlock,iSym)+J_inBlock-1
                  kL=kOffI+nDim_Batch(iSym)*(J-1)
                  Call dScal_(nDim_Batch(iSym),
     &                       Work(kOffZ+iTri(J_inBlock,J_inBlock)),
     &                       Work(kL),1)
                  ! Subtract from subsequent columns in current block
                  ! using BLAS1
                  ! (uv|K) <- (uv|K) - L(uv,J)*Z(K,J), K>J (in jBlock)
                  Do K_inBlock=J_inBlock+1,nV(jBlock,iSym)
                     K=iV1(jBlock,iSym)+K_inBlock-1
                     Call dAXPY_(nDim_Batch(iSym),
     &                    -Work(kOffZ+iTri(K_inBlock,J_inBlock)),
     &                    Work(kL),1,
     &                    Work(kOffI+nDim_Batch(iSym)*(K-1)),1)
                  End Do
               End Do
               ! Subtract from subsequent blocks using BLAS3
               ! (uv|K) <- (uv|K) - sum[J] L(uv,J)*Z(K,J),
               ! K>J (K in kBlock>jBlock containing J)
               kL=kOffI+nDim_Batch(iSym)*(iV1(jBlock,iSym)-1)
               Do kBlock=jBlock+1,nBlock(iSym)
                  kI=kOffI+nDim_Batch(iSym)*(iV1(kBlock,iSym)-1)
                  kZ=ip_Z(iTri(kBlock,jBlock),iSym)
                  ldZ=max(nV(kBlock,iSym),1)
                  Call dGeMM_('N','T',
     &                        nDim_Batch(iSym),nV(kBlock,iSym),
     &                        nV(jBlock,iSym),
     &                        -1.0d0,Work(kL),ldL,Work(kZ),ldZ,
     &                        1.0d0,Work(kI),ldL)
               End Do
            End Do
            nDGM_Call=nDGM_Call+nBlock(iSym)*(nBlock(iSym)-1)/2
            kOffI=kOffI+nDim_Batch(iSym)*NVT(iSym)
         End Do
         Call Cho_Timer(C1,W1)
         tDecom(1,3)=tDecom(1,3)+(C1-C0)
         tDecom(2,3)=tDecom(2,3)+(W1-W0)
         ! Write vectors to temp files
         Call Cho_Timer(C0,W0)
         kL=ip_Int
         Do iSym=1,nSym
            lTot=nDim_Batch(iSym)*NVT(iSym)
            If (lTot.gt.0) Then
               Call DDAFile(LuTmp(iSym),1,Work(kL),lTot,iAdr(iSym))
               kL=kL+lTot
            End If
         End Do
         Call Cho_Timer(C1,W1)
         tDecom(1,2)=tDecom(1,2)+(C1-C0)
         tDecom(2,2)=tDecom(2,2)+(W1-W0)
         ! Print
         If (iPrint.ge.Inf_Pass) Then
            Call Cho_Timer(X1,Y1)
            PDone=1.0d2*DBLE(iSP_2)/DBLE(iCountSP)
            lTot=nDim_Batch(1)*NVT(1)
            Do iSym=2,nSym
               lTot=lTot+nDim_Batch(iSym)*NVT(iSym)
            End Do
            Call Cho_Word2Byte(lTot,8,Byte,Unt)
            PMem=1.0d2*DBLE(lTot)/DBLE(l_Int)
            Write(LuPri,
     &'(I6,1X,I8,1X,I8,1X,F6.1,1X,F10.3,1X,A,A,F7.2,A,1X,F9.2,1X,F9.2)')
     &      nBatch+1,iSP_1,iSP_2,PDone,Byte,Unt,' (',PMem,'%)',
     &      (X1-X0)/6.0d1,(Y1-Y0)/6.0d1
            Call Cho_Flush(LuPri)
            TotMem=TotMem+DBLE(lTot)
            TotCPU=TotCPU+(X1-X0)/6.0d1
            TotWall=TotWall+(Y1-Y0)/6.0d1
         End If
         ! Update counters and save batch dimension
         iSP_1=iSP_1+nSP_this_batch
         iWork(ip_BatchDim+nBatch)=nSP_this_batch
         nBatch=nBatch+1
      End Do
      If (iPrint.ge.Inf_Pass) Then
         Write(LuPri,'(A,A)')
     &   '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ',
     &   '- - - - - - - - '
         Call Cho_RWord2Byte(TotMem,Byte,Unt)
         Write(LuPri,'(32X,F10.3,1X,A,12X,F9.2,1X,F9.2)')
     &   Byte,Unt,TotCPU,TotWall
         Write(LuPri,'(A,A)')
     &   '------------------------------------------------------------',
     &   '----------------'
         Call Cho_Flush(LuPri)
      End If

      ! De-allocate Z matrix or reconstruct diagonal of Z matrix
      ! Z(J,J) <- 1/Z(J,J)
      If (Free_Z) Then
         lTot=NVT(1)*(NVT(1)+1)/2
         Do iSym=2,nSym
            lTot=lTot+NVT(iSym)*(NVT(iSym)+1)/2
         End Do
         Call GetMem('XCVFZ','Free','Real',ip_Z(1,1),lTot)
      Else
         kZd=ip_Zd-1
         Do iSym=1,nSym
            Do jBlock=1,nBlock(iSym)
               kOffZ=ip_Z(iTri(jBlock,jBlock),iSym)-1
               Do J_inBlock=1,nV(jBlock,iSym)
                  J=iV1(jBlock,iSym)+J_inBlock-1
                  Work(kOffZ+iTri(J_inBlock,J_inBlock))=Work(kZd+J)
               End Do
            End Do
            kZd=kZd+NVT(iSym)
         End Do
      End If

      ! Deallocations
      Call xRlsMem_Ints()
      Call GetMem('XCVInt','Free','Real',ip_Int,l_Int)
      Call GetMem('XCVOFB','Free','Inte',ip_iOff_Batch,l_iOff_Batch)
      l_iOff_Batch=0
      Call GetMem('XCVZd','Free','Real',ip_Zd,l_Zd)
      iQuAB => Null()
      Call mma_deallocate(iQuAB_here)
      Call GetMem('XCVLSP','Free','Inte',ip_ListSP,l_ListSP)

      ! Reset qualification array pointers and MaxQual
      iQuAB => pTemp
      MaxQual=MaxQual_SAVE

      ! Parallel runs: distribute vectors across nodes (store on files)
      ! Serial runs: write vectors to permanent files
      Call Cho_Timer(C0,W0)
      Call Cho_XCV_DistributeVectors(irc,iWork(ip_BatchDim),nBatch,
     &                               iWork(ip_Tmp),iCountSP,NVT,l_NVT)
      If (irc .ne. 0) Then
         Write(LuPri,'(A,A,I8)')
     &   SecNam,': Cho_XCV_DistributeVectors returned code',irc
         Call Cho_Quit(SecNam//': Error in Cho_XCV_DistributeVectors',
     &                 104)
      End If
      Call Cho_Timer(C1,W1)
      tDecom(1,2)=tDecom(1,2)+(C1-C0)
      tDecom(2,2)=tDecom(2,2)+(W1-W0)

      ! Deallocations
      Call GetMem('XCVnBt','Free','Inte',ip_BatchDim,l_BatchDim)
      Call GetMem('XCVTMP','Free','Inte',ip_Tmp,l_Tmp)

      ! Parallel runs: close and erase tmp vector files
      Call Cho_XCV_TmpFiles(irc,3)
      If (irc .ne. 0) Then
         Write(LuPri,'(A,A,I8,A)')
     &   SecNam,': [3] Error in Cho_XCV_TmpFiles! (Return code:',irc,')'
         Call Cho_Quit(SecNam//': Error in Cho_XCV_TmpFiles',104)
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C TMP FILE HANDLERS (OPEN, CLOSE, DELETE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SubRoutine Cho_XCV_TmpFiles(irc,iOpt)
      Implicit None
      Integer irc, iOpt
      irc=0
      If (iOpt.eq.1) Then
         Call Cho_XCV_OpenTmpFiles()
      Else If (iOpt.eq.2) Then
         Call Cho_XCV_CloseAndKeepTmpFiles()
      Else If (iOpt.eq.3) Then
         Call Cho_XCV_CloseAndEraseTmpFiles()
      Else
         irc=1
         Return
      End If
      End
      SubRoutine Cho_XCV_OpenTmpFiles()
      Implicit None
#include "cholesky.fh"
      Integer iSym
      Character*6 Filename
      Do iSym=1,nSym
         LuTmp(iSym)=7
         Write(Filename,'(A4,I2.2)') 'VTMP',iSym
         Call DAName_MF_WA(LuTmp(iSym),Filename)
      End Do
      End
      SubRoutine Cho_XCV_CloseAndKeepTmpFiles()
      Implicit None
#include "cholesky.fh"
      Integer iSym
      Do iSym=1,nSym
         If (LuTmp(iSym).gt.0) Then
            Call DAClos(LuTmp(iSym))
            LuTmp(iSym)=0
         End If
      End Do
      End
      SubRoutine Cho_XCV_CloseAndEraseTmpFiles()
      Implicit None
#include "cholesky.fh"
      Integer iSym
      Do iSym=1,nSym
         If (LuTmp(iSym).gt.0) Then
            Call DAEras(LuTmp(iSym))
            LuTmp(iSym)=0
         End If
      End Do
      End
