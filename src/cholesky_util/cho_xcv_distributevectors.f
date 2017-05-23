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
*               2012,2014, Victor P. Vysotskiy                         *
************************************************************************
      SubRoutine Cho_XCV_DistributeVectors(irc,SP_BatchDim,nSP_Batch,
     &                                     idSP,n_idSP,NVT,l_NVT)
C
C     Thomas Bondo Pedersen, April 2010.
C
C     Parallel execution: distribute vectors across nodes.
C     Serial execution: reorder vectors on tmp files and write them to
C     permanent vector files.
C
C     Victor P. Vysotskiy, 2012:
C     Number of 'ga_put' has been remarkably reduced.
C     Victor P. Vysotskiy, 2014:
C     Number of 'ga_get' has been remarkably reduced by using the stripped mode

      Implicit None
      Integer irc
      Integer nSP_Batch
      Integer SP_BatchDim(nSP_Batch)
      Integer n_idSP
      Integer idSP(n_idSP)
      Integer l_NVT
      Integer NVT(l_NVT)
#include "cho_para_info.fh"
#include "choprint.fh"
#include "cholesky.fh"

      Real*8 C0, C1, W0, W1

      irc=0
      If (Cho_Real_Par) Then
         If (iPrint.ge.Inf_Pass) Then
            Call Cho_Timer(C0,W0)
         End If
         Call Cho_XCV_DV_P(irc,SP_BatchDim,nSP_Batch,idSP,n_idSP,
     &                     NVT,l_NVT)
         If (iPrint.ge.Inf_Pass) Then
            Call Cho_Timer(C1,W1)
            Write(LuPri,'(/,1X,A)')
     &      'Timing of vector distribution:'
            Call Cho_PrtTim(' ',C1,C0,W1,W0,-1)
         End If
      Else
         If (iPrint.ge.Inf_Pass) Then
            Call Cho_Timer(C0,W0)
         End If
         Call Cho_XCV_DV_S(irc,SP_BatchDim,nSP_Batch,idSP,n_idSP)
         If (iPrint.ge.Inf_Pass) Then
            Call Cho_Timer(C1,W1)
            Write(LuPri,'(/,1X,A)')
     &      'Timing of vector write:'
            Call Cho_PrtTim(' ',C1,C0,W1,W0,-1)
         End If
      End If

      End
      SubRoutine Cho_XCV_DV_P(irc,SP_BatchDim,nSP_Batch,
     &                        id_mySP,n_mySP,NVT,l_NVT)
      Implicit None
      Integer irc
      Integer nSP_Batch
      Integer SP_BatchDim(nSP_Batch)
      Integer n_mySP
      Integer id_mySP(n_mySP)
      Integer l_NVT
      Integer NVT(l_NVT)
#if defined (_MOLCAS_MPP_)
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"
#include "choprint.fh"
#include "mafdecls.fh"
#include "cho_para_info.fh"
      Character*12 SecNam
      Parameter (SecNam='Cho_XCV_DV_P')

      Integer iOpt
      Parameter (iOpt=1)

      Logical  ga_create, ga_destroy
      External ga_create, ga_destroy
      Logical  ok

      Real*8 X0, X1, Y0, Y1

      Integer g_a
      Integer iSym
      Integer ip_myNumCho, l_myNumCho
      Integer ip_Mem, l_Mem
      Integer max_vector_dim
      Integer nVec_per_batch, nVec_this_batch
      Integer iBatch, nBatch
      Integer ip_V, l_V
      Integer ip_numV, l_numV
      Integer ip_IDV, l_IDV
      Integer kV, myStart, myEnd
      Integer my_nV
      Integer lTot, iAdr, iAdr0, nDim
      Integer J0, J1, J2
#if !defined(_GA_)
      Integer Jst,Jen
#endif
      Integer iSP, iSP_, iSP1, iSP2
      Integer iSP_Batch
      Integer nSP_this_batch

      Integer i, j, k
      Integer iiBstRsh, nnBstRSh
      iiBStRsh(i,j,k)=iWork(ip_iiBstRSh-1+nSym*nnShl*(k-1)+nSym*(j-1)+i)
      nnBStRsh(i,j,k)=iWork(ip_nnBstRSh-1+nSym*nnShl*(k-1)+nSym*(j-1)+i)

      ! Init return code
      irc=0

      ! Find max vector dimension
      max_vector_dim=nnBstR(1,2)
      Do iSym=1,nSym
         max_vector_dim=max(max_vector_dim,nnBstR(iSym,2))
      End Do
      If (max_vector_dim .lt. 1) Then
         irc=-2
         Return
      End If

      ! Allocate my vector counter array
      l_myNumCho=nSym
      Call GetMem('GAMNCH','Allo','Inte',ip_myNumCho,l_myNumCho)
      Call iZero(iWork(ip_myNumCho),l_myNumCho)

      ! Allocate index arrays for GA use
      l_numV=1
      Call GetMem('GANUMV','Allo','Inte',ip_numV,l_numV)

      ! Figure out how much 1/3 of total memory is
      Call GetMem('GAMX','Max ','Real',ip_Mem,l_Mem)
      l_Mem=l_Mem/3
      If (l_Mem.lt.max_vector_dim) Then
         ! OK, let us try 2/3 then...
         l_Mem=l_Mem*2
         If (l_Mem.lt.max_vector_dim) Then
            irc=-1
            Return
         End If
      End If

      ! distribute vectors, one symmetry at a time
      Do iSym=1,nSym
         If (iPrint.ge.Inf_Progress) Then
            Write(LuPri,'(/,A,I2,/,A)')
     &      'Distributing vectors, symmetry',iSym,
     &      '--------------------------------'
            Write(LuPri,'(3X,A,I8)')
     &      'Total number of vectors:',NVT(iSym)
            Write(LuPri,'(3X,A,I8)')
     &      'Local number of vectors:',NumCho(iSym)
            Write(LuPri,'(3X,A,I8)')
     &      'Vector dimension       :',nnBstR(iSym,2)
            Write(LuPri,'(3X,A,I8)')
     &      'Shell pair batches     :',nSP_Batch
            Call Cho_Flush(LuPri)
         End If
         If (NVT(iSym).gt.0 .and. nnBstR(iSym,2).gt.0) Then
            ! Set up batching, ensuring that the number of vectors is
            ! the same on all nodes
            iWork(ip_numV)=min(l_Mem/nnBstR(iSym,2),NVT(iSym))
            Call Cho_GAIGOp(iWork(ip_numV),1,'min')
            nVec_per_batch=iWork(ip_numV)
            If (nVec_per_batch .lt. 1) Then
               Call Cho_Quit(
     &               'Insufficient memory for batching in '//SecNam,101)
            End If
            nBatch=(NVT(iSym)-1)/nVec_per_batch+1
            ! Allocate memory for vectors (will hold local as well as
            ! full vectors)
            l_V=nnBstR(iSym,2)*nVec_per_batch
            Call GetMem('GAVEC','Allo','Real',ip_V,l_V)
            ! Allocate vector ID array for distribution
            l_IDV=nVec_per_batch
            Call GetMem('GADIST','Allo','Inte',ip_IDV,l_IDV)
            ! Create global array with evenly distributed chunks
            ok=ga_create(mt_dbl,nnBstR(iSym,2),nVec_per_batch,'GA_XCV',
     &                   0,0,g_a)
            If (.not. ok) Then
               Call Cho_Quit(SecNam,': ga_create() failed!',101)
            End If
            If (iPrint.ge.Inf_Progress) Then
               Write(LuPri,'(3X,A,I8)')
     &         'Vector batches         :',nBatch
               Call Cho_Flush(LuPri)
            End If
            ! Distribute in batches
            Do iBatch=1,nBatch
               ! Sync: all must be at the same batch
               ! (so that we do not put while another is getting)
               Call Cho_GASync()
               ! Determine number of vectors in this batch
               If (iBatch.eq.nBatch) Then
                  nVec_this_batch=NVT(iSym)-nVec_per_batch*(nBatch-1)
               Else
                  nVec_this_batch=nVec_per_batch
               End If
               ! First and last vector in this batch
               J1=nVec_per_batch*(iBatch-1)+1
               J2=J1+nVec_this_batch-1
               If (iPrint.ge.Inf_Progress) Then
                  Write(LuPri,'(3X,A,I8,/,3X,A)')
     &            'Vector batch number:',iBatch,
     &            '++++++++++++++++++++++++++++'
                  Write(LuPri,'(6X,A,I8)')
     &            'Number of vectors in this batch:',nVec_this_batch
                  Write(LuPri,'(6X,A,I8,1X,I8)')
     &            'First and last vector          :',J1,J2
                  Call Cho_Flush(LuPri)
                  Call Cho_Timer(X0,Y0)
               End If
               ! Loop through SP blocks of vectors in this vector batch
               iAdr0=0
               iSP1=1
               Do iSP_Batch=1,nSP_Batch
                  nSP_this_batch=SP_BatchDim(iSP_Batch)
                  iSP2=iSP1+nSP_this_batch-1
                  nDim=0
                  Do iSP=iSP1,iSP2
                     nDim=nDim+nnBstRSh(iSym,id_mySP(iSP),2)
                  End Do
                  If (nDim.gt.0) Then
                     ! Read vector block
                     lTot=nDim*nVec_this_batch
                     iAdr=iAdr0+nDim*(J1-1)
                     Call DDAFile(LuTmp(iSym),2,Work(ip_V),lTot,iAdr)
                     iAdr0=iAdr0+nDim*NVT(iSym)
                     ! Put vector block into global array
                     kV=ip_V
                     Do iSP_=iSP1,iSP2
                        iSP=id_mySP(iSP_)
                        If (nnBstRSh(iSym,iSP,2).gt.0) Then
                           myStart=iiBstRSh(iSym,iSP,2)+1
                           myEnd=myStart+nnBstRSh(iSym,iSP,2)-1
C VPV:
                           Call ga_put(g_a,myStart,myEnd,1,
     &                          nVec_this_batch,Work(kV),nDim)
                           kV=kV+nnBstRSh(iSym,iSP,2)
                        End If
                     End Do
CCC ORIGINAL_CODE
C                     Do J=1,nVec_this_batch
C                        Do iSP_=iSP1,iSP2
C                           iSP=id_mySP(iSP_)
C                           If (nnBstRSh(iSym,iSP,2).gt.0) Then
C                              myStart=iiBstRSh(iSym,iSP,2)+1
C                              myEnd=myStart+nnBstRSh(iSym,iSP,2)-1
C                              write(6,*) "ga_put:", myStart,myEnd,J,kv,
C     &                        nnBstRSh(iSym,iSP,2),nDim
C                              Call ga_put(g_a,myStart,myEnd,J,J,
C     &                                    Work(kV),nnBstRSh(iSym,iSP,2))
C                              kV=kV+nnBstRSh(iSym,iSP,2)
C                           End If
C                        End Do
C                     End Do
CCC END OF ORIGINAL
                  End If
                  iSP1=iSP1+nSP_this_batch
               End Do
               If (iPrint.ge.Inf_Progress) Then
                  Call Cho_Timer(X1,Y1)
                  Write(LuPri,'(6X,A,F12.2,1X,F12.2)')
     &            'Time for read/ga_put (sec)     :',(X1-X0),(Y1-Y0)
                  Call Cho_Flush(LuPri)
               End If
               ! Sync: all must be done putting data into g_a
               Call Cho_GASync()
               ! Compute vector distribution for this batch
               my_nV=0
               Call Cho_P_Distrib_Vec(J1,J2,iWork(ip_IDV),my_nV)
               If (iPrint.ge.Inf_Progress) Then
                  Call Cho_Timer(X0,Y0)
               End If
               ! Get vectors from global array
               J0=J1-1
               kV=ip_V
#if defined(_GA_)
               Do i=1,my_nV
                  J=iWork(ip_IDV-1+i)-J0
                  Call ga_get(g_a,1,nnBstR(iSym,2),J,J,
     &                        Work(kV),nnBstR(iSym,2))
                  kV=kV+nnBstR(iSym,2)
               End Do
#else
               if(my_nV.gt.0) Then
                  Jst=iWork(ip_IDV)-J0
                  Jen=iWork(ip_IDV-1+my_nV)-J0
                  Call ga_get_striped(g_a,1,nnBstR(iSym,2),Jst,
     &                 Jen,Work(kV),nnBstR(iSym,2),nProcs)
                  kV=kV+nnBstR(iSym,2)*my_nV
               End If
               ! VVP: First performs RMA and only then I/O
               Call Cho_GASync()
#endif
               ! Write vectors to disk
               lTot=nnBstR(iSym,2)*my_nV
               If (lTot .gt. 0) Then
                  iAdr=nnBstR(iSym,2)*iWork(ip_myNumCho-1+iSym)
                  Call DDAFile(LuCho(iSym),iOpt,Work(ip_V),lTot,iAdr)
               End If
               If (iPrint.ge.Inf_Progress) Then
                  Call Cho_Timer(X1,Y1)
                  Write(LuPri,'(6X,A,F12.2,1X,F12.2)')
     &            'Time for ga_get/write (sec)    :',(X1-X0),(Y1-Y0)
                  Call Cho_Flush(LuPri)
               End If
               ! Update my vector counter
               iWork(ip_myNumCho-1+iSym)=iWork(ip_myNumCho-1+iSym)+my_nV
            End Do
            ! Destroy global array
            ok = ga_destroy(g_a)
            ! Deallocations
            Call GetMem('GADIST','Free','Inte',ip_IDV,l_IDV)
            Call GetMem('GAVEC','Free','Real',ip_V,l_V)
         End If
      End Do

      ! Deallocations
      Call GetMem('GANUMV','Free','Inte',ip_numV,l_numV)
      Call GetMem('GAMNCH','Free','Inte',ip_myNumCho,l_myNumCho)
#else
      irc=999 ! should never be called in serial installation
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(SP_BatchDim)
         Call Unused_integer(id_mySP)
         Call Unused_integer(NVT)
      End If
#endif

      End
      SubRoutine Cho_XCV_DV_S(irc,SP_BatchDim,nSP_Batch,id_mySP,n_mySP)
      Implicit None
      Integer irc
      Integer nSP_Batch
      Integer SP_BatchDim(nSP_Batch)
      Integer n_mySP
      Integer id_mySP(n_mySP)
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"
#include "choprint.fh"

      Character*12 SecNam
      Parameter (SecNam='Cho_XCV_DV_S')

      Real*8 X0, X1, Y0, Y1

      Integer iSym, iSP, iSP1, iSP2, J1, nDim
      Integer ip_Mem, l_Mem
      Integer ip_V, ip_T
      Integer max_vector_dim, max_block_dim
      Integer nVec_per_batch, nVec_this_batch
      Integer iBatch, nBatch
      Integer kV, kT, kOffV, kOffT
      Integer lTot, iAdr, iAdr0
      Integer iSP_Batch, nSP_this_batch

      Integer i, j, k
      Integer iiBstRsh, nnBstRSh
      iiBStRsh(i,j,k)=iWork(ip_iiBstRSh-1+nSym*nnShl*(k-1)+nSym*(j-1)+i)
      nnBStRsh(i,j,k)=iWork(ip_nnBstRSh-1+nSym*nnShl*(k-1)+nSym*(j-1)+i)

      ! Init return code
      irc=0

      ! Find max vector dimension
      max_vector_dim=nnBstR(1,2)
      Do iSym=2,nSym
         max_vector_dim=max(max_vector_dim,nnBstR(iSym,2))
      End Do
      If (max_vector_dim .lt. 1) Then
         irc=-2
         Return
      End If

      ! Find max block dimension
      max_block_dim=0
      iSP1=1
      Do iSP_Batch=1,nSP_Batch
         nSP_this_batch=SP_BatchDim(iSP_Batch)
         iSP2=iSP1+nSP_this_batch-1
         Do iSym=1,nSym
            nDim=0
            Do iSP=iSP1,iSP2
               nDim=nDim+nnBstRSh(iSym,id_mySP(iSP),2)
            End Do
            max_block_dim=max(max_block_dim,nDim)
         End Do
         iSP1=iSP1+nSP_this_batch
      End Do
#if defined (_DEBUG_)
      If ((iSP1-1).ne.n_mySP) Then
         Call Cho_Quit(SecNam//': SP batch dimension error',103)
      End If
#endif

      ! Get largest memory block
      Call GetMem('DVSMX','Max ','Real',ip_Mem,l_Mem)
      lTot=max_vector_dim+max_block_dim
      If (l_Mem .lt. lTot) Then
         irc=-1
         Return
      End If
      Call GetMem('DVSVEC','Allo','Real',ip_Mem,l_Mem)

      ! read, reorder, and write vectors, one symmetry at a time
      Do iSym=1,nSym
         If (iPrint.ge.Inf_Progress) Then
            Write(LuPri,'(/,A,I2,/,A)')
     &      'Writing vectors, symmetry',iSym,
     &      '---------------------------'
            Write(LuPri,'(3X,A,I8)')
     &      'Total number of vectors:',NumCho(iSym)
            Write(LuPri,'(3X,A,I8)')
     &      'Vector dimension       :',nnBstR(iSym,2)
            Write(LuPri,'(3X,A,I8)')
     &      'Shell pair batches     :',nSP_Batch
            Call Cho_Flush(LuPri)
         End If
         If (NumCho(iSym).gt.0 .and. nnBstR(iSym,2).gt.0) Then
            ! Set up batching
            lTot=max_block_dim+nnBstR(iSym,2)
            nVec_per_batch=min(l_Mem/lTot,NumCho(iSym))
            If (nVec_per_batch .lt. 1) Then
               Call Cho_Quit(
     &               'Insufficient memory for batching in '//SecNam,101)
            End If
            nBatch=(NumCho(iSym)-1)/nVec_per_batch+1
            If (iPrint.ge.Inf_Progress) Then
               Write(LuPri,'(3X,A,I8)')
     &         'Vector batches         :',nBatch
               Call Cho_Flush(LuPri)
            End If
            ! Read and write vectors in batches
            Do iBatch=1,nBatch
               ! Determine number of vectors in this batch
               If (iBatch.eq.nBatch) Then
                  nVec_this_batch=NumCho(iSym)-nVec_per_batch*(nBatch-1)
               Else
                  nVec_this_batch=nVec_per_batch
               End If
               ! First vector in this batch
               J1=nVec_per_batch*(iBatch-1)+1
               If (iPrint.ge.Inf_Progress) Then
                  Write(LuPri,'(3X,A,I8,/,3X,A)')
     &            'Vector batch number:',iBatch,
     &            '++++++++++++++++++++++++++++'
                  Write(LuPri,'(6X,A,I8)')
     &            'Number of vectors in this batch:',nVec_this_batch
                  Write(LuPri,'(6X,A,I8,1X,I8)')
     &            'First and last vector          :',
     &            J1,J1+nVec_this_batch-1
                  Call Cho_Flush(LuPri)
                  Call Cho_Timer(X0,Y0)
               End If
               ! Set memory pointers
               ip_T=ip_Mem
               ip_V=ip_T+max_block_dim*nVec_this_batch
               ! Read blocked vectors and copy into full vector array
               iAdr0=0
               iSP1=1
               Do iSP_Batch=1,nSP_Batch
                  nSP_this_batch=SP_BatchDim(iSP_Batch)
                  iSP2=iSP1+nSP_this_batch-1
                  nDim=0
                  Do iSP=iSP1,iSP2
                     nDim=nDim+nnBstRSh(iSym,id_mySP(iSP),2)
                  End Do
                  If (nDim.gt.0) Then
                     lToT=nDim*nVec_this_batch
                     iAdr=iAdr0+nDim*(J1-1)
                     Call DDAFile(LuTmp(iSym),2,Work(ip_T),lTot,iAdr)
                     kT=ip_T-1
                     kV=ip_V-1
                     Do J=1,nVec_this_batch
                        kOffT=kT+nDim*(J-1)
                        kOffV=kV+nnBstR(iSym,2)*(J-1)
     &                          +iiBstRSh(iSym,id_mySP(iSP1),2)
                        Do i=1,nDim
                           Work(kOffV+i)=Work(kOffT+i)
                        End Do
                     End Do
                     iAdr0=iAdr0+nDim*NumCho(iSym)
                  End If
                  iSP1=iSP1+nSP_this_batch
               End Do
               If (iPrint.ge.Inf_Progress) Then
                  Call Cho_Timer(X1,Y1)
                  Write(LuPri,'(6X,A,F12.2,1X,F12.2)')
     &            'Time for read/reorder (sec)    :',(X1-X0),(Y1-Y0)
                  Call Cho_Flush(LuPri)
               End If
               ! Write full vectors to disk
               lTot=nnBstR(iSym,2)*nVec_this_batch
               iAdr=nnBstR(iSym,2)*(J1-1)
               Call DDAFile(LuCho(iSym),1,Work(ip_V),lTot,iAdr)
               If (iPrint.ge.Inf_Progress) Then
                  Call Cho_Timer(X0,Y0)
                  Write(LuPri,'(6X,A,F12.2,1X,F12.2)')
     &            'Time for write (sec)           :',(X0-X1),(Y0-Y1)
                  Call Cho_Flush(LuPri)
               End If
            End Do
         End If
      End Do

      ! Deallocation
      Call GetMem('DVSVEC','Free','Real',ip_Mem,l_Mem)

      End
