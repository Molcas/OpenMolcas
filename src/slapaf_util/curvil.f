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
* Copyright (C) 2004, Roland Lindh                                     *
************************************************************************
      Subroutine CurviL(nAtoms,nDim,Cx,Gx,nIter,iIter,iRef,nStab,
     &                  jStab,Degen,Smmtrc,mTR,TRVec,
     &                  ip_rInt,ip_drInt,HSet,BSet,ipBMx,Numerical,iANr,
     &                  HWRS,Analytic_Hessian,iOptC,Name,PrQ,
     &                  dMass,iCoSet,iTabBonds,
     &                  iTabAtoms,nBonds,nMax,iTabAI,mAtoms,lOld,
     &                  ip_KtB_Hessian,nQQ,nqInt,MaxItr,nWndw)
************************************************************************
*                                                                      *
*     Objective: to handle curvilinear internal coordinates.           *
*                                                                      *
*                                                                      *
*     Authors: R. Lindh, Dept. of Theoretical Chemistry                *
*              University of Lund, SWEDEN.                             *
*              2004                                                    *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "warnings.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "db.fh"
#include "print.fh"
      Real*8 Cx(3*nAtoms,nIter), Degen(3*nAtoms), Gx(3*nAtoms,nIter),
     &       dMass(nAtoms), TRVec(nDim,mTR)
      Integer nStab(nAtoms), iCoSet(0:7,nAtoms), jStab(0:7,nAtoms),
     &        iANr(nAtoms), iDum(6), iTabBonds(3,nBonds),
     &        iTabAtoms(0:nMax,nAtoms), iTabAI(2,mAtoms)
      Logical Smmtrc(3*nAtoms), HSet, BSet, Proc, Numerical, HWRS,
     &        Analytic_Hessian, PrQ, Proc_dB, lOld, Proc_H
      Character(LEN=32) filnam
      Character(LEN=LENIN) Name(nAtoms)
      Character(LEN=14) cDum
      Real*8 Dum(1)
      Logical, Save:: g12K=.False.
      Real*8, Allocatable:: Proj(:), Temp2(:), KtM(:,:), Degen2(:),
     &                      EVal(:), G(:), GxR(:,:), qVal(:,:),
     &                      F_c(:), KtB(:), K(:), GRef(:), Mult(:)
      Character(LEN=14), Allocatable:: qLbl(:)
      Integer, Allocatable:: Ind(:,:)
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      iPrint=99
#else
      iRout=128
      iPrint=nPrint(iRout)+1
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_allocate(Proj,nDim,Label='Proj')
*                                                                      *
************************************************************************
*                                                                      *
      ip_KtB_Hessian= ip_Dummy
      Thr_raw=3.0D-2
      If (HWRS) Then
         Thr_ElRed=Thr_raw**2
      Else
         Thr_ElRed=Thr_raw
      End If
*
      i=0
      Do iX = 1, 3*nAtoms
         If (Smmtrc(iX)) Then
            i = i + 1
            Proj(i)=One/Degen(iX)
         End If
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*.... First some words about the handling of the symmetry. In this code
*     I have selected to use only the unique components of the vectors.
*     Hence, some care has to be taken to account for the degeneracy of
*     some of the elements. This is handled by inserting the u matrix
*     at appropriate places.
*
      nQQ=nDim-mTR
      If (ip_rInt.eq.ip_Dummy) Then
         nqInt=nQQ*MaxItr
         Call GetMem(' qInt','Allo','Real',ip_rInt, nqInt)
         Call GetMem('dqInt','Allo','Real',ip_drInt,nqInt)
         Call FZero(Work(ip_rInt),nqInt)
         Call FZero(Work(ip_drInt),nqInt)
      End If
*
      Call mma_allocate(Degen2,nDim)
      i=0
      Do ix = 1, 3*nAtoms
         If (Smmtrc(ix)) Then
            i = i + 1
            Degen2(i) = Degen(ix)
         End If
      End Do
#ifdef _DEBUGPRINT_
      Call RecPrt('Degen2',' ',Degen2,nDim,1)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      filnam='INTCOR'
      LuIC=isfreeunit(4)
      call molcas_open(LuIC,filnam)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     nq: Number of redundant internal coordinates.
*
      Proc=.False.     ! Flag for processing B
      Proc_dB=.False.  ! Flag for processing dB
*
      Thr_small=(30.0D0/180.0D0)*Pi
      Do While (Thr_small.gt.1.0D-6)
         Call Get_Curvil
     &          (nq,nqRF,nqB,nqA,nqT,nqO,
     &           nAtoms,iIter,nIter,Cx,jStab,
     &           nStab,Smmtrc,Proc,Dum,1,iANr,cDum,
     &           iRef,Dum,Dum,iOptC,LuIC,
     &           Name,iDum,iIter,dMass,iCoSet,Dum,
     &           iDum(1),iDum(1),
     &           Proc_dB,
     &           iTabBonds,iTabAtoms,nBonds,nMax,iTabAI,mAtoms,
     &           mB_Tot,mdB_Tot,Dum,Dum,iDum,iDum,1,1,
     &           iDum,Thr_small)
         If (nq.ge.nQQ) Exit
         Thr_small=Thr_small-(5.0D0/180.0D0)*Pi
      End Do

      Rewind(LuIC)
*
      If (nq.eq.0) Then
         Call WarningMessage(2,' Curvil: nq.eq.0')
         Call Quit(_RC_INTERNAL_ERROR_)
      End If
*
      iGlow=1+nqRF+nqB
      iGhi =nqRF+nqB+nqA
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      Write (6,*) 'nq, nqB, nqA, nqT, nqO=',
     &             nq, nqB, nqA, nqT, nqO
#endif
*
*---- Now allocate some arrays which depend on nq
*
      If (ip_B.ne.ip_Dummy) Then
         Call Free_Work(ip_B)
         Call Free_iWork(ip_iB)
         Call Free_iWork(ip_nqB)
      End If
      Call GetMem('BM','Allo','Real',ip_B,mB_Tot)
      Call GetMem('iBM','Allo','Inte',ip_iB,mB_Tot)
      Call GetMem('nqB','Allo','Inte',ip_nqB,nq)
      mq=nq
*
      Call mma_allocate(qVal,nq,nIter,Label='qVal')
      Call mma_allocate(qLbl,nq,Label='qLbl')
      Call mma_allocate(F_c,nq,Label='F_c')
      Call mma_allocate(Mult,nq,Label='Mult')
      Call mma_allocate(Ind,3,nq,Label='Ind')
      Call mma_allocate(GRef,9*nqA*nIter,Label='GRef')
*
      F_c(:) = Zero
      Mult(:)= Zero
      qVal(:,:) = Zero
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*---- Process the redundant internal coordinates for the current
*     structure.
*
      Proc=.True.
      Proc_dB=.False.
*
      Call Get_Curvil
     &          (iq,iqRF,iqR,iqA,iqT,iqO,
     &           nAtoms,iIter,nIter,Cx,jStab,
     &           nStab,Smmtrc,Proc,
     &           qVal,nq,iANr,qLbl,
     &           iRef,F_c,Mult,iOptC,
     &           LuIC,Name,Ind,iIter,dMass,iCoSet,GRef,
     &           iGlow,iGHi,
     &           Proc_dB,
     &           iTabBonds,iTabAtoms,nBonds,nMax,iTabAI,mAtoms,
     &           nB_Tot,ndB_Tot,
     &           Work(ip_B),Dum,iWork(ip_iB),iDum,
     &           mB_Tot,mdB_Tot,iWork(ip_nqB),Thr_small)
      Rewind(LuIC)
*
      If (iq.ne.nq) Then
         Call WarningMessage(2,' Error in Curvil')
         Write (6,*) 'In Curvil: iq.ne.nq'
         Write (6,*) 'iq=',iq
         Write (6,*) 'nq=',nq
         Call Abend
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      If (iPrint.ge.49) Then
         Write (6,*) 'nq, nqB, nqA, nqT, nqO=',
     &             nq, nqB, nqA, nqT, nqO
         Call RecPrt('q-values',' ',qVal,nq,nIter)
         Call RecPrt('Force Constant matrix in redundant basis',' ',
     &               F_c,1,nq)
         Call RecPrt('Multiplicity factors',' ',Mult,1,nq)
         Call RecPrt('Cx',' ',Cx,3*nAtoms,nIter)
         Call RecPrt('Gx',' ',Gx,3*nAtoms,nIter)
      End If
#endif
*
*     Notation:
*     X: cartesian coordinates
*     q: redundant internal coordinates
*     Q: nonredundant internal coordinates
*     u: the degeneracy matrix
*
*     dq = B u dx
*
*---- Start processing the B matrix for the current structure to
*     generate the K matrix (dQ/dq).
*
      Call mma_allocate(K,nq**2,Label='K')
      Call mma_allocate(KtB,nDim**2,Label='KtB')
*                                                                      *
************************************************************************
*                                                                      *
*---- Stick in the metric of the force constants in the redundant
*     space.
*
      If (HWRS) Then
*------- Scale each coordinate with the force constant
*
         i = 0
         Do iq = 0, nq-1
            nB = iWork(ip_nqB+iq)
            Call DScal_(nB,F_c(1+iq),Work(ip_B+i),1)
            i = i + nB
         End Do
*ifdef _DEBUGPRINT_
         If (iPrint.ge.99) Then
            i = 0
            Do iq = 0, nq-1
               nB = iWork(ip_nqB+iq)
               Call RecPrt('fcB',' ',Work(ip_B+i),1,nB)
               i = i + nB
            End Do
         End If
*endif
      End If
*                                                                      *
************************************************************************
*                                                                      *
*---- Eliminate redundancy to produce the nonredundant internal
*     coordinates.
*
*           t        t
*     dQ = K M dq = K M B u dx
*
      i = 0
      Do iq = 0, nq-1
         nB = iWork(ip_nqB+iq)
         Call DScal_(nB,Mult(1+iq),Work(ip_B+i),1)
         i = i + nB
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      If (BSet) Then
         Call mma_allocate(G,nq*nq,Label='G')
         Call mma_allocate(EVal,nq*(nq+1)/2,Label='EVal')
         Call ElRed2(nq,ndim,G,EVal,K,nK,Proj,
     &               g12K,Thr_ElRed,Work(ip_B),iWork(ip_iB),mB_Tot,
     &               iWork(ip_nqB))

         If (nK.gt.nQQ) Then
            Call Remove_TR(nq,ndim,nQQ,K,nK,TRVec,mTR,
     &                     Work(ip_B),iWork(ip_iB),iWork(ip_nqB),mB_Tot)
         End If
         Call mma_deallocate(EVal)
         Call mma_deallocate(G)
*
         Call Put_dArray('K',K,nq*nQQ)
      Else
         Call Get_dArray('K',K,nq*nQQ)
      End If
#ifdef _DEBUGPRINT_
      If (iPrint.ge.99) Then
         Call RecPrt('K',' ',K,nq,nQQ)
      End If
#endif
*                                                                      *
************************************************************************
*                                                                      *
*---- Print the nonreduntant coordinates.
*
      If (PrQ) Then
         Call mma_allocate(Temp2,nq*nQQ,Label='Temp2')
         call dcopy_(nQQ*nq,K,1,Temp2,1)
         Do iq = 1, nq
            Call DScal_(nQQ,One/Mult(iq),Temp2(iq),nq)
         End Do
         Call PrintQ(Temp2,qLbl,nq,nQQ,LuIC,Mult)
         Call mma_deallocate(Temp2)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (PrQ.and.iPrint.ge.6) Then
         Write (6,*)
         Write (6,'(A)')   ' ******************************************'
         Write (6,'(A)')   ' * Statistics of the internal coordinates *'
         Write (6,'(A)')   ' ******************************************'
         Write (6,'(A,I5)')' Translations and Rotations:     ',nqRF
         Write (6,'(A,I5)')' Bonds                     :     ',nqB
         Write (6,'(A,I5)')' Angles                    :     ',nqA
         Write (6,'(A,I5)')' Torsions                  :     ',nqT
         Write (6,'(A,I5)')' Out-of-plane angles       :     ',nqO
         Write (6,*)
      End If
      If (nq.lt.nQQ) Then
         Call WarningMessage(2,' Error in Curvil')
         Write (6,*) 'In Curvil: nq.lt.nQQ'
         Write (6,*) 'nq=',nq
         Write (6,*) 'nQQ=',nQQ
         Call Abend
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     The nonredundant internal coordinates are now defined!           *
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*---- Compute the values and the gradients of the nonredundant internal
*     curvilinear coordinates for all the iterations.
*
*     Loop over
*     1) the last to the first iteration if normal execusion
*     2) only the last if in transformation from internal to
*        cartesian
*
      If (BSet) Then
*------- Produce list of compressed cartesian gradients.
         Call mma_allocate(GxR,nDim,nIter,Label='GxR')
         Do jIter = 1, nIter
            Call NRed(Gx(:,jIter),GxR(:,jIter),3*nAtoms,nDim,Smmtrc)
         End Do
*
         iSt = nIter
         iEnd = iSt - Min(nIter,nWndw+1) + 1
C        iEnd = 1
         iOff = nIter
      Else
         If (Numerical) Then
            Call mma_allocate(GxR,nDim,1,Label='GxR')
            iOff = 1
            Call NRed(Gx(1,nIter),GxR(:,1),3*nAtoms,nDim,Smmtrc)
         Else
            Call mma_allocate(GxR,1,1,Label='GxR')
            iOff = 1
         End If
         iEnd=nIter
         iSt =nIter
      End If
*
*     Note that the loop is in reverse order.
*
      Do jIter = iSt, iEnd, -1
         Proc_H=HSet.and.jIter.eq.iRef.and..Not.lOld
*        iOptC(256) = constrained optimization
         Proc_dB=Proc_H.and.
     &           (Analytic_Hessian.or.Numerical.or.
     &            iAnd(iOptC,256).eq.256)
*        Compute and store dBQQ in the reference structure
         If (Proc_dB) Then
            If (ip_dB.ne.ip_Dummy) Then
               Call Free_Work(ip_dB)
               Call Free_iWork(ip_idB)
            End If
            Call GetMem('dBM','Allo','Real',ip_dB,mdB_Tot)
            Call GetMem('idBM','Allo','Inte',ip_idB,mdB_Tot*2)
         End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
         Call Get_Curvil
     &             (iq,iqRF,iqR,iqA,iqT,iqO,
     &              nAtoms,jIter,nIter,Cx,jStab,
     &              nStab,Smmtrc,Proc,
     &              qVal,nq,iANr,qLbl,
     &              iRef, F_c,Mult,
     &              iOptC,LuIC,Name,Ind,iIter,dMass,iCoSet,
     &              GRef,iGlow,iGHi,
     &              Proc_dB,
     &              iTabBonds,iTabAtoms,nBonds,nMax,iTabAI,mAtoms,
     &              nB_Tot,ndB_Tot,
     &              Work(ip_B),Work(ip_dB),iWork(ip_iB),iWork(ip_idB),
     &              mB_Tot,mdB_Tot,iWork(ip_nqB),Thr_small)
         Rewind(LuIC)
*
         If (iq.ne.nq) Then
            Write (6,*) 'In Curvil: iq.ne.nq'
            Write (6,*) 'iq=',iq
            Write (6,*) 'nq=',nq
            Call Abend
         End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*------- Form the gradients of Q
*                                         T
*------- Form the K(t)MB matrix for dQ = K M B u dx
*
         i = 0
         Do iq = 0, nq-1
            nB = iWork(ip_nqB+iq)
            Call DScal_(nB,Mult(1+iq),Work(ip_B+i),1)
            i = i + nB
         End Do
         Ktb(1:nQQ*nDim)=Zero
         Do iQQ = 0, nQQ-1
            i = 0
            Do iq = 0, nq-1
               nB = iWork(ip_nqB+iq)
               Do iB = 0, nB-1
                  iDim = iWork(ip_iB+i)
                  KtB((iDim-1)*nQQ + iQQ + 1) =
     &                 KtB((iDim-1)*nQQ + iQQ + 1)
     &               + K(iQQ*nq+iq+1)
     &               * Work(i+ip_B)
                  i = i + 1
               End Do
            End Do
         End Do
#ifdef _DEBUGPRINT_
         If (iPrint.ge.99) Then
            Call RecPrt(' The K matrix',' ',K,nq,nQQ)
            Call RecPrt(' The K(t)B matrix',' ',KtB,nQQ,nDim)
         End If
#endif
*                                                                      *
************************************************************************
*                                                                      *
*------- Store the B matrix for the last structure
*
         If (jIter.eq.nIter) Then
            nX = 3*nAtoms
            Call Allocate_Work(ipBmx,nX**2)
            Call FZero(Work(ipBMx),nX**2)
*
*           modify from compact to full Cartesian storage.
*
            iDim = 0
            Do iX = 1, nX
               If (Smmtrc(iX)) Then
                  iDim = iDim + 1
                  Do iQQ = 1, nQQ
                     iXQ = (iQQ-1)*nX + iX + ipBMx - 1
                     iQD = (iDim-1)*nQQ + iQQ
                     Work(iXQ) = KtB(iQD)
                  End Do
               Else
                  Do iQQ = 1, nQQ
                     iXQ = (iQQ-1)*nX + iX + ipBMx - 1
                     Work(iXQ) = Zero
                  End Do
               End If
            End Do

         End If
*                                                                      *
************************************************************************
*                                                                      *
*------- Branch out if only values are to be computed.
*
         If (BSet .or. Numerical) Then
*                                                                      *
************************************************************************
*                                                                      *
*---------- Form the gradient for iteration jIter for the new definition
*           of the K matrix.
*
*           dq/dx dE/dq = dE/dq
*
*           The B-matrix is stored (3*natom x nQQ)
*           KtB is stored nQQ, nDim
*
            Call Allocate_Work(ip_KtB,nDim*nQQ)
            Call TRNSPS(nQQ,nDim,KtB,Work(ip_KtB))
*
*           Strip KtB of the degeneracy factor (full).
*
            Do iQQ = 1, nQQ
               Do iDim = 1, nDim
                  ij = (iQQ-1)*nDim + iDim -1 + ip_KtB
                  Work(ij) = Work(ij) / Degen2(iDim)
               End Do
            End Do
*
            ip = ip_drInt + (jIter-1)*nQQ
            M = nDim
            N = nQQ
            NRHS=1
            Call Eq_Solver('N',M,N,NRHS,Work(ip_KtB   ),.False.,
     &                     Degen2,GxR(:,iOff),Work(ip))
*           Call RecPrt('GxR(:,iSt',' ',GxR(:,iOff),nDim,1)
*           Call RecPrt('KtB   ',' ',KtB,nQQ,nDim)
*           Call RecPrt('drInt',' ',Work(ip),nQQ,1)

            iOff = iOff - 1
*                                                                      *
************************************************************************
*                                                                      *
*---------- Save pointer to KtB to be used in backtransforming the
*           Hessian to internal coordinates.
*
            If (Proc_H) Then
               ip_KtB_Hessian=ip_KtB
*              Write (*,*) 'ip_KtB_Hessian=',ip_KtB_Hessian
*              Call RecPrt('KtB',' ',Work(ip_KtB),nDim,nQQ)
            Else
               Call Free_Work(ip_KtB)
            End If
*                                                                      *
************************************************************************
*                                                                      *
         End If
*                                                                      *
************************************************************************
*                                                                      *
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*---- Process the values
*          t
*     Q = K M q
*
      jIter = 1
      mIter = nIter
      If (.Not.BSet) Then
         jIter = nIter
         mIter = 1
      End If
      Call mma_allocate(KtM,nQQ,nq,Label='KtM')
      Do iq = 1, nq
         Alpha=Mult(iq)
         Do iQQ = 1, nQQ
            temp = K(iq+(iQQ-1)*nq)*Alpha
            KtM(iQQ,iq)=temp
         End Do
      End Do
      ip = ip_rInt + (jIter-1)*(nDim-mTR)
      Call DGEMM_('N','N',
     &            nQQ,mIter,nq,
     &            1.0d0,KtM,nQQ,
     &                  qVal(:,jIter),nq,
     &            0.0d0,Work(ip),nQQ)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      If (iPrint.ge.49) Then
         Call RecPrt(' The K Matrix',' ',K,nq,nQQ)
         Call RecPrt(' q-values',' ',qVal,nq,nIter)
         Call RecPrt('Q-values',' ',Work(ip_rInt),nQQ,nIter)
         Call RecPrt('Cx',' ',Cx,3*nAtoms,nIter)
      End If
      If (BSet.and.iPrint.ge.49) Then
          Call RecPrt('Q-gradients',' ',Work(ip_drInt),nQQ,nIter)
          Call RecPrt('Gx',' ',Gx,3*nAtoms,nIter)
      End If
#endif
*                                                                      *
************************************************************************
*                                                                      *
*---- Deallocate memory
*
      Call mma_deallocate(KtM)
      If (Allocated(GxR)) Call mma_deallocate(GxR)
      Call mma_deallocate(KtB)
      Call mma_deallocate(K)
      Call mma_deallocate(GRef)
      Call mma_deallocate(Ind)
      Call mma_deallocate(Mult)
      Call mma_deallocate(F_c)
      Call mma_deallocate(qLbl)
      Call mma_deallocate(qVal)
*
      Call mma_deallocate(Degen2)
      Call mma_deallocate(Proj)
*
      Close (LuIC)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
