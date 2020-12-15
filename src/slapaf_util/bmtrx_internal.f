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
      Subroutine BMtrx_Internal(nAtoms,nDimBC,
     &                          nIter,mAtoms,Numerical,
     &                          iIter,mTR,TRVec,iTabAI,iTabAtoms,
     &                          iTabBonds,nBonds,nMax,iRef,nQQ,nWndw)
************************************************************************
*                                                                      *
*     Objective: to handle curvilinear internal coordinates.           *
*                                                                      *
*                                                                      *
*     Authors: R. Lindh, Dept. of Theoretical Chemistry                *
*              University of Lund, SWEDEN.                             *
*              2004                                                    *
************************************************************************
      use Slapaf_Info, only: qInt, dqInt, BM, dBM, iBM, idBM, nqBM, KtB,
     &                       Cx, Gx, BMx, Degen, Smmtrc
      use Slapaf_Parameters, only: HWRS, Analytic_Hessian, MaxItr,
     &                             iOptC, BSet, HSet, PrQ, lOld
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "warnings.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "db.fh"
#include "print.fh"
      Integer, Intent(In):: nAtoms, nDimBC
      Integer, Intent(In):: nIter, mAtoms
      Logical, Intent(In):: Numerical
      Integer, Intent(In):: iIter, mTR
      Real*8, Intent(In):: TRVec(nDimBC,mTR)
      Integer, Intent(In):: iTabBonds(3,nBonds),
     &                      iTabAtoms(0:nMax,nAtoms),
     &                      iTabAI(2,mAtoms)
      Integer, Intent(In):: nBonds,nMax,iRef
      Integer, Intent(InOut):: nQQ
      Integer, Intent(In):: nWndW

      Integer iDum(6)
      Logical Proc, Proc_dB, Proc_H
      Character(LEN=32) filnam
      Character(LEN=14) cDum
      Real*8 Dum(1)
      Logical, Save:: g12K=.False.
      Real*8, Allocatable:: Proj(:), Temp2(:), KtM(:,:), Degen2(:),
     &                      EVal(:), G(:), GxR(:,:), qVal(:,:),
     &                      F_c(:), K(:), GRef(:), Mult(:)
      Real*8, Allocatable:: KtBu(:), KtBt(:,:)
      Character(LEN=14), Allocatable:: qLbl(:)
      Integer, Allocatable:: Ind(:,:)
*                                                                      *
************************************************************************
*                                                                      *
*#define _DEBUGPRINT_
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
      Call mma_allocate(Proj,nDimBC,Label='Proj')
*                                                                      *
************************************************************************
*                                                                      *
      Thr_raw=3.0D-2
      If (HWRS) Then
         Thr_ElRed=Thr_raw**2
      Else
         Thr_ElRed=Thr_raw
      End If
*
      i=0
      Do iX = 1, 3*nAtoms
         iAtom = (iX+2)/3
         ixyz = iX - (iAtom-1)*3
         If (Smmtrc(ixyz,iAtom)) Then
            i = i + 1
            Proj(i)=One/Degen(ixyz,iAtom)
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
      nQQ=nDimBC-mTR
      If (.NOT.Allocated(qInt)) Then
         Call mma_allocate( qInt,nQQ,MaxItr,Label=' qInt')
         Call mma_allocate(dqInt,nQQ,MaxItr,Label='dqInt')
          qInt(:,:)=Zero
         dqInt(:,:)=Zero
      End If
*
      Call mma_allocate(Degen2,nDimBC)
      i=0
      Do ix = 1, 3*nAtoms
         iAtom = (ix+2)/3
         ixyz = ix - (iAtom-1)*3
         If (Smmtrc(ixyz,iAtom)) Then
            i = i + 1
            Degen2(i) = Degen(ixyz,iAtom)
         End If
      End Do
#ifdef _DEBUGPRINT_
      Call RecPrt('Degen2',' ',Degen2,nDimBC,1)
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
     &           nAtoms,iIter,nIter,Cx,
     &           Proc,Dum,1,cDum,
     &           iRef,Dum,Dum,LuIC,
     &           iDum,iIter,Dum,
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
      If (Allocated(BM)) Call mma_deallocate(BM)
      If (Allocated(iBM)) Call mma_deallocate(iBM)
      If (Allocated(nqBM)) Call mma_deallocate(nqBM)
      Call mma_allocate(BM,mB_Tot,Label='BM')
      Call mma_allocate(iBM,mB_Tot,Label='iBM')
      Call mma_allocate(nqBM,nq,Label='nqBM')
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
     &           nAtoms,iIter,nIter,Cx,
     &           Proc,
     &           qVal,nq,qLbl,
     &           iRef,F_c,Mult,
     &           LuIC,Ind,iIter,GRef,
     &           iGlow,iGHi,
     &           Proc_dB,
     &           iTabBonds,iTabAtoms,nBonds,nMax,iTabAI,mAtoms,
     &           nB_Tot,ndB_Tot,
     &           BM,Dum,iBM,iDum,
     &           mB_Tot,mdB_Tot,nqBM,Thr_small)
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
      Call mma_allocate(KtBu,nDimBC**2,Label='KtBu')
*                                                                      *
************************************************************************
*                                                                      *
*---- Stick in the metric of the force constants in the redundant
*     space.
*
      If (HWRS) Then
*------- Scale each coordinate with the force constant
*
         i = 1
         Do iq = 1, nq
            nB = nqBM(iq)
            Call DScal_(nB,F_c(iq),BM(i),1)
            i = i + nB
         End Do
*ifdef _DEBUGPRINT_
         If (iPrint.ge.99) Then
            i = 1
            Do iq = 1, nq
               nB = nqBM(iq)
               Call RecPrt('fcB',' ',BM(i),1,nB)
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
      i = 1
      Do iq = 1, nq
         nB = nqBM(iq)
         Call DScal_(nB,Mult(iq),BM(i),1)
         i = i + nB
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      If (BSet) Then
         Call mma_allocate(G,nq*nq,Label='G')
         Call mma_allocate(EVal,nq*(nq+1)/2,Label='EVal')
         Call ElRed2(nq,nDimBC,G,EVal,K,nK,Proj,
     &               g12K,Thr_ElRed,BM,iBM,mB_Tot,nqBM)

         If (nK.gt.nQQ) Then
            Call Remove_TR(nq,nDimBC,nQQ,K,nK,TRVec,mTR,
     &                     BM,iBM,nqBM,mB_Tot)
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
         Call mma_allocate(GxR,nDimBC,nIter,Label='GxR')
         Do jIter = 1, nIter
            Call NRed(Gx(:,:,jIter),GxR(:,jIter),3*nAtoms,nDimBC,Smmtrc)
         End Do
*
         iSt = nIter
         iEnd = iSt - Min(nIter,nWndw+1) + 1
C        iEnd = 1
         iOff = nIter
      Else
         If (Numerical) Then
            Call mma_allocate(GxR,nDimBC,1,Label='GxR')
            iOff = 1
            Call NRed(Gx(:,:,nIter),GxR(:,1),3*nAtoms,nDimBC,Smmtrc)
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
            If (Allocated(dBM)) Call mma_deallocate(dBM)
            If (Allocated(idBM)) Call mma_deallocate(idBM)
            Call mma_allocate(dBM,mdB_Tot,Label='dBM')
            Call mma_allocate(idBM,mdB_Tot*2,Label='idBM')
         Else
           If (.NOT.Allocated(dBM)) Call mma_allocate(dBM,1,Label='dBM')
           If (.NOT.Allocated(idBM)) Call mma_allocate(idBM,1*2,
     &                                                 Label='idBM')
         End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
         Call Get_Curvil
     &             (iq,iqRF,iqR,iqA,iqT,iqO,
     &              nAtoms,jIter,nIter,Cx,
     &              Proc,
     &              qVal,nq,qLbl,
     &              iRef, F_c,Mult,
     &              LuIC,Ind,iIter,
     &              GRef,iGlow,iGHi,
     &              Proc_dB,
     &              iTabBonds,iTabAtoms,nBonds,nMax,iTabAI,mAtoms,
     &              nB_Tot,ndB_Tot,
     &              BM,dBM,iBM,idBM,
     &              mB_Tot,mdB_Tot,nqBM,Thr_small)
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
         i = 1
         Do iq = 1, nq
            nB = nqBM(iq)
            Call DScal_(nB,Mult(iq),BM(i),1)
            i = i + nB
         End Do
         KtBu(1:nQQ*nDimBC)=Zero
         Do iQQ = 0, nQQ-1
            i = 1
            Do iq = 1, nq
               nB = nqBM(iq)
               Do iB = 0, nB-1
                  iDim = iBM(i)
                  KtBu((iDim-1)*nQQ + iQQ + 1) =
     &                 KtBu((iDim-1)*nQQ + iQQ + 1)
     &               + K(iQQ*nq+iq) * BM(i)
                  i = i + 1
               End Do
            End Do
         End Do
#ifdef _DEBUGPRINT_
         If (iPrint.ge.99) Then
            Call RecPrt(' The K matrix',' ',K,nq,nQQ)
            Call RecPrt(' The K(t)B matrix',' ',KtBu,nQQ,nDimBC)
         End If
#endif
*                                                                      *
************************************************************************
*                                                                      *
*------- Store the B matrix for the last structure
*
         If (jIter.eq.nIter) Then
            nX = 3*nAtoms
            Call mma_allocate(BMx,nX,nX,Label='BMx')
            BMx(:,:)=Zero
*
*           modify from compact to full Cartesian storage.
*
            iDim = 0
            Do iX = 1, nX
               iAtom = (iX+2)/3
               ixyz = iX - (iAtom-1)*3
               If (Smmtrc(ixyz,iAtom)) Then
                  iDim = iDim + 1
                  Do iQQ = 1, nQQ
                     iQD = (iDim-1)*nQQ + iQQ
                     BMx(iX,iQQ) = KtBu(iQD)
                  End Do
               Else
                  Do iQQ = 1, nQQ
                     BMx(iX,iQQ) = Zero
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
*           KtBu is stored nQQ, nDimBC
*
            Call mma_allocate(KtBt,nDimBC,nQQ,Label='KtBt')
            Call TRNSPS(nQQ,nDimBC,KtBu,KtBt)
*
*           Strip KtB of the degeneracy factor (full).
*
            Do iQQ = 1, nQQ
               Do iDim = 1, nDimBC
                  KtBt(iDim,iQQ) = KtBt(iDim,iQQ) / Degen2(iDim)
               End Do
            End Do
*
            M = nDimBC
            N = nQQ
            NRHS=1
            Call Eq_Solver('N',M,N,NRHS,KtBt,.False.,
     &                     Degen2,GxR(:,iOff),dqInt(:,jIter))
*           Call RecPrt('GxR(:,iSt',' ',GxR(:,iOff),nDimBC,1)
*           Call RecPrt('KtB   ',' ',KtBt,nQQ,nDimBC)
*           Call RecPrt('drInt',' ',dqInt(:,jIter),nQQ,1)

            iOff = iOff - 1
*                                                                      *
************************************************************************
*                                                                      *
*---------- Save pointer to KtB to be used in backtransforming the
*           Hessian to internal coordinates.
*
            If (Proc_H) Then
               Call mma_allocate(KtB,nDimBC,nQQ,Label='KtB')
               Call DCopy_(nDimBC*nQQ,KtBt,1,KtB,1)
*              Call RecPrt('KtB',' ',KtB,nDimBC,nQQ)
            End If
            Call mma_deallocate(KtBt)
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
      Call DGEMM_('N','N',
     &            nQQ,mIter,nq,
     &            1.0d0,KtM,nQQ,
     &                  qVal(:,jIter),nq,
     &            0.0d0,qInt(:,jIter),nQQ)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      If (iPrint.ge.49) Then
         Call RecPrt(' The K Matrix',' ',K,nq,nQQ)
         Call RecPrt(' q-values',' ',qVal,nq,nIter)
         Call RecPrt('Q-values',' ',qInt,nQQ,nIter)
         Call RecPrt('Cx',' ',Cx,3*nAtoms,nIter)
      End If
      If (BSet.and.iPrint.ge.49) Then
          Call RecPrt('Q-gradients',' ',dqInt,nQQ,nIter)
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
      Call mma_deallocate(KtBu)
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
