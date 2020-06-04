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
      Subroutine CurviL(nAtoms,nDim,Cx,Gx,nIter,iIter,iRef,nStab,iOper,
     &                  nSym,jStab,Degen,Smmtrc,mTR,TRVec,
     &                  ip_rInt,ip_drInt,HSet,BSet,ipBMx,Numerical,iANr,
     &                  HWRS,Analytic_Hessian,iOptC,Name,PrQ,Proj,
     &                  dMass,iCoSet,iTabBonds,
     &                  iTabAtoms, nBonds,nMax,iTabAI,mAtoms,lOld,
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
#include "db.fh"
#include "print.fh"
      Real*8 Cx(3*nAtoms,nIter), Degen(3*nAtoms),
     &       Gx(3*nAtoms,nIter), Proj(nDim), dMass(nAtoms),
     &       TRVec(nDim,mTR)
      Integer nStab(nAtoms), iOper(0:nSym-1), iCoSet(0:7,nAtoms),
     &        jStab(0:7,nAtoms), iANr(nAtoms), iDum(6),
     &        iTabBonds(3,nBonds), iTabAtoms(0:nMax,nAtoms),
     &        iTabAI(2,mAtoms)
      Logical Smmtrc(3*nAtoms), HSet, BSet, Proc, Numerical, g12K,
     &        HWRS, Analytic_Hessian, PrQ,
     &        Proc_dB, lOld, Proc_H
      Character filnam*32
      Character*(LENIN) Name(nAtoms)
      Character*14 cDum
      Dimension Dum(1)
      Save g12K
      Data g12K/.False./
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
      iPrint=99
#else
      iRout=128
      iPrint=nPrint(iRout)+1
#endif
      Call QEnter('Curvil')
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
      Call Allocate_Work(ipDegen,nDim)
      i=0
      Do ix = 1, 3*nAtoms
         If (Smmtrc(ix)) Then
            Work(ipDegen+i) = Degen(ix)
            i = i + 1
         End If
      End Do
#ifdef _DEBUG_
      Call RecPrt('Work(ipDegen)',' ',Work(ipDegen),nDim,1)
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
     &           nAtoms,iIter,nIter,Cx,iOper,nSym,jStab,
     &           nStab,nDim,Smmtrc,Proc,Dum,1,iANr,cDum,
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
#ifdef _DEBUG_
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
      Call GetMem('q-Values','Allo','Real',ipqVal,nq*nIter)
      Call GetMem('q-Labels','Allo','Char',ipqLbl,14*nq)
      Call GetMem('F_const','Allo','Real',ipf_c,nq)
      Call GetMem('Mult','Allo','Real',ipMult,nq)
      Call GetMem('Ind','Allo','Inte',ipInd,3*nq)
      Call GetMem('GRef','Allo','Real',ipGRef,9*nqA*nIter)
*
      Call FZero(Work(ipf_c),nq)
      Call FZero(Work(ipMult),nq)
      Call FZero(Work(ipqVal),nq*nIter)
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
     &           nAtoms,iIter,nIter,Cx,iOper,nSym,jStab,
     &           nStab,nDim,Smmtrc,Proc,
     &           Work(ipqVal),nq,iANr,cWork(ipqLbl),
     &           iRef,Work(ipf_c),Work(ipMult),iOptC,
     &           LuIC,Name,iWork(ipInd),iIter,dMass,iCoSet,Work(ipGRef),
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
#ifdef _DEBUG_
      If (iPrint.ge.49) Then
         Write (6,*) 'nq, nqB, nqA, nqT, nqO=',
     &             nq, nqB, nqA, nqT, nqO
         Call RecPrt('q-values',' ',Work(ipqVal),nq,nIter)
         Call RecPrt('Force Constant matrix in redundant basis',' ',
     &            Work(ipf_c),1,nq)
         Call RecPrt('Multiplicity factors',' ',Work(ipMult),1,nq)
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
      Call GetMem('EVec','Allo','Real',ipK,nq**2)
      Call GetMem('K(t)B','Allo','Real',ipKtB,nDim**2)
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
            Call DScal_(nB,Work(ipf_c+iq),Work(ip_B+i),1)
            i = i + nB
         End Do
*ifdef _DEBUG_
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
         Call DScal_(nB,Work(ipMult+iq),Work(ip_B+i),1)
         i = i + nB
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      If (BSet) Then
         Call GetMem('G','Allo','Real',ipG,nq*nq)
         Call GetMem('EVal','Allo','Real',ipEVal,nq*(nq+1)/2)
         Call ElRed2(nq,ndim,Work(ipG),Work(ipEVal),Work(ipK),nK,Proj,
     &               g12K,Thr_ElRed,Work(ip_B),iWork(ip_iB),mB_Tot,
     &               iWork(ip_nqB))

         If (nK.gt.nQQ) Then
            Call Remove_TR(nq,ndim,nQQ,Work(ipK),nK,TRVec,mTR,
     &                     Work(ip_B),iWork(ip_iB),iWork(ip_nqB),mB_Tot)
         End If
         Call GetMem('EVal','Free','Real',ipEVal,nq*(nq+1)/2)
         Call GetMem('G','Free','Real',ipG,nq*nq)
*
         Call Put_dArray('K',Work(ipK),nq*nQQ)
      Else
         Call Get_dArray('K',Work(ipK),nq*nQQ)
      End If
#ifdef _DEBUG_
      If (iPrint.ge.99) Then
         Call RecPrt('K',' ',Work(ipK),nq,nQQ)
      End If
#endif
*                                                                      *
************************************************************************
*                                                                      *
*---- Print the nonreduntant coordinates.
*
      If (PrQ) Then
         Call GetMem('Temp2','Allo','Real',ipTemp2,nq*nQQ)
         call dcopy_(nQQ*nq,Work(ipK),1,Work(ipTemp2),1)
         Do iq = 0, nq-1
            Call DScal_(nQQ,One/Work(ipMult+iq),Work(ipTemp2+iq),nq)
         End Do
         Call PrintQ(Work(ipTemp2),cWork(ipqLbl),nq,nQQ,LuIC,
     &               Work(ipMult))
         Call GetMem('Temp2','Free','Real',ipTemp2,nq*nQQ)
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
         Call GetMem('GxR','Allo','Real',ipGxR,nDim*nIter)
         iOff = ipGxR
         Do jIter = 1, nIter
            Call NRed(Gx(1,jIter),Work(iOff),3*nAtoms,nDim,Smmtrc)
            iOff = iOff + nDim
         End Do
*
         iSt = nIter
         iEnd = iSt - Min(nIter,nWndw+1) + 1
C        iEnd = 1
         iOff = (nIter-1)*nDim + ipGxR
      Else
         If (Numerical) Then
            Call GetMem('GxR','Allo','Real',ipGxR,nDim)
            iOff = ipGxR
            Call NRed(Gx(1,nIter),Work(iOff),3*nAtoms,nDim,Smmtrc)
         Else
            iOff = ip_Dummy
         End If
         iEnd=nIter
         iSt =nIter
      End If
*
*     Note that the loop is in reverse order.
*
      Do jIter = iSt, iEnd, -1
         Proc_H=HSet.and.jIter.eq.iRef.and..Not.lOld
         Proc_dB=Proc_H.and.Analytic_Hessian
*        Compute and store dBQQ in the reference structure
         Proc_dB=Proc_dB.or.(Proc_H.and.Numerical)
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
     &              nAtoms,jIter,nIter,Cx,iOper,nSym,jStab,
     &              nStab,nDim,Smmtrc,Proc,
     &              Work(ipqVal),nq,iANr,cWork(ipqLbl),
     &              iRef, Work(ipf_c),Work(ipMult),
     &              iOptC,LuIC,Name,iWork(ipInd),iIter,dMass,iCoSet,
     &              Work(ipGRef),iGlow,iGHi,
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
            Call DScal_(nB,Work(ipMult+iq),Work(ip_B+i),1)
            i = i + nB
         End Do
         Call FZero(Work(ipKtB),nQQ*nDim)
         Do iQQ = 0, nQQ-1
            i = 0
            Do iq = 0, nq-1
               nB = iWork(ip_nqB+iq)
               Do iB = 0, nB-1
                  iDim = iWork(ip_iB+i)
                  Work((iDim-1)*nQQ + iQQ + ipKtB) =
     &                 Work((iDim-1)*nQQ + iQQ + ipKtB)
     &               + Work(iQQ*nq+iq+ipK)
     &               * Work(i+ip_B)
                  i = i + 1
               End Do
            End Do
         End Do
#ifdef _DEBUG_
         If (iPrint.ge.99) Then
            Call RecPrt(' The K matrix',' ',Work(ipK),nq,nQQ)
            Call RecPrt(' The K(t)B matrix',' ',Work(ipKtB),nQQ,nDim)
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
                     iQD = (iDim-1)*nQQ + iQQ + ipKtB -1
                     Work(iXQ) = Work(iQD)
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
            Call TRNSPS(nQQ,nDim,Work(ipKtB   ),Work(ip_KtB))
*
*           Strip KtB of the degeneracy factor (full).
*
            Do iQQ = 1, nQQ
               Do iDim = 1, nDim
                  ij = (iQQ-1)*nDim + iDim -1 + ip_KtB
                  Work(ij) = Work(ij) / Work(ipDegen+iDim-1)
               End Do
            End Do
*
            ip = ip_drInt + (jIter-1)*nQQ
            M = nDim
            N = nQQ
            NRHS=1
            Call Eq_Solver('N',M,N,NRHS,Work(ip_KtB   ),.False.,
     &                     Work(ipDegen),Work(iOff),Work(ip))
*           Call RecPrt('Work(iOff)',' ',Work(iOff),nDim,1)
*           Call RecPrt('KtB   ',' ',Work(ipKtB   ),nQQ,nDim)
*           Call RecPrt('drInt',' ',Work(ip),nQQ,1)

            iOff = iOff - nDim
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
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*---- Process the values
*          t
*     Q = K M q
*
      ipqV=ipqVal
      jIter = 1
      mIter = nIter
      If (.Not.BSet) Then
         ipqV = (nIter-1)*nq + ipqVal
         jIter = nIter
         mIter = 1
      End If
      Call GetMem('KtM','Allo','Real',ipKtM,nQQ*nq)
      Do iq = 0, nq-1
         Alpha=Work(ipMult+iq)
         Do iQQ = 0, nQQ-1
            temp = Work(ipK+iQQ*nq+iq)*Alpha
            Work(ipKtM+iq*nQQ+iQQ)=temp
         End Do
      End Do
      ip = ip_rInt + (jIter-1)*(nDim-mTR)
      Call DGEMM_('N','N',
     &            nQQ,mIter,nq,
     &            1.0d0,Work(ipKtM),nQQ,
     &                  Work(ipqV),nq,
     &            0.0d0,Work(ip),nQQ)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
      If (iPrint.ge.49) Then
         Call RecPrt(' The K Matrix',' ',Work(ipK),nq,nQQ)
         Call RecPrt(' q-values',' ',Work(ipqVal),nq,nIter)
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
      Call GetMem('KtM','Free','Real',ipKtM,nQQ*nq)
      If (BSet .or. Numerical)
     &   Call GetMem('GxR','Free','Real',ipGxR,nDim*nIter)
      Call GetMem('K(t)B','Free','Real',ipKtB,nQQ*nDim)
      Call GetMem('EVec','Free','Real',ipK,nq**2)
      Call GetMem('GRef','Free','Real',ipGRef,9*nqA*nIter)
      Call GetMem('Ind','Free','Inte',ipInd,3*nq)
      Call GetMem('Mult','Free','Real',ipMult,nq)
      Call GetMem('F_const','Free','Real',ipf_c,nq)
      Call GetMem('q-Labels','Free','Char',ipqLbl,14*nq)
      Call GetMem('Q-Values','Free','Real',ipQVal,nq*nIter)
*
      Call Free_Work(ipDegen)
*
      Close (LuIC)
*                                                                      *
************************************************************************
*                                                                      *
      Call QExit ('Curvil')
      Return
      End
