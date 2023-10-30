!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine IndSft2(iCmp,iShell,iBas,jBas,kBas,lBas,
     &                   Shijij, iAO, iAOst, ijkl,SOint,nSOint,
     &                   iSOSym,nSOs)
!***********************************************************************
!  object: to sift and index the SO integrals.                         *
!                                                                      *
!          the indices has been scrambled before calling this routine. *
!          Hence we must take special care in order to regain the can- *
!          onical order.                                               *
!                                                                      *
!  Author: Roland Lindh, IBM Almaden Research Center, San Jose, Ca     *
!          april '90                                                   *
!***********************************************************************
      use k2_arrays, only: Sew_Scr
      use SOAO_Info, only: iAOtSO, iOffSO
      use lw_Info, only: lwInt, lwSqn, lwSyb
      use Gateway_Info, only: ThrInt
      use Symmetry_Info, only: nIrrep
      use sort_data, only: DimSyB, iStBin, lSll, mxSyP, nSkip, Square,
     &                     TriSyB
#ifdef _DEBUGPRINT_
      use Constants, only: Zero, One
#endif
      Implicit None
!
      Integer ijkl, nSOInt, ibas, jBas, kBas, lBas, nSOs
      Real*8 SOint(ijkl,nSOint)
      Integer iCmp(4), iShell(4), iAO(4), iAOst(4), iSOSym(2,nSOs)
      Logical Shijij

!     local array
      Logical Shij, Shkl, qijij, qij, qkl
      Integer iSym(0:7), jSym(0:7), kSym(0:7), lSym(0:7)
      Integer k12, k34, MemSO2, nUt, i1, i2, i3, i4, j1, j2, j3, j4,
     &        jCmpMx, lCmpMx, iSymi, jSymj, kSymk, lSyml,
     &        iSO, jSO, kSO, lSO, i12, i34, iSq1, iSq2, iSq3, iSq4,
     &        iqq1, iqq2, iqq3, iqq4, iSym12, iSym34, iSyBlk, jSyBlk,
     &        nij, nkl, ipP1, ipP2, ipP3, ipP4, iSOi, jSOj, kSOk, lSOl,
     &        ij, kl, iSqNum, jSqNum, j, ix, j2max, j12, nijkl, ipD,
     &        iBin, jBin
      Real*8 AInt
      Logical dupli
#ifdef _DEBUGPRINT_
      Real*8, Save :: Tr1=Zero, Tr2=Zero
      Real*8, External:: DDot_
      Real*8 r1, r2
#endif

      k12=0
      k34=0
#ifdef _DEBUGPRINT_
      r1=DDot_(ijkl*nSOInt,SOInt,1,[One],0)
      r2=DDot_(ijkl*nSOInt,SOInt,1,SOInt,1)
      tr1=tr1+r1
      tr2=tr2+r2
      Write (6,*) ' Sum=',r1,tr1
      Write (6,*) ' Dot=',r2,tr2
      Call RecPrt(' in indsft:SOint ',' ',SOint,ijkl,nSOint)
#endif
      memSO2 = 0
!
!     allocate space to store integrals together with their
!     Symmetry batch and sequence number
!     To avoid conflicts in using memory this is done in the
!     subroutine PSOAO
!
      nUt=-1
!
!     quadruple loop over elements of the basis functions angular
!     description. loops are reduced to just produce unique SO integrals
!     observe that we will walk through the memory in AOint in a
!     sequential way.
!
      Shij = iShell(1).eq.iShell(2)
      Shkl = iShell(3).eq.iShell(4)
      Do 100 i1 = 1, iCmp(1)
         Do 101 j = 0, nIrrep-1
            ix = 0
            If (iAOtSO(iAO(1)+i1,j)>0) ix = 2**j
            iSym(j) = ix
101      Continue
         jCmpMx = iCmp(2)
         If (Shij) jCmpMx = i1
         Do 200 i2 = 1, jCmpMx
            Do 201 j = 0, nIrrep-1
               ix = 0
               If (iAOtSO(iAO(2)+i2,j)>0) ix = 2**j
               jSym(j) = ix
201         Continue
            qij = i1.eq.i2
            If (iShell(2).gt.iShell(1)) then
               i12 = iCmp(2)*(i1-1) + i2
            else
               i12 = iCmp(1)*(i2-1) + i1
            End If
            Do 300 i3 = 1, iCmp(3)
               Do 301 j = 0, nIrrep-1
                  ix = 0
                  If (iAOtSO(iAO(3)+i3,j)>0) ix = 2**j
                  kSym(j) = ix
301            Continue
               lCmpMx = iCmp(4)
               If (Shkl) lCmpMx = i3
               Do 400 i4 = 1, lCmpMx
                  Do 401 j = 0, nIrrep-1
                     ix = 0
                     If (iAOtSO(iAO(4)+i4,j)>0) ix = 2**j
                     lSym(j) = ix
401               Continue
                  qkl = i3.eq.i4
                  If (iShell(4).gt.iShell(3)) then
                     i34 = iCmp(4)*(i3-1) + i4
                  else
                     i34 = iCmp(3)*(i4-1) + i3
                  End If
                  If (Shijij .and. i34.gt.i12) go to 400
                  qijij = Shijij .and. i12.eq.i34
!          Write (6,*) 'i1,i2,i3,i4=',i1,i2,i3,i4
!
!      loop over Irreps which are spanned by the basis function.
!      again, the loop structure is restricted to ensure unique
!      integrals.
!
       Do 110 j1 = 0, nIrrep-1
          If (iSym(j1).eq.0) go to 110
          j2max = nIrrep-1
          If (Shij .and. qij) j2max = j1
          Do 210 j2 = 0, j2max
             If (jSym(j2).eq.0) go to 210
             j12 = ieor(j1,j2)
             If (qijij) then
                If (Shij .and. qij) then
                    k12 = j1*(j1+1)/2 + j2+1
                else If (Shij) then
                    k12 = nIrrep*j1 + j2+1
                else If (iShell(1).gt.iShell(2)) then
                    k12 = nIrrep*j1 + j2+1
                else
                    k12 = nIrrep*j2 + j1+1
                End If
             End If
!
             iSymi=max(j1,j2)+1
             jSymj=min(j1,j2)+1
!
             Do 310 j3 = 0, nIrrep-1
                If (kSym(j3).eq.0) go to 310
                j4 = ieor(j12,j3)
                If (lSym(j4).eq.0) go to 310
                If (Shkl .and. qkl .and. j4.gt.j3) go to 310
                If (qijij) then
                   If (Shkl .and. qkl) then
                      k34 = j3*(j3+1)/2 + j4+1
                   else If (Shkl) then
                      k34 = nIrrep*j3 + j4+1
                   else If (iShell(3).gt.iShell(4)) then
                      k34 = nIrrep*j3 + j4+1
                   else
                      k34 = nIrrep*j4 + j3+1
                   End If
                   If (k34.gt.k12) go to 310
                End If
!               Write (6,*) 'j1,j2,j3,j4=',j1,j2,j3,j4
!
                memSO2 = memSO2 + 1
                If ( (nSkip(j1+1)+nSkip(j2+1)+
     &                nSkip(j3+1)+nSkip(j4+1) ).ne.0 ) GoTo 310
!
!               Compute absolute starting SO index
                iSO = iAOtSO(iAO(1)+i1,j1)+iAOst(1)+iOffSO(j1)
                jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)+iOffSO(j2)
                kSO = iAOtSO(iAO(3)+i3,j3)+iAOst(3)+iOffSO(j3)
                lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)+iOffSO(j4)
!               Write (6,*) 'iSO,jSO,kSO,lSO=',iSO,jSO,kSO,lSO
!
                kSymk=max(j3,j4)+1
                lSyml=min(j3,j4)+1
!
                iSym12=TriSyB(iSymi,jSymj)
                iSym34=TriSyB(kSymk,lSyml)
!
!               Order Irrep index canonical
                If ( iSym34.gt.iSym12 ) then
                   iSyBlk= mxSyP*(iSym34-1)+iSym12
                   jSyBlk= mxSyP*(iSym12-1)+iSym34
                   nij=DimSyB(kSymk,lSyml)
                   nkl=DimSyB(iSymi,jSymj)
                   iSq1=nkl
                   iSq2=1
                   iSq3=1
                   iSq4=nij
                   iQQ1=lSll(iSyBlk)/nkl
                   iQQ2=nkl+1
                   iQQ3=nij+1
                   iQQ4=lSll(jSyBlk)/nij
                   iPP1=iQQ1
                   iPP2=0
                   iPP3=0
                   iPP4=iQQ4
                Else
                   iSyBlk= mxSyP*(iSym12-1)+iSym34
                   jSyBlk= mxSyP*(iSym34-1)+iSym12
                   nij=DimSyB(iSymi,jSymj)
                   nkl=DimSyB(kSymk,lSyml)
                   iSq1=1
                   iSq2=nkl
                   iSq3=nij
                   iSq4=1
                   iQQ1=nkl+1
                   iQQ2=lSll(iSyBlk)/nkl
                   iQQ3=lSll(jSyBlk)/nij
                   iQQ4=nij+1
                   iPP1=0
                   iPP2=iQQ2
                   iPP3=iQQ3
                   iPP4=0
                End If
!               Write (*,*) '<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
!               Write (*,*) 'iSyBlk,jSyBlk=',iSyBlk,jSyBlk
!               Write (*,*) 'iSq1,iSq2,iSq3,iSq4=',
!    &                       iSq1,iSq2,iSq3,iSq4
!               Write (*,*) 'iQQ1,iQQ2,iQQ3,iQQ4=',
!    &                       iQQ1,iQQ2,iQQ3,iQQ4
!               Write (*,*) 'nkl,nij=',
!    &                       nkl,nij
!
!               Duplicate integral if permuted integral is
!               in the same irrep or if all blocks should
!               be duplicated.
                dupli=(iSyBlk.eq.jSyBlk).or.Square
!
                If (.Not.dupli) Then
!
                nijkl = 0
                Do lSOl = lSO, lSO+lBas-1
                   Do kSOk = kSO, kSO+kBas-1
                      kl=iPD(kSOk,lSOl,iSOSym,nSOs)
                      Do jSOj = jSO, jSO+jBas-1
                         Do iSOi = iSO, iSO+iBas-1
                            nijkl = nijkl + 1
                            AInt=SOint(nijkl,memSO2)
                            If (Abs(AInt).lt.ThrInt) Go To 199
                            ij=iPD(iSOi,jSOj,iSOSym,nSOs)
!                           Write (*,*)
!                           Write (*,*) 'iSOi,jSOj,kSOk,lSOl=',
!    &                                   iSOi,jSOj,kSOk,lSO
!                           Write (*,*) 'ij,kl=',ij,kl
!
                            nUt=nUt+1
                            Sew_Scr(lwInt+nUt)=AInt
                            iBin=(kl-1)/iQQ1 +(ij-1)/iQQ2
                            iSqNum = (kl-iBin*iPP1)*iSq1
     &                             + (ij-iBin*iPP2)*iSq2
     &                             - nkl
                            Sew_Scr(lwSqN+nUt)=DBLE(iSqNum)
                            Sew_Scr(lwSyB+nUt)=DBLE(iBin+iStBin(iSyBlk))
!                           Write (*,*) 'iSqNum,iBin=',iSqNum,iBin+
!    &                                   iStBin(iSyBlk)
!
 199                        Continue
                         End Do
                      End Do
                   End Do
                End Do
!
                Else
!
                nijkl = 0
                Do lSOl = lSO, lSO+lBas-1
                   Do kSOk = kSO, kSO+kBas-1
                      kl=iPD(kSOk,lSOl,iSOSym,nSOs)
                      Do jSOj = jSO, jSO+jBas-1
                         Do iSOi = iSO, iSO+iBas-1
                            nijkl = nijkl + 1
                            AInt=SOint(nijkl,memSO2)
                            If (Abs(AInt).lt.ThrInt) Go To 299
                            ij=iPD(iSOi,jSOj,iSOSym,nSOs)
!
!                           Write (*,*)
!                           Write (*,*) 'iSOi,jSOj,kSOk,lSOl=',
!    &                                   iSOi,jSOj,kSOk,lSO
!                           Write (*,*) 'ij,kl=',ij,kl
!
                            nUt=nUt+1
                            Sew_Scr(lwInt+nUt)=AInt
                            iBin=(kl-1)/iQQ1 +(ij-1)/iQQ2
                            iSqNum = (kl-iBin*iPP1)*iSq1
     &                             + (ij-iBin*iPP2)*iSq2
     &                             - nkl
                            Sew_Scr(lwSqN+nUt)=DBLE(iSqNum)
                            Sew_Scr(lwSyB+nUt)=DBLE(iBin+iStBin(iSyBlk))
!                           Write (*,*) 'iSqNum,iBin=',iSqNum,iBin+
!    &                                   iStBin(iSyBlk)
!
                            nUt=nUt+1
                            Sew_Scr(lwInt+nUt)=AInt
                            jBin=(kl-1)/iQQ3 +(ij-1)/iQQ4
                            jSqNum = (kl-jBin*iPP3)*iSq3
     &                             + (ij-jBin*iPP4)*iSq4
     &                             - nij
                            Sew_Scr(lwSqN+nUt)=DBLE(jSqNum)
                            Sew_Scr(lwSyB+nUt)=DBLE(jBin+iStBin(jSyBlk))
!                           Write (*,*) 'jSqNum,jBin=',jSqNum,jBin+
!    &                                   iStBin(jSyBlk)
!
 299                        Continue
                         End Do
                      End Do
                   End Do
                End Do
!
                End If
!
310          Continue
210       Continue
110    Continue
!
400            Continue
300         Continue
200      Continue
100   Continue
!
!     pass the integral to phase 1 of the bin sorting algorithm
!
      Call SORT1A(nUt+1,Sew_Scr(lwInt),Sew_Scr(lwSqN),Sew_Scr(lwSyB))
      nUt=0
      Return
      End SubRoutine IndSft2
