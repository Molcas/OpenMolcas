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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      SubRoutine IndSft_Cho(TInt,lInt,
     &                   iCmp,iShell,iBas,jBas,kBas,lBas,
     &                   Shijij, iAO, iAOst, ijkl,SOint,nSOint,
     &                   iSOSym,nSOs)
************************************************************************
*  object: to sift and index the SO integrals.                         *
*                                                                      *
*          the indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
*  Author: Roland Lindh, IBM Almaden Research Center, San Jose, Ca     *
*          april '90                                                   *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "cholesky.fh"
#include "choptr.fh"
#include "real.fh"
#include "print.fh"
#include "srt0.fh"
#include "WrkSpc.fh"
*
      Real*8 SOint(ijkl,nSOint), TInt(lInt)
      Integer iCmp(4), iShell(4), iAO(4), iAOst(4), iSOSym(2,nSOs)
      Logical Shijij, Shij, Shkl, qijij, qij, qkl
      INTEGER A, B, C, D, AB, CD, CDAB, ABCD
*     local array
      Integer iSym(0:7), jSym(0:7), kSym(0:7), lSym(0:7)
      Data tr1,tr2/0.0d0,0.0d0/
      Save tr1,tr2

      external ddot_
*
*     statement function
*
      iTri(i,j)=Max(i,j)*(Max(i,j)-3)/2 + i + j
      ISOSHL(I)=IWORK(ip_iSOShl-1+I)
      ISHLSO(I)=IWORK(ip_iShlSO-1+I)
      NBSTSH(I)=IWORK(ip_NBSTSH-1+I)
*
#if defined (_DEBUG_)
      Call qEnter('IndSft_Cho')
#endif
      irout = 39
      jprint = nprint(irout)
      k12=0
      k34=0
      If (jPrint.ge.49) Then
         r1=DDot_(ijkl*nSOInt,SOInt,1,One,0)
         r2=DDot_(ijkl*nSOInt,SOInt,1,SOInt,1)
         tr1=tr1+r1
         tr2=tr2+r2
         Write (6,*) ' Sum=',r1,tr1
         Write (6,*) ' Dot=',r2,tr2
      End If
      If (jprint.ge.99)
     &   Call RecPrt(' in indsft:SOint ',' ',SOint,ijkl,nSOint)
      memSO2 = 0

      NUMC = NBSTSH(SHC)
      NUMD = NBSTSH(SHD)
      NUMA = NBSTSH(SHA)
      NUMB = NBSTSH(SHB)

      If (nSOs .gt. 0) Then ! to make some compilers happy
         iDummy_2 = iSOSym(1,1)
      End If

      IF (SHC .EQ. SHD) THEN
         NUMCD = NUMC*(NUMC + 1)/2
      ELSE
         NUMCD = NUMC*NUMD
      END IF
      IF (SHA .EQ. SHB) THEN
         NUMAB = NUMA*(NUMA + 1)/2
      ELSE
         NUMAB = NUMA*NUMB
      END IF
      NCDAB = NUMCD*NUMAB
      IF (NCDAB .NE. LINT) THEN
         WRITE(LUPRI,*)
     &   'IndSft_Cho: dimension of integral array: ',NCDAB
         WRITE(LUPRI,*)
     &   'IndSft_Cho: expected                   : ',LINT
         WRITE(LUPRI,*)
     &   'IndSft_Cho: YOU HAVE A DIMENSION PROBLEM!'
      END IF

      ISHLCD = ITRI(SHC,SHD)
      ISHLAB = ITRI(SHA,SHB)

C to avoid stupid compiler warnings:

      C = 0
      D = 0
      A = 0
      B = 0

      NTELM = 0
*
*     allocate space to store integrals to gether with their
*     Symmetry batch and sequence number
*     To avoid conflicts in using memory this is done in the
*     subroutine PSOAO
*
*
*     quadruple loop over elements of the basis functions angular
*     description. loops are reduced to just produce unique SO integrals
*     observe that we will walk through the memory in AOint in a
*     sequential way.
*
      Shij = iShell(1).eq.iShell(2)
      Shkl = iShell(3).eq.iShell(4)
      Do 100 i1 = 1, iCmp(1)
         Do 101 j = 0, nIrrep-1
            iSym(j) = iand(IrrCmp(inds(iShell(1))+i1),2**j)
101      Continue
         jCmpMx = iCmp(2)
         If (Shij) jCmpMx = i1
         Do 200 i2 = 1, jCmpMx
            Do 201 j = 0, nIrrep-1
               jSym(j) = iand(IrrCmp(inds(iShell(2))+i2),2**j)
201         Continue
            qij = i1.eq.i2
            If (iShell(2).gt.iShell(1)) then
               i12 = iCmp(2)*(i1-1) + i2
            else
               i12 = iCmp(1)*(i2-1) + i1
            End If
            Do 300 i3 = 1, iCmp(3)
               Do 301 j = 0, nIrrep-1
                  kSym(j) = iand(IrrCmp(inds(iShell(3))+i3),2**j)
301            Continue
               lCmpMx = iCmp(4)
               If (Shkl) lCmpMx = i3
               Do 400 i4 = 1, lCmpMx
                  Do 401 j = 0, nIrrep-1
                     lSym(j) = iand(IrrCmp(inds(iShell(4))+i4),2**j)
401               Continue
                  qkl = i3.eq.i4
                  If (iShell(4).gt.iShell(3)) then
                     i34 = iCmp(4)*(i3-1) + i4
                  else
                     i34 = iCmp(3)*(i4-1) + i3
                  End If
                  If (Shijij .and. i34.gt.i12) go to 400
                  qijij = Shijij .and. i12.eq.i34
*
*      loop over Irreps which are spanned by the basis function.
*      again, the loop structure is restricted to ensure unique
*      integrals.
*
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
*
             iSymi=max(j1,j2)+1
             jSymj=min(j1,j2)+1
*
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
*
                memSO2 = memSO2 + 1
                If ( (nSkip(j1+1)+nSkip(j2+1)+
     &                nSkip(j3+1)+nSkip(j4+1) ).ne.0 ) GoTo 310
*
*               Compute absolute starting SO index
                iSO = iAOtSO(iAO(1)+i1,j1)+iAOst(1)+iOffSO(j1)
                jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)+iOffSO(j2)
                kSO = iAOtSO(iAO(3)+i3,j3)+iAOst(3)+iOffSO(j3)
                lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)+iOffSO(j4)
*
                kSymk=max(j3,j4)+1
                lSyml=min(j3,j4)+1
*
                nijkl = 0
                Do lSOl = lSO, lSO+lBas-1
                   Do kSOk = kSO, kSO+kBas-1
                      Do jSOj = jSO, jSO+jBas-1
                         Do iSOi = iSO, iSO+iBas-1
                            nijkl = nijkl + 1
*
                            NTELM = NTELM + 1
*
                            ISHLI = ISOSHL(ISOI)
                            ISHLJ = ISOSHL(JSOJ)
                            ISHLK = ISOSHL(KSOK)
                            ISHLL = ISOSHL(LSOL)

                            IF ((ISHLI.EQ.SHC).AND.(ISHLJ.EQ.SHD)
     &                      .AND.(ISHLK.EQ.SHA).AND.(ISHLL.EQ.SHB)) THEN
                               C = ISHLSO(ISOI)
                               D = ISHLSO(JSOJ)
                               A = ISHLSO(KSOK)
                               B = ISHLSO(LSOL)
                            ELSE IF ((ISHLJ.EQ.SHC).AND.(ISHLI.EQ.SHD)
     &                      .AND.(ISHLK.EQ.SHA).AND.(ISHLL.EQ.SHB)) THEN
                               C = ISHLSO(JSOJ)
                               D = ISHLSO(ISOI)
                               A = ISHLSO(KSOK)
                               B = ISHLSO(LSOL)
                            ELSE IF ((ISHLI.EQ.SHC).AND.(ISHLJ.EQ.SHD)
     &                      .AND.(ISHLL.EQ.SHA).AND.(ISHLK.EQ.SHB)) THEN
                               C = ISHLSO(ISOI)
                               D = ISHLSO(JSOJ)
                               A = ISHLSO(LSOL)
                               B = ISHLSO(KSOK)
                            ELSE IF ((ISHLJ.EQ.SHC).AND.(ISHLI.EQ.SHD)
     &                      .AND.(ISHLL.EQ.SHA).AND.(ISHLK.EQ.SHB)) THEN
                               C = ISHLSO(JSOJ)
                               D = ISHLSO(ISOI)
                               A = ISHLSO(LSOL)
                               B = ISHLSO(KSOK)
                            ELSE IF ((ISHLK.EQ.SHC).AND.(ISHLL.EQ.SHD)
     &                      .AND.(ISHLI.EQ.SHA).AND.(ISHLJ.EQ.SHB)) THEN
                               C = ISHLSO(KSOK)
                               D = ISHLSO(LSOL)
                               A = ISHLSO(ISOI)
                               B = ISHLSO(JSOJ)
                            ELSE IF ((ISHLL.EQ.SHC).AND.(ISHLK.EQ.SHD)
     &                      .AND.(ISHLI.EQ.SHA).AND.(ISHLJ.EQ.SHB)) THEN
                               C = ISHLSO(LSOL)
                               D = ISHLSO(KSOK)
                               A = ISHLSO(ISOI)
                               B = ISHLSO(JSOJ)
                            ELSE IF ((ISHLK.EQ.SHC).AND.(ISHLL.EQ.SHD)
     &                      .AND.(ISHLJ.EQ.SHA).AND.(ISHLI.EQ.SHB)) THEN
                               C = ISHLSO(KSOK)
                               D = ISHLSO(LSOL)
                               A = ISHLSO(JSOJ)
                               B = ISHLSO(ISOI)
                            ELSE IF ((ISHLL.EQ.SHC).AND.(ISHLK.EQ.SHD)
     &                      .AND.(ISHLJ.EQ.SHA).AND.(ISHLI.EQ.SHB)) THEN
                               C = ISHLSO(LSOL)
                               D = ISHLSO(KSOK)
                               A = ISHLSO(JSOJ)
                               B = ISHLSO(ISOI)
                            ELSE
                               WRITE(LUPRI,*)
     &                         'Shell quadruple requested: ',
     &                         SHC,SHD,SHA,SHB
                               WRITE(LUPRI,*)
     &                         'Shell quadruple of element ',NTELM,':',
     &                         ISHLI,ISHLJ,ISHLK,ISHLL
                               CALL CHO_QUIT(
     &                                    'Logical error in IndSft_Cho',
     &                                    103)
                            END IF

                            IF (SHA .EQ. SHB) THEN
                               AB = ITRI(A,B)
                            ELSE
                               AB = NUMA*(B - 1) + A
                            END IF
                            IF (SHC .EQ. SHD) THEN
                               CD = ITRI(C,D)
                            ELSE
                               CD = NUMC*(D - 1) + C
                            END IF

                            CDAB = NUMCD*(AB - 1) + CD
#if defined (_DEBUG_)
                            IF ((CDAB.GT.LINT) .OR. (CDAB.LT.1)) THEN
                               WRITE(LUPRI,*) 'CDAB: ',CDAB
                               WRITE(LUPRI,*) 'Dimension: ',LINT
                               CALL CHO_QUIT('IndSft_Cho: Error!',104)
                            END IF
#endif
                            TINT(CDAB) = SOint(nijkl,memSO2)
                            IF (ISHLCD .EQ. ISHLAB) THEN
                               IF (SHC.EQ.SHD .OR. SHC.EQ.SHA) THEN
                                  ABCD = NUMAB*(CD - 1) + AB
                                  TINT(ABCD) = SOint(nijkl,memSO2)
                               ELSE IF (SHC.EQ.SHB) THEN
                                  CD   = NUMD*(C - 1) + D
                                  AB   = NUMB*(A - 1) + B
                                  ABCD = NUMAB*(CD - 1) + AB
                                  TINT(ABCD) = SOint(nijkl,memSO2)
                               END IF
                            END IF

                         End Do
                      End Do
                   End Do
                End Do
*
310          Continue
210       Continue
110    Continue
*
400            Continue
300         Continue
200      Continue
100   Continue
*
#if defined (_DEBUG_)
      Call qExit('IndSft_Cho')
#endif
      Return
      End
