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
      Subroutine PLF_Cho(TInt,lInt,
     &                AOint,ijkl,iCmp,jCmp,kCmp,lCmp,iShell,
     &                iAO,iAOst,Shijij,iBas,jBas,kBas,lBas,kOp)
************************************************************************
*                                                                      *
*  object: to sift and index the petite list format integrals.         *
*                                                                      *
*          the indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
*                                                                      *
*  Author: Roland Lindh, IBM Almaden Research Center, San Jose, Ca     *
*          May '90                                                     *
*                                                                      *
************************************************************************
      use SOAO_Info, only: iAOtSO
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
      Real*8 AOint(ijkl,iCmp,jCmp,kCmp,lCmp), TInt(lInt)
      Integer iShell(4), iAO(4), kOp(4),
     &        iAOst(4), iSOs(4)
      Logical Shijij

      external ddot_

      INTEGER ABCD, CDAB, CD, AB, A, B, C, D
*
      iTri(i,j)=Max(i,j)*(Max(i,j)-3)/2 + i + j
      iSOShl(i)=iWork(ip_iSOShl-1+i)
      iShlSO(i)=iWork(ip_iShlSO-1+i)
      nBstSh(i)=iWork(ip_nBstSh-1+i)
*
#if defined (_DEBUG_)
      Call qEnter('Plf_Cho')
#endif
      irout = 109
      jprint = nprint(irout)
      If (jPrint.ge.49) Then
         r1=DDot_(ijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,[One],0)
         r2=DDot_(ijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,AOInt,1)
         Write (6,*) ' Sum=',r1
         Write (6,*) ' Dot=',r2
      End If
      If (jPrint.ge.99) Call RecPrt(' In Plf_CD: AOInt',' ',
     &                              AOInt,ijkl,iCmp*jCmp*kCmp*lCmp)

      If (Shijij) Then ! avoid compiler warnings about unused variables
         iDummy_1 = iShell(1)
      End If

      NUMC = NBSTSH(SHC)
      NUMD = NBSTSH(SHD)
      NUMA = NBSTSH(SHA)
      NUMB = NBSTSH(SHB)

      IF (SHC .EQ. SHD) THEN
         NUMCD = NUMC*(NUMD + 1)/2
      ELSE
         NUMCD = NUMC*NUMD
      END IF
      IF (SHA .EQ. SHB) THEN
         NUMAB = NUMA*(NUMB + 1)/2
      ELSE
         NUMAB = NUMA*NUMB
      END IF
      NCDAB = NUMCD*NUMAB
      IF (NCDAB .NE. LINT) THEN
         WRITE(LUPRI,*) 'PLF_Cho: dimension of integral array: ',NCDAB
         WRITE(LUPRI,*) 'PLF_Cho: expected                   : ',LINT
         WRITE(LUPRI,*) 'PLF_Cho: YOU HAVE A DIMENSION PROBLEM!'
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
*     Allocate space to store integrals to gether with their
*     Symmetry batch and sequence number.
*     To avoid conflicts in using memory this is done in the
*     subroutine PSOAO
*
*
*     quadruple loop over elements of the basis functions angular
*     description. loops are reduced to just produce unique SO integrals
*     observe that we will walk through the memory in AOint in a
*     sequential way.
*
      iAOsti=iAOst(1)
      iAOstj=iAOst(2)
      iAOstk=iAOst(3)
      iAOstl=iAOst(4)
      iAOi=iAO(1)
      iAOj=iAO(2)
      iAOk=iAO(3)
      iAOl=iAO(4)
*
      ijklCmp=iCmp*jCmp*kCmp*lCmp
*
      Do 100 i1 = 1, iCmp
         iSOs(1)=iAOtSO(iAOi+i1,kOp(1))+iAOsti
         Do 200 i2 = 1, jCmp
            iSOs(2)=iAOtSO(iAOj+i2,kOp(2))+iAOstj
            Do 300 i3 = 1, kCmp
               iSOs(3)=iAOtSO(iAOk+i3,kOp(3))+iAOstk
               Do 400 i4 = 1, lCmp
                  iSOs(4)=iAOtSO(iAOl+i4,kOp(4))+iAOstl
*
                iSO =iSOs(1)
                jSO =iSOs(2)
                kSO =iSOs(3)
                lSO =iSOs(4)
*
                nijkl = 0
                Do 120 lSOl = lSO, lSO+lBas-1
                   Do 220 kSOk = kSO, kSO+kBas-1
                      Do 320 jSOj = jSO, jSO+jBas-1
                         Do 420 iSOi = iSO, iSO+iBas-1
*
                            nijkl = nijkl + 1
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
                               CALL CHO_QUIT('Logical error in PLF_Cho',
     &                                       103)
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
                               CALL CHO_QUIT('PLF_Cho: Error!',104)
                            END IF
#endif
                            TINT(CDAB) = AOint(nijkl,i1,i2,i3,i4)
                            IF (ISHLCD .EQ. ISHLAB) THEN
                               IF (SHC.EQ.SHD .OR. SHC.EQ.SHA) THEN
                                  ABCD = NUMAB*(CD - 1) + AB
                                  TINT(ABCD) = AOint(nijkl,
     &                                               i1,i2,i3,i4)
                               ELSE IF (SHC.EQ.SHB) THEN
                                  CD   = NUMD*(C - 1) + D
                                  AB   = NUMB*(A - 1) + B
                                  ABCD = NUMAB*(CD - 1) + AB
                                  TINT(ABCD) = AOint(nijkl,
     &                                               i1,i2,i3,i4)
                               END IF
                            END IF

420                      Continue
320                   Continue
220                Continue
120             Continue

400            Continue
300         Continue
200      Continue
100   Continue

#if defined (_DEBUG_)
      Call qExit('Plf_Cho')
#endif
      Return
      End
