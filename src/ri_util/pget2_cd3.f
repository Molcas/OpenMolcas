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
* Copyright (C) 1992,2007, Roland Lindh                                *
************************************************************************
      SubRoutine PGet2_CD3(iCmp,IndShl,iBas,jBas,kBas,lBas,
     &                  Shijij, iAO, iAOst, nijkl,PSO,nPSO,
     &                  DSO,DSSO,nDSO,ExFac,CoulFac,PMax,V_k,mV_k)
************************************************************************
*  Object: to assemble the 2nd order density matrix of a SCF wave      *
*          function from the 1st order density matrix.                 *
*                                                                      *
*          The indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
* Called from: PGet0                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             January '92.                                             *
*                                                                      *
*             Modified for Cholesky 1-center gradients May 2007 by     *
*             R. Lindh                                                 *
*                                                                      *
************************************************************************
      use SOAO_Info, only: iAOtSO
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "lundio.fh"
#include "print.fh"
#include "WrkSpc.fh"
      Real*8 PSO(nijkl,nPSO), DSO(nDSO), DSSO(nDSO), V_k(mV_k)
      Integer iCmp(4), iAO(4), iAOst(4), IndShl(4)
      Logical Shijij
*     Local Array
      Integer iSym(0:7), jSym(0:7), kSym(0:7), lSym(0:7)
      Integer iTwoj(0:7)
      Data iTwoj/1,2,4,8,16,32,64,128/
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 39
      iPrint = nPrint(iRout)
#ifdef _DEBUG_
      Call qEnter('PGet2_CD3')
      iPrint=99
      If (iPrint.ge.99) Then
         iComp = 1
         Call PrMtrx(' In PGet2_CD3:DSO ',[iD0Lbl],iComp,1,D0)
         Call RecPrt('V_K',' ',V_K,1,mV_K)
      End If
#endif
*                                                                      *
************************************************************************
*                                                                      *
C     Fac = One / Four
      Fac = One / Two
      lOper=1
      PMax=Zero
*
*-----Quadruple loop over elements of the basis functions angular
*     description.
*     Observe that we will walk through the memory in AOInt in a
*     sequential way.
*
      MemSO2 = 0
      Do 100 i1 = 1, iCmp(1)
         niSym = 0
         Do 101 j = 0, nIrrep-1
            If (iAnd(IrrCmp(IndShl(1)+i1),
     &          iTwoj(j)).ne.0) Then
               iSym(niSym) = j
               niSym = niSym + 1
            End if
101      Continue
         Do 200 i2 = 1, iCmp(2)
            njSym = 0
            Do 201 j = 0, nIrrep-1
               If (iAnd(IrrCmp(IndShl(2)+i2),
     &             iTwoj(j)).ne.0) Then
                  jSym(njSym) = j
                  njSym = njSym + 1
               End If
201         Continue
            Do 300 i3 = 1, iCmp(3)
               nkSym = 0
               Do 301 j = 0, nIrrep-1
                  If (iAnd(IrrCmp(IndShl(3)+i3),
     &                iTwoj(j)).ne.0) Then
                     kSym(nkSym) = j
                     nkSym = nkSym + 1
                  End If
301            Continue
               Do 400 i4 = 1, iCmp(4)
                  nlSym = 0
                  Do 401 j = 0, nIrrep-1
                     If (iAnd(IrrCmp(IndShl(4)+i4),
     &                   iTwoj(j)).ne.0) Then
                        lSym(nlSym) = j
                        nlSym = nlSym + 1
                     End If
401               Continue
*
*------Loop over irreps which are spanned by the basis function.
*
       Do 110 is = 0, niSym-1
          j1 = iSym(is)
*
          Do 210 js = 0, njSym-1
             j2 = jSym(js)
             j12 = iEor(j1,j2)
*
             Do 310 ks = 0, nkSym-1
                j3 = kSym(ks)
                j123 = iEor(j12,j3)
                Do 410 ls = 0, nlSym-1
                   j4 = lSym(ls)
                   If (j123.ne.j4) Go To 410
*
                MemSO2 = MemSO2 + 1
*
*               Unfold the way the eight indicies have been reordered.
                iSO = iAOtSO(iAO(1)+i1,j1)+iAOst(1)
                jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)
                kSO = iAOtSO(iAO(3)+i3,j3)+iAOst(3)
                lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)
*
                If (j1.ne.j2 .and. j1.ne.j3 .and. j1.ne.j4) Then
*------------------all irreps are different and the 2nd order density
*                  matrix will be identical to zero for a SCF type wave
*                  function.
                   call dcopy_(nijkl,[Zero],0,PSO(1,MemSO2),1)
                   Go To 310
                End If
*
                mijkl = 0
                Do 120 lAOl = 0, lBas-1
                   lSOl = lSO + lAOl
                   Do 220 kAOk = 0, kBas-1
                      kSOk = kSO + kAOk
                      Do 320 jAOj = 0, jBas-1
                         jSOj = jSO + jAOj
                         Do 420 iAOi = 0, iBas-1
                            iSOi = iSO + iAOi
                            mijkl = mijkl + 1
*
*---------------------------Contribution V_k(ij)*D(kl) to P(ijkl)
                            If (j1.eq.j2) Then
*------------------------------j3.eq.j4 also
                               Indi=Max(iSOi,jSOj)
                               Indj=iSOi+jSOj-Indi
                               Indk=Max(kSOk,lSOl)
                               Indl=kSOk+lSOl-Indk
                               iPntij=iPntSO(j1,j2,lOper,nbas)
                               iPntkl=iPntSO(j3,j4,lOper,nbas)
                               Indij=iPntij+(Indi-1)*Indi/2+Indj
                               Indkl=iPntkl+(Indk-1)*Indk/2+Indl
                               temp=V_k(Indij)*DSO(Indkl)*Coulfac
                            Else
                               temp = Zero
                            End If
*
*---------------------------Contribution -1/4*D(ik)*D(jl) to P(ijkl)
C                           If (j1.eq.j3) Then
*------------------------------j2.eq.j4 also
C                           End If
*
*---------------------------Contribution -1/4*D(il)*D(jk) to P(ijkl)
C                           If (j1.eq.j4) Then
*------------------------------j2.eq.j3 also
C                           End If
*
                            PMax=Max(PMax,Abs(Temp))
                            PSO(mijkl,MemSO2) =  Fac * temp
*
 420                     Continue
 320                  Continue
 220               Continue
 120            Continue
*
 410            Continue
 310         Continue
 210      Continue
 110   Continue
*
*
 400           Continue
 300        Continue
 200     Continue
 100  Continue
      If (nPSO.ne.MemSO2) Then
        Write (6,*) ' PGet2_CD3: nPSO.ne.MemSO2'
        Write (6,*) nPSO, MemSO2
        Call Abend
      End If
*
#ifdef _DEBUG_
      If (iPrint.ge.99) Then
         Call RecPrt(' In PGet2_CD3:PSO ',' ',PSO,nijkl,nPSO)
      End If
      Call GetMem(' Exit PGet2_CD3','CHECK','REAL',iDum,iDum)
      Call qExit('PGet2_CD3')
#endif
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_logical(Shijij)
         Call Unused_real_array(DSSO)
         Call Unused_real(ExFac)
      End If
      End
