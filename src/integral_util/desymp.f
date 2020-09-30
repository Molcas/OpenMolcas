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
* Copyright (C) 1990,1992, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      Subroutine DesymP(iAng, iCmp, jCmp, kCmp, lCmp, Shijij, iShll,
     &                  iShell, iAO, kOp, ijkl,Aux,nAux,PAO,PSO,nPSO)
************************************************************************
*                                                                      *
*  Object: to transform the integrals in AO basis to symmetry adapted  *
*          orbitals , SO. This is done by accumulating the AO inte-    *
*          grals onto the SO integrals.                                *
*          Observe that one of the operator is the Unit operation      *
*          performed on center A. However, since we scramble the order *
*          we do not really know which center is which. However, the   *
*          Unit operator will always give a factor of one. Hence, this *
*          is a convenient way to do the symmetry transformation with- *
*          out having to know the order.                               *
*                                                                      *
*          The indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             Modified to desymmetrization of the second order density *
*             matrix, January '92.                                     *
************************************************************************
      use Basis_Info
      use Symmetry_Info, only: nIrrep, iChTbl, iOper, iChBas
      use SOAO_Info, only: iAOtSO
      use Real_Spherical, only: iSphCr
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 PAO(ijkl,iCmp,jCmp,kCmp,lCmp), PSO(ijkl,nPSO), Aux(nAux)
      Logical Shij, Shkl, Shijij
      Integer iAng(4), iShell(4), kOp(4), iShll(4), iAO(4)
*     Local Array
      Integer iSym(0:7), jSym(0:7), kSym(0:7), lSym(0:7)
      Integer iTwoj(0:7)
      Real*8 Prmt(0:7)
      Data iTwoj/1,2,4,8,16,32,64,128/
      Data Prmt/1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0/
*
*     Statement Function
*
      xPrmt(i,j) = Prmt(iAnd(i,j))
*
      iRout = 38
      iPrint = nPrint(iRout)
      Shij = iShell(1).eq.iShell(2)
      Shkl = iShell(3).eq.iShell(4)
      MemSO2 = 1
*
#ifdef _DEBUG_
      If (iPrint.ge.99) Then
         Call RecPrt(' In DesymP: PSO ',' ',PSO,ijkl,nPSO)
         Call WrCheck(' In DesymP: PSO ',PSO,ijkl*nPSO)
         Write (6,*) 'iCmp,jCmp,kCmp,lCmp,nPSO=',
     &                iCmp,jCmp,kCmp,lCmp,nPSO
         Write (6,*) Shij, Shkl, Shijij
         Write (6,*) 'kOp=',kOp
      End If
#endif
      Fact=Eight
      If (Shij)   Fact=Fact*Half
      If (Shkl)   Fact=Fact*Half
      If (Shijij) Fact=Fact*Half
*
*-----Initialize second order density matrix in AO basis
*
      call dcopy_(ijkl*iCmp*jCmp*kCmp*lCmp,[Zero],0,PAO,1)
*
*
*     Quadruple loop over elements of the basis functions angular
*     description. Loops are reduced to just produce unique SO integrals
*     Observe that we will walk through the memory in PAO in a
*     sequential way.
*
      ii = iAng(1)*(iAng(1)+1)*(iAng(1)+2)/6
      jj = iAng(2)*(iAng(2)+1)*(iAng(2)+2)/6
      kk = iAng(3)*(iAng(3)+1)*(iAng(3)+2)/6
      ll = iAng(4)*(iAng(4)+1)*(iAng(4)+2)/6
      Do 100 i1 = 1, iCmp
         iChBs = iChBas(ii+i1)
         If (Shells(iShll(1))%Transf) iChBs = iChBas(iSphCr(ii+i1))
         pa = xPrmt(iOper(kOp(1)),iChBs)
         niSym=0
         Do 101 j = 0, nIrrep-1
            If (iAOtSO(iAO(1)+i1,j)>0) Then
               iSym(niSym)=j
               niSym=niSym+1
            End If
101      Continue
         Do 200 i2 = 1, jCmp
            jChBs = iChBas(jj+i2)
            If (Shells(iShll(2))%Transf) jChBs = iChBas(iSphCr(jj+i2))
            pb = xPrmt(iOper(kOp(2)),jChBs)
            njSym=0
            Do 201 j = 0, nIrrep-1
               If (iAOtSO(iAO(2)+i2,j)>0) Then
                  jSym(njSym) = j
                  njSym=njSym+1
               End If
201         Continue
            Do 300 i3 = 1, kCmp
               kChBs = iChBas(kk+i3)
               If (Shells(iShll(3))%Transf)
     &            kChBs = iChBas(iSphCr(kk+i3))
               pc = xPrmt(iOper(kOp(3)),kChBs)
               nkSym=0
               Do 301 j = 0, nIrrep-1
                  If (iAOtSO(iAO(3)+i3,j)>0) Then
                     kSym(nkSym) = j
                     nkSym=nkSym+1
                   End If
301            Continue
               Do 400 i4 = 1, lCmp
                  lChBs = iChBas(ll+i4)
                  If (Shells(iShll(4))%Transf)
     &               lChBs = iChBas(iSphCr(ll+i4))
*-----------------Parity factor due to symmetry operations applied to the
*                 angular part of the basis functions.
                  FactNs = pa*pb*pc * xPrmt(iOper(kOp(4)),lChBs)
                  nlSym=0
                  Do 401 j = 0, nIrrep-1
                     If (iAOtSO(iAO(4)+i4,j)>0) Then
                        lSym(nlSym)=j
                        nlSym=nlSym+1
                     End If
401               Continue
*
*------Loop over irreps which are spanned by the basis functions.
*      Again, the loop structure is restricted to ensure unique
*      integrals.
*
*
      If (nIrrep.eq.1) Then
*--------FactNs=1
         Call DYaX(ijkl,Fact,PSO(1,MemSO2),1,PAO(1,i1,i2,i3,i4),1)
         MemSO2 = MemSO2 + 1
         Go To 400
      End If
*
      iAux = 0
      Do 110 is = 0, niSym-1
         j1=iSym(is)
         Xa = DBLE(iChTbl(j1,kOp(1))) * FactNs
         Do 210 js = 0, njSym-1
            j2=jSym(js)
            Xb = DBLE(iChTbl(j2,kOp(2))) * Xa
            j12 = iEor(j1,j2)
            Do 310 ks = 0, nkSym-1
               j3=kSym(ks)
               Xg = DBLE(iChTbl(j3,kOp(3))) * Xb
               j123=iEor(j12,j3)
               Do 320 ls = 0, nlSym-1
                  j4=lSym(ls)
                  If (j123.ne.j4) Go To 320
*
                  iAux = iAux + 1
                  Aux(iAux) = DBLE(iChTbl(j4,kOp(4))) * Xg *
     &                        Fact
                  Go To 310
*
 320           Continue
 310        Continue
 210     Continue
 110  Continue
*
      If (iAux.ne.0) Then
#ifdef _DEBUG_
         If (iPrint.ge.99) Call RecPrt(' Aux',' ',Aux,iAux,1)
#endif
         If (iAux.ne.1) Then
            Call DNaXpY(iAux,ijkl,Aux,1,PSO(1,MemSO2),1,ijkl,
     &                  PAO(1,i1,i2,i3,i4),1,0)
         Else
            Call DaXpY_(ijkl,Aux(1),PSO(1,MemSO2),1,
     &                 PAO(1,i1,i2,i3,i4),1)
         End If
         MemSO2 = MemSO2 + iAux
      End If
*
 400           Continue
 300        Continue
 200     Continue
 100  Continue
*
#ifdef _DEBUG_
      If (iPrint.ge.99) Then
         Call RecPrt(' On exit from DesymP: PAO ',' ', PAO,
     &            ijkl,iCmp*jCmp*kCmp*lCmp)
         Do 3333 i = 1, ijkl
            Write (6,*) DDot_(iCmp*jCmp*kCmp*lCmp,
     &                  PAO(i,1,1,1,1),ijkl,
     &                  PAO(i,1,1,1,1),ijkl)
 3333    Continue
        Call WrCheck('DesymP: PAO ',PAO,ijkl*iCmp*jCmp*kCmp*lCmp)
        write (6,*) ddot_(ijkl*iCmp*jCmp*kCmp*lCmp,PAO,1,[One],0)
        write (6,*) ddot_(ijkl*iCmp*jCmp*kCmp*lCmp,PAO,1,PAO,1)
      End If
#endif
      Return
      End
