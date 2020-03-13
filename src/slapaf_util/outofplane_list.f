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
      Subroutine OutOfPlane_List(
     &                 nq,
     &                 nAtoms,iIter,nIter,Cx,iOper,nSym,jStab,
     &                 nStab,nDim,Smmtrc,Process,Value,
     &                 nB,iANr,qLbl,iRef,
     &                 fconst,rMult,LuIC,Name,Indq,iPrv,Proc_dB,
     &                 iTabBonds,nBonds,iTabAI,mAtoms,iTabAtoms,nMax,
     &                 mB_Tot,mdB_Tot,BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,
     &                 nqB)
************************************************************************
*     This is a quick and possibly dirty implementation of the out-    *
*     of-plane angle. RL, Tokyo June, 2004.                            *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "print.fh"
      Parameter (mB = 4*3)
      Real*8 Cx(3,nAtoms,nIter), A(3,4), Grad(mB), Hess(mB**2),
     &       fconst(nB), Value(nB,nIter),
     &       Ref(3,4), Prv(3,4), rMult(nB),
     &       Grad_ref(9), RX4Y(3,3), BM(nB_Tot), dBM(ndB_Tot)
      Integer   nStab(nAtoms), iOper(0:nSym-1), iANr(nAtoms),
     &          iDCRR(0:7), jStab(0:7,nAtoms), iPhase(3,0:7),
     &          iStabM(0:7), Ind(4), iDCR(4), iDCRT(0:7),
     &          iDCRS(0:7), iStabN(0:7), iStabO(0:7), iChOp(0:7),
     &          Indq(3,nB), iDCRX(0:7), iDCRY(0:7), nqB(nB),
     &          iTabBonds(3,nBonds), iTabAI(2,mAtoms),
     &          iTabAtoms(2,0:nMax,mAtoms), iBM(nB_Tot), idBM(2,ndB_Tot)
      Logical Smmtrc(3,nAtoms), Process, PSPrint,
     &        MinBas, Help, Proc_dB, R_Stab_A
      Character*14 Label, qLbl(nB)
      Character*3 ChOp(0:7)
#include "Molcas.fh"
      Character*(LENIN) Name(nAtoms)
      Character*(LENIN4) Lbls(4)
#include "bondtypes.fh"
#define _FMIN_
#include "ddvdt.fh"
#include "ddvdt_outofp.fh"
      Data iPhase/ 1, 1, 1,   -1, 1, 1,   1,-1, 1,  -1,-1, 1,
     &             1, 1,-1,   -1, 1,-1,   1,-1,-1,  -1,-1,-1/
      Data ChOp/'E  ','X ','Y ','XY ','Z  ','XZ ','YZ ','XYZ'/
      Data iChOp/1,1,1,2,1,2,2,3/
#include "constants.fh"
*                                                                      *
************************************************************************
*                                                                      *
*#define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
*
      If (nBonds.lt.3) Return
      iRout=152
      iPrint=nPrint(iRout)
      Call QEnter('OutOfPs')
      nqO=0
      PSPrint=.False.
      Call FZero(Hess,144)
#ifdef _DEBUG_
      iPrint=99
      If (iPrint.ge.99) PSPrint=.True.
#endif
*
*---- Loop over out-of-plane angles.
*
      bohr=CONST_BOHR_RADIUS_IN_SI_ * 1.0D+10
      MinBas=.False.
      If (MinBas) Then
         Fact=1.3d0
      Else
         Fact=One
      End If
      nCent=4
*
************************************************************************
*     Notes relevant to the out-of-plane implementation                *
*                                                                      *
*    Connection is 1-4, 2-4, and 3-4. We renumber them according to    *
*    j-i, k-i, and l-i                                                 *
*                                                                      *
************************************************************************
*
*     Order will play a role here. That is, the torsion
*     A-R(B)-T(C)-TS(D) is NOT identical to
*     A-R(B)-TS(D)-T(C). Hence we put no restriction on the
*     pairs AB and CD. However, for the pair of pairs we have
*     that order is irrelevant, i.e. ABCD is identical to
*     DCBA. To garantee this we limit the pairs to the unique
*     combinations.
*
      Do jBond = 1, nBonds
         jBondType=iTabBonds(3,jBond)
         If (jBondType.eq.vdW_Bond) Go To 101
         If (jBondType.gt.Magic_Bond) Go To 101
*
         Do iCase = 1, 2
*
            If (iCase.eq.1) Then
               iAtom_ = iTabBonds(1,jBond)
               jAtom_ = iTabBonds(2,jBond)
            Else
               iAtom_ = iTabBonds(2,jBond)
               jAtom_ = iTabBonds(1,jBond)
            End If
            iAtom=iTabAI(1,iAtom_)
            jAtom=iTabAI(1,jAtom_)
            ir = iTabRow(iANr(iAtom))
            jr = iTabRow(iANr(jAtom))
            Ind(1) = jAtom
            Ind(4) = iAtom
*
            Help = ir.gt.3.or.jr.gt.3
            iDCR(4)=iTabAI(2,iAtom_)
            iDCR(1)=iTabAI(2,jAtom_)
#ifdef _DEBUG_
            Write (6,*)
            Write (6,*) 'E,R=',Name(iAtom),ChOp(iDCR(4)),
     &                         Name(jAtom),ChOp(iDCR(1))
#endif
            nCoBond_j=nCoBond(jAtom_,mAtoms,nMax,iTabBonds,nBonds,
     &                        nBonds,iTabAtoms)
            nFgBond_j=nFgBond(jAtom_,mAtoms,nMax,iTabBonds,nBonds,
     &                        nBonds,iTabAtoms)
            If (nCoBond_j.gt.1.and.nFgBond_j.eq.0) Go To 201
            If (iDCR(4).ne.iOper(0)) Go To 201
*
*           R
*
            If (R_Stab_A(iDCR(1),jStab(0,iAtom),nStab(iAtom)) .and.
     &          iDCR(1).ne.iOper(0)) Go To 201
*
            call dcopy_(3,Cx(1,iAtom,iIter),1,A(1,4),  1)
            call dcopy_(3,Cx(1,iAtom,iRef), 1,Ref(1,4),1)
            call dcopy_(3,Cx(1,iAtom,iPrv), 1,Prv(1,4),1)
*
*---------- Form double coset representatives for (iAtom,jAtom)
*
            Call DCR(Lambda,iOper,nSym,
     &               jStab(0,iAtom),nStab(iAtom),
     &               jStab(0,jAtom),nStab(jAtom),
     &               iDCRR,nDCRR)
            kDCRR=iDCR(1)
#ifdef _DEBUG_
            If (PSPrint) Then
               Write (6,'(10A)') 'U={',(ChOp(jStab(i,iAtom)),
     &                            i=0,nStab(iAtom)-1),'}  '
               Write (6,'(10A)') 'V={',(ChOp(jStab(i,jAtom)),
     &                            i=0,nStab(jAtom)-1),'}  '
               Write (6,'(10A)') 'R={',(ChOp(iDCRR(i)),
     &                            i=0,nDCRR-1),'}  '
               Write (6,'(2A)') 'R=',ChOp(kDCRR)
            End If
#endif
*
            A(1,1)   = DBLE(iPhase(1,kDCRR))*Cx(1,jAtom,iIter)
            A(2,1)   = DBLE(iPhase(2,kDCRR))*Cx(2,jAtom,iIter)
            A(3,1)   = DBLE(iPhase(3,kDCRR))*Cx(3,jAtom,iIter)
            Ref(1,1) = DBLE(iPhase(1,kDCRR))*Cx(1,jAtom,iRef)
            Ref(2,1) = DBLE(iPhase(2,kDCRR))*Cx(2,jAtom,iRef)
            Ref(3,1) = DBLE(iPhase(3,kDCRR))*Cx(3,jAtom,iRef)
            Prv(1,1) = DBLE(iPhase(1,kDCRR))*Cx(1,jAtom,iPrv)
            Prv(2,1) = DBLE(iPhase(2,kDCRR))*Cx(2,jAtom,iPrv)
            Prv(3,1) = DBLE(iPhase(3,kDCRR))*Cx(3,jAtom,iPrv)
*
*---------- Form stabilizer for (iAtom,jAtom)
*
            Call Inter(jStab(0,iAtom),nStab(iAtom),
     &                 jStab(0,jAtom),nStab(jAtom),
     &                     iStabM,nStabM)
#ifdef _DEBUG_
            If (PSPrint) Then
               Write (6,'(10A)') 'M={',
     &               (ChOp(iStabM(i)),i=0,nStabM-1),'}  '
            End If
#endif
*
            If (Help) Then
               f_Const_ij_Ref=rko
               f_Const_ij    =rko
            Else
               r0=rAV(ir,jr)
               Alpha=aAv(ir,jr)
               rij2_Ref=(Ref(1,4)-Ref(1,1))**2
     &                 +(Ref(2,4)-Ref(2,1))**2
     &                 +(Ref(3,4)-Ref(3,1))**2
               f_Const_ij_Ref=rko*Exp(Alpha*(r0**2-rij2_Ref))
               rij2=(A(1,4)-A(1,1))**2
     &             +(A(2,4)-A(2,1))**2
     &             +(A(3,4)-A(3,1))**2
               f_Const_ij=rko*Exp(Alpha*(r0**2-rij2))
            End If
*
*
            nNeighbor_i = iTabAtoms(1,0,iAtom_)
            Do kNeighbor = 1, nNeighbor_i
               kAtom_=iTabAtoms(1,kNeighbor,iAtom_)
               If (kAtom_.eq.jAtom_) Go To 301
               kBond =iTabAtoms(2,kNeighbor,iAtom_)
               kBondType=iTabBonds(3,kBond)
#ifdef _DEBUG_
               Write (6,*) 'kBond,kBondType=',
     &                      kBond,kBondType
#endif
               If (kBondType.eq.vdW_Bond) Go To 301
               If (kBondType.gt.Magic_Bond) Go To 301
               If (kBond.eq.jBond) Go To 301
*
               kAtom=iTabAI(1,kAtom_)
               ik_=nAtoms*(kAtom-1)+iAtom
               kr = iTabRow(iANr(kAtom))
               Ind(2) = kAtom
               iDCR(2)=iTabAI(2,kAtom_)
*
               If (iDCR(1).eq.iOper(0)) Then
                  If (R_Stab_A(iDCR(2),jStab(0,iAtom),nStab(iAtom)).and.
     &                R_Stab_A(iDCR(2),jStab(0,jAtom),nStab(jAtom)).and.
     &                iDCR(2).ne.iOper(0)) Go To 301
               Else
                  If (R_Stab_A(iDCR(2),jStab(0,iAtom),nStab(iAtom)).and.
     &                iDCR(2).ne.iOper(0)) Go To 301
               End If
#ifdef _DEBUG_
               Write (6,*)
               Write (6,*) 'T=',Name(kAtom),ChOp(iDCR(2))
               Write (6,*) 'kAtom=', kAtom
#endif
*
               Do lNeighbor = 1, nNeighbor_i
                  lAtom_=iTabAtoms(1,lNeighbor,iAtom_)
                  If (lAtom_.eq.jAtom_) Go To 401
                  If (lAtom_.le.kAtom_) Go To 401
                  lBond =iTabAtoms(2,lNeighbor,iAtom_)
                  lBondType=iTabBonds(3,lBond)
                  If (lBondType.eq.vdW_Bond)   Go To 401
                  If (lBondType.gt.Magic_Bond)   Go To 401
                  If (lBond.eq.jBond)   Go To 401
                  If (lBond.eq.kBond)   Go To 401
                  lAtom=iTabAI(1,lAtom_)
#ifdef _DEBUG_
                  Write (6,*) 'lBond,lBondType=',
     &                         lBond,lBondType
                  Write (6,*) 'lAtom=', lAtom
#endif
*
                  il_=nAtoms*(lAtom-1)+iAtom
                  lr = iTabRow(iANr(lAtom))
                  Ind(3) = lAtom
                  iDCR(3)=iTabAI(2,lAtom_)
C                 If (kAtom.gt.lAtom) Go To 401
*
               If (iDCR(1).eq.iOper(0)) Then
                  If (R_Stab_A(iDCR(3),jStab(0,iAtom),nStab(iAtom)).and.
     &                R_Stab_A(iDCR(3),jStab(0,jAtom),nStab(jAtom)).and.
     &                iDCR(3).ne.iOper(0).and.iDCR(2).ne.iOper(0))
     &                Go To 401
               Else
                  If (R_Stab_A(iDCR(2),jStab(0,iAtom),nStab(iAtom)).and.
     &                iDCR(3).ne.iOper(0).and.iDCR(2).ne.iOper(0))
     &                Go To 401
               End If

               If (kAtom.eq.lAtom) Then
                  If (R_Stab_A(iDCR(2),jStab(0,iAtom),nStab(iAtom)).and.
     &                R_Stab_A(iDCR(2),jStab(0,jAtom),nStab(jAtom)).and.
     &                iDCR(2).ne.iOper(0)) Go To 401
                  If (iDCR(3).eq.iOper(0)) Go To 401
                  If (R_Stab_A(iDCR(3),jStab(0,iAtom),nStab(iAtom)).and.
     &                R_Stab_A(iDCR(3),jStab(0,jAtom),nStab(jAtom)).and.
     &                iDCR(3).ne.iOper(0).and.iDCR(2).ne.iOper(0))
     &                Go To 401
               End If
#ifdef _DEBUG_
                  Write (6,*)
                  Write (6,*) 'TS=',Name(lAtom),ChOp(iDCR(3))
#endif
*
                  Help = ir.gt.3.or.jr.gt.3.or.kr.gt.3.or.lr.gt.3
*
                  Write (Label,'(A,I2,A,I2,A,I2,A,I2,A)')
     &                   'D(',iAtom,',',jAtom,',',kAtom,',',lAtom,')'
*
*---------------- Form double coset representatives for (kAtom,lAtom)
*
                  Call DCR(Lambda,iOper,nSym,
     &                     jStab(0,kAtom),nStab(kAtom),
     &                     jStab(0,lAtom),nStab(lAtom),
     &                     iDCRS,nDCRS)
                  kDCRS=iEor(iDCR(2),iDCR(3))
*
#ifdef _DEBUG_
                  If (PSPrint) Then
                     Write (6,'(10A)') 'W={',(ChOp(jStab(i,kAtom)),
     &                                  i=0,nStab(kAtom)-1),'}  '
                     Write (6,'(10A)') 'X={',(ChOp(jStab(i,lAtom)),
     &                                  i=0,nStab(lAtom)-1),'}  '
                     Write (6,'(10A)') 'S={',(ChOp(iDCRS(i)),
     &                                  i=0,nDCRS-1),'}  '
                     Write (6,'(2A)') 'S=',ChOp(kDCRS)
                  End If
#endif
*
                  Ref(1,2) =                        Cx(1,kAtom,iRef)
                  Ref(2,2) =                        Cx(2,kAtom,iRef)
                  Ref(3,2) =                        Cx(3,kAtom,iRef)
                  Ref(1,3) = DBLE(iPhase(1,kDCRS))* Cx(1,lAtom,iRef)
                  Ref(2,3) = DBLE(iPhase(2,kDCRS))* Cx(2,lAtom,iRef)
                  Ref(3,3) = DBLE(iPhase(3,kDCRS))* Cx(3,lAtom,iRef)
                  Prv(1,2) =                        Cx(1,kAtom,iPrv)
                  Prv(2,2) =                        Cx(2,kAtom,iPrv)
                  Prv(3,2) =                        Cx(3,kAtom,iPrv)
                  Prv(1,3) = DBLE(iPhase(1,kDCRS))* Cx(1,lAtom,iPrv)
                  Prv(2,3) = DBLE(iPhase(2,kDCRS))* Cx(2,lAtom,iPrv)
                  Prv(3,3) = DBLE(iPhase(3,kDCRS))* Cx(3,lAtom,iPrv)
*
*---------------- Form stabilizer for (kAtom,lAtom)
*
                  Call Inter(jStab(0,kAtom),nStab(kAtom),
     &                       jStab(0,lAtom),nStab(lAtom),
     &                       iStabN,nStabN)
*
#ifdef _DEBUG_
                  If (PSPrint) Then
                     Write (6,'(10A)') 'N={',
     &                     (ChOp(iStabN(i)),i=0,nStabN-1),'}  '
                  End If
#endif
*
*---------------- Form double coset representatives for
*                 ((iAtom,jAtom),(kAtom,lAtom))
*
                  Call DCR(Lambda,iOper,nSym,
     &                     iSTabM,nStabM,
     &                     iStabN,nStabN,
     &                     iDCRT,nDCRT)
*
*---------------- Take care of some special cases which normally
*                 are not included. If A=B we will normally exclude
*                 the pairs R(A)-A and TS(C)-T(C).
*
                  Call iCopy(nDCRT,iDCRT,1,iDCRX,1)
                  Call iCopy(nDCRT,iDCRT,1,iDCRY,1)
                  nDCRX=nDCRT
                  nDCRY=nDCRT
                  If (iAtom.eq.jAtom) Then
*                    Write (*,*) ' Special fix'
                     Call Union(iDCRX,nDCRX,iDCRY,nDCRY,
     &                          kDCRR,iDCRT,nDCRT)
                  Else If (kAtom.eq.lAtom) Then
*                    Write (*,*) ' Special fix'
                     Call Union(iDCRX,nDCRX,iDCRY,nDCRY,
     &                          kDCRS,iDCRT,nDCRT)
                  End If
*
                  kDCRT =iDCR(2)
                  kDCRTS=iDCR(3)
*
#ifdef _DEBUG_
                  If (PSPrint) Then
                     Write (6,'(10A)') 'T={',
     &                     (ChOp(iDCRT(i)),i=0,nDCRT-1),'}  '
                     Write (6,'(2A)') 'T=',ChOp(kDCRT)
                  End If
#endif
*
                  A(1,2)   = DBLE(iPhase(1,kDCRT ))*Cx(1,kAtom,iIter)
                  A(2,2)   = DBLE(iPhase(2,kDCRT ))*Cx(2,kAtom,iIter)
                  A(3,2)   = DBLE(iPhase(3,kDCRT ))*Cx(3,kAtom,iIter)
                  Ref(1,2) = DBLE(iPhase(1,kDCRT ))*Cx(1,kAtom,iRef)
                  Ref(2,2) = DBLE(iPhase(2,kDCRT ))*Cx(2,kAtom,iRef)
                  Ref(3,2) = DBLE(iPhase(3,kDCRT ))*Cx(3,kAtom,iRef)
                  Prv(1,2) = DBLE(iPhase(1,kDCRT ))*Cx(1,kAtom,iPrv)
                  Prv(2,2) = DBLE(iPhase(2,kDCRT ))*Cx(2,kAtom,iPrv)
                  Prv(3,2) = DBLE(iPhase(3,kDCRT ))*Cx(3,kAtom,iPrv)
                  A(1,3)   = DBLE(iPhase(1,kDCRTS))*Cx(1,lAtom,iIter)
                  A(2,3)   = DBLE(iPhase(2,kDCRTS))*Cx(2,lAtom,iIter)
                  A(3,3)   = DBLE(iPhase(3,kDCRTS))*Cx(3,lAtom,iIter)
                  Ref(1,3) = DBLE(iPhase(1,kDCRTS))*Cx(1,lAtom,iRef)
                  Ref(2,3) = DBLE(iPhase(2,kDCRTS))*Cx(2,lAtom,iRef)
                  Ref(3,3) = DBLE(iPhase(3,kDCRTS))*Cx(3,lAtom,iRef)
                  Prv(1,3) = DBLE(iPhase(1,kDCRTS))*Cx(1,lAtom,iPrv)
                  Prv(2,3) = DBLE(iPhase(2,kDCRTS))*Cx(2,lAtom,iPrv)
                  Prv(3,3) = DBLE(iPhase(3,kDCRTS))*Cx(3,lAtom,iPrv)
*
*---------------- Form the stabilizer for the out-of-plane
*
                  If (iAtom.eq.lAtom .and.
     &                jAtom.eq.kAtom.and.kDCRR.eq.kDCRS) Then
                     Call Union(iStabM,nStabM,
     &                          iStabN,nStabN,
     &                          kDCRTS,iStabO,nStabO)
                  Else
                     Call Inter(iStabM,nStabM,
     &                          iStabN,nStabN,
     &                          iStabO,nStabO)
                  End If
*
#ifdef _DEBUG_
                  If (PSPrint) Then
                     Write (6,'(10A)') 'M={',
     &                     (ChOp(iStabM(i)),i=0,nStabM-1),'}  '
                     Write (6,'(10A)') 'N={',
     &                     (ChOp(iStabN(i)),i=0,nStabN-1),'}  '
                     Write (6,'(10A)') 'O={',
     &                     (ChOp(iStabO(i)),i=0,nStabO-1),'}  '
                  End If
                  Write (6,*) 'jAtom,iAtom,kAtom,lAtom=',
     &                         jAtom,iAtom,kAtom,lAtom
#endif
*
*-----------------Compute the degeneracy of the torsion
*
                  iDeg=nSym/nStabO
                  Deg=Sqrt(DBLE(iDeg))
*
*-----------------Test if coordinate should be included
*
                  If (Help) Then
                     f_Const_ijk_Ref=f_Const_ij_Ref
                     f_Const_ijk=f_Const_ij
                     f_Const_Ref=f_Const_ijk_Ref
                     f_Const=f_Const_ijk
                  Else
*
*                    Test the ik-pair
*
                     r0=rAV(ir,kr)
                     Alpha=aAv(ir,kr)
                     rik2_Ref=(Ref(1,4)-Ref(1,2))**2
     &                       +(Ref(2,4)-Ref(2,2))**2
     &                       +(Ref(3,4)-Ref(3,2))**2
                     f_Const_ijk_Ref=f_Const_ij_Ref
     &                          *Exp(Alpha*(r0**2-rik2_Ref))
                     rik2=(A(1,4)-A(1,2))**2
     &                   +(A(2,4)-A(2,2))**2
     &                   +(A(3,4)-A(3,2))**2
                     f_Const_ijk=f_Const_ij
     &                          *Exp(Alpha*(r0**2-rik2))
*
*                    Test the il-pair
*
                              r0=rAV(ir,lr)
                     Alpha=aAv(ir,lr)
                     ril2_Ref=(Ref(1,4)-Ref(1,3))**2
     &                       +(Ref(2,4)-Ref(2,3))**2
     &                       +(Ref(3,4)-Ref(3,3))**2
                     f_Const_Ref=f_Const_ijk_Ref
     &                      *Exp(Alpha*(r0**2-ril2_Ref))
                     ril2=(A(1,4)-A(1,3))**2
     &                   +(A(2,4)-A(2,3))**2
     &                   +(A(3,4)-A(3,3))**2
                     f_Const=f_Const_ijk
     &                      *Exp(Alpha*(r0**2-ril2))
                  End If
                  If (f_Const_Ref.lt.f_Const_Min .and.
     &                jBondtype.ne.Fragments_Bond .and.
     &                kBondtype.ne.Fragments_Bond .and.
     &                lBondtype.ne.Fragments_Bond) Go To 401
*
*---------------- Check that valence angles are above threshold
*
                  mCent=3
                  delta0 = (45.0D0/180.D0)*Pi
*                 If (jBondType.eq.Fragments_Bond .or.
*    &                kBondType.eq.Fragments_Bond .or.
*    &                lBondType.eq.Fragments_Bond) delta=0.0D0
*
*---------------- 1-4-2
*
                  call dcopy_(3,Ref(1,1),1,RX4Y(1,1),1)
                  call dcopy_(3,Ref(1,4),1,RX4Y(1,2),1)
                  call dcopy_(3,Ref(1,2),1,RX4Y(1,3),1)
                  Call Bend(RX4Y,mCent,Fi2,Grad_ref,
     &                     .False.,
     &                     .False.,'        ',Hess,.False.)
#ifdef _DEBUG_
                  Write (6,*) '1-4-2: Fi2=',Fi2
#endif
                  delta = delta0
                  If (jBondType.eq.Fragments_Bond .or.
     &                kBondType.eq.Fragments_Bond) delta=0.0D0
                  If (Fi2.gt.Pi-delta) Go To 401
                  If (Fi2.lt.delta)    Go To 401
*
*---------------- 1-4-3
*
                  call dcopy_(3,Ref(1,3),1,RX4Y(1,3),1)
                  Call Bend(RX4Y,mCent,Fi3,Grad_ref,
     &                      .False.,
     &                      .False.,'        ',Hess,.False.)
#ifdef _DEBUG_
                  Write (6,*) '1-4-3: Fi3=',Fi3
#endif
                  delta = delta0
                  If (jBondType.eq.Fragments_Bond .or.
     &                lBondType.eq.Fragments_Bond) delta=0.0D0
                  If (Fi3.gt.Pi-delta) Go To 401
                  If (Fi3.lt.delta)    Go To 401
*
*---------------- 2-4-3
*
                  call dcopy_(3,Ref(1,2),1,RX4Y(1,1),1)
                  Call Bend(RX4Y,mCent,Fi4,Grad_ref,
     &                      .False.,
     &                      .False.,'        ',Hess,.False.)
#ifdef _DEBUG_
                  Write (6,*) '2-4-3: Fi4=',Fi4
#endif
                  delta = delta0
                  If (kBondType.eq.Fragments_Bond .or.
     &                lBondType.eq.Fragments_Bond) delta=0.0D0
                  If (Fi4.gt.Pi-delta) Go To 401
                  If (Fi4.lt.delta)    Go To 401
*
                  Call OutofP(Ref,nCent,Val,Grad,.False.,
     &                       .False.,
     &                       '        ',Hess,.False.)
#ifdef _DEBUG_
                  Write (6,*) 'Val=',Val*180.D0/Pi
#endif
*
                  If (Abs(Val).gt.35.D0*(Pi/180.D0))  Go To 401
*
                  Call OutofP(A,nCent,Val,Grad,.False.,.False.,
     &                      '        ',Hess,Proc_dB)
*
                  nq = nq + 1
                  If (.Not.Process) mB_Tot = mB_Tot + mB
                  If (.Not.Proc_dB) mdB_Tot = mdB_Tot + mB**2
#ifdef _DEBUG_
                  Write (6,*) 'nq=',nq
#endif
*
                  nqO = nqO + 1
                  iF1=1
                  Call NxtWrd(Name(iAtom),iF1,iE1)
                  Lbls(1)=Name(iAtom)(iF1:iE1)
*
                  iF2=1
                  Call NxtWrd(Name(jAtom),iF2,iE2)
                  Lbls(2)=Name(jAtom)(iF2:iE2)
                  If (kDCRR.ne.0) Then
                     Lbls(2)(iE2+1:iE2+1)='('
                     Lbls(2)(iE2+2:iE2+1+iChOp(kDCRR))=
     &                     ChOp(kDCRR)(1:iChOp(kDCRR))
                     Lbls(2)(iE2+2+iChOp(kDCRR):
     &                       iE2+2+iChOp(kDCRR))=')'
                     Call NxtWrd(Lbls(2),iF2,iE2)
                  End If
*
                  iF3=1
                  Call NxtWrd(Name(kAtom),iF3,iE3)
                  Lbls(3)=Name(kAtom)(iF3:iE3)
                  If (kDCRT.ne.0) Then
                     Lbls(3)(iE3+1:iE3+1)='('
                     Lbls(3)(iE3+2:iE3+1+iChOp(kDCRT))=
     &                     ChOp(kDCRT)(1:iChOp(kDCRT))
                     Lbls(3)(iE3+2+iChOp(kDCRT):
     &                       iE3+2+iChOp(kDCRT))=')'
                     Call NxtWrd(Lbls(3),iF3,iE3)
                  End If
*
                  iF4=1
                  Call NxtWrd(Name(lAtom),iF4,iE4)
                  Lbls(4)=Name(lAtom)(iF4:iE4)
                  If (kDCRTS.ne.0) Then
                     Lbls(4)(iE4+1:iE4+1)='('
                     Lbls(4)(iE4+2:iE4+1+iChOp(kDCRTS))=
     &                     ChOp(kDCRTS)(1:iChOp(kDCRTS))
                     Lbls(4)(iE4+2+iChOp(kDCRTS):
     &                       iE4+2+iChOp(kDCRTS))=')'
                     Call NxtWrd(Lbls(4),iF4,iE4)
                  End If
                  Write (LuIC,'(A,I3.3,8A)')
     &                   'o',nqO,' = Outofp   ',
     &                       Lbls(2)(iF2:iE2),
     &                   ' ',Lbls(3)(iF3:iE3),
     &                   ' ',Lbls(4)(iF4:iE4),
     &                   ' ',Lbls(1)(iF1:iE1)
#ifdef _DEBUG_
                  If (iPrint.ge.49)
     &            Write (6,'(A,I3.3,8A)')
     &                   'o',nqO,' = Outofp   ',
     &                       Lbls(2)(iF2:iE2),
     &                   ' ',Lbls(3)(iF3:iE3),
     &                   ' ',Lbls(4)(iF4:iE4),
     &                   ' ',Lbls(1)(iF1:iE1)
                  Write (6,*) 'iDeg=',iDeg
#endif
                  Label=' '
                  Write (Label,'(A,I3.3)') 'o',nqO
*
*
                  If (Process) Then
*
                     Indq(1,nq)=6
                     ij = (jAtom-1)*nAtoms + iAtom
                     kl = (lAtom-1)*nAtoms + kAtom
                     Indq(2,nq) = (kl-1)*nAtoms**2 + ij
                     ijDCR = kDCRT*8 + kDCRR+1
                     Indq(3,nq) = kDCRS*8**2 + ijDCR
*
*                    f_Const=Max(f_Const,f_Const_Min)
#ifdef _DEBUG_
                     Write (6,*) 'f_const=',f_const
#endif
                     fconst(nq)=Sqrt(f_Const)
                     rMult(nq)=Deg
*
                     Value(nq,iIter)=Val
                     qLbl(nq) = Label
*
*------------------- Project the gradient vector
*
                     Call ProjSym(nAtoms,nCent,Ind,nStab,
     &                            jStab,A,iDCR,Grad,
     &                            Smmtrc,nDim,PSPrint,Hess,
     &                            mB_Tot,mdB_Tot,
     &                            BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,
     &                            Proc_dB,nqB,nB,nq,rMult(nq))
*
                  End If
*
  401             Continue
               End Do            ! lNeighbor
  301          Continue
            End Do               ! kNeighbor
  201       Continue
         End Do                  ! iCase
  101    Continue
      End Do                     ! jBond
*
      Call QExit ('OutOfPs')
      Return
      End
