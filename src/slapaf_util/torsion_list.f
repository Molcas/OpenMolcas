************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Torsion_List(
     &                 nq,
     &                 nAtoms,iIter,nIter,Cx,jStab,
     &                 nStab,Smmtrc,Process,Value,
     &                 nB,iANr,qLbl,iRef,
     &                 fconst,rMult,LuIC,Name,Indq,iPrv,Proc_dB,
     &                 iTabBonds,nBonds,iTabAI,mAtoms,iTabAtoms,nMax,
     &                 mB_Tot,mdB_Tot,
     &                 BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,nqB)
      use Symmetry_Info, only: nIrrep, iOper
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Parameter (mB = 4*3)
      Real*8 Cx(3,nAtoms,nIter), A(3,4), Grad(mB), Hess(mB**2),
     &       fconst(nB), Value(nB,nIter),
     &       Ref(3,4), Prv(3,4), rMult(nB), Grad_ref(9),
     &       BM(nB_Tot), dBM(ndB_Tot)
      Integer   nStab(nAtoms), iANr(nAtoms),
     &          iDCRR(0:7), jStab(0:7,nAtoms),
     &          iStabM(0:7), Ind(4), iDCR(4), iDCRT(0:7),
     &          iDCRS(0:7), iStabN(0:7), iStabO(0:7), iChOp(0:7),
     &          Indq(3,nB), iDCRX(0:7), iDCRY(0:7), nqB(nB),
     &          iTabBonds(3,nBonds), iTabAI(2,mAtoms),
     &          iTabAtoms(2,0:nMax,mAtoms), iBM(nB_Tot),idBM(2,ndB_Tot)
      Logical Smmtrc(3,nAtoms), Process,
     &        MinBas, Help, Proc_dB, R_Stab_A, Torsion_Check
      Character*14 Label, qLbl(nB)
      Character*3 ChOp(0:7)
#include "Molcas.fh"
      Character*(LENIN) Name(nAtoms)
      Character*(LENIN4) Lbls(4)
#include "bondtypes.fh"
#define _FMIN_
#include "ddvdt.fh"
#include "ddvdt_trsn.fh"
      Data ChOp/'E  ','X  ','Y  ','XY ','Z  ','XZ ','YZ ','XYZ'/
      Data iChOp/1,1,1,2,1,2,2,3/
      Data f_Const_Min2/1.0D-1/
#include "constants.fh"
*                                                                      *
************************************************************************
*                                                                      *
!#define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
      If (nBonds.lt.3) Return

      nqT=0
      Call FZero(Hess,144)
*
*---- Loop over dihedrals
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
*     Order will play a role here. That is, the torsion
*     A-R(B)-T(C)-TS(D) is NOT identical to
*     A-R(B)-TS(D)-T(C). Hence we put no restriction on the
*     pairs AB and CD. However, for the pair of pairs we have
*     that order is irrelevant, i.e. ABCD is identical to
*     DCBA. To guarantee this we limit the pairs to the unique
*     combinations.
*
*     Start with the center bond: B-C
*
#ifdef _DEBUGPRINT_
      Write (6,*)
      Write (6,*) "List of all available bonds"
      Write (6,*)
      Do iBond = 1, nBonds
         iBondType=iTabBonds(3,iBond)
         Write (6,*)
         jAtom_ = iTabBonds(1,iBond)
         kAtom_ = iTabBonds(2,iBond)
         jAtom = iTabAI(1,jAtom_)
         kAtom = iTabAI(1,kAtom_)
         Write (6,*) 'Atoms-pair:  ',Name(jAtom),Name(kAtom)
         Write (6,*) 'iBond,iBondType=',iBond,Bondtype(Min(3,iBondType))
      End Do
#endif
      Do iBond = 1, nBonds
*
*        The center bond may be a "magic" bond
         iBondType=iTabBonds(3,iBond)
#ifdef _DEBUGPRINT_
         Write (6,*)
         jAtom_ = iTabBonds(1,iBond)
         kAtom_ = iTabBonds(2,iBond)
         jAtom = iTabAI(1,jAtom_)
         kAtom = iTabAI(1,kAtom_)
         Write (6,*) 'Atoms-pair:',Name(jAtom),Name(kAtom)
         Write (6,*) 'iBond,iBondType=',iBond,Bondtype(Min(3,iBondType))
#endif
*
*        Center bond should not be a van der Waals bond,
*        anything else goes!
*
         If (iBondType.eq.vdW_Bond) Go To 201
*
*        Extract index to the center atom in an "Magic" bond if this
*        is a magic bond.
*
         If (iBondType.gt.Magic_Bond) Then
            iMagic=iBondType-3
         Else
            iMagic=0
         End If
*
*        cases: BC or CB
*
         Do iCase = 1, 2
*
            If (iCase.eq.1) Then
               jAtom_ = iTabBonds(1,iBond)
               kAtom_ = iTabBonds(2,iBond)
            Else
               jAtom_ = iTabBonds(2,iBond)
               kAtom_ = iTabBonds(1,iBond)
            End If
*
            jAtom = iTabAI(1,jAtom_)
            kAtom = iTabAI(1,kAtom_)
            jr = iTabRow(iANr(jAtom))
            kr = iTabRow(iANr(kAtom))
            Ind(2) = jAtom
            Ind(3) = kAtom
            iDCR(2) = iTabAI(2,jAtom_)
            If (R_Stab_A(iDCR(2),jStab(0,jAtom),nStab(jAtom)))
     &          iDCR(2)=iOper(0)
            iDCR(3) = iTabAI(2,kAtom_)
            If (R_Stab_A(iDCR(3),jStab(0,kAtom),nStab(kAtom)))
     &          iDCR(3)=iOper(0)
#ifdef _DEBUGPRINT_
            Write (6,*)
            Write (6,*) 'R,T=',Name(jAtom),ChOp(iDCR(2)),
     &                         Name(kAtom),ChOp(iDCR(3))
            Write (6,*)
#endif
*
            nNeighbor_j = iTabAtoms(1,0,jAtom_)
            nCoBond_j=nCoBond(jAtom_,mAtoms,nMax,iTabBonds,
     &                        nBonds,iTabAtoms)
            nFgBond_j=nFgBond(jAtom_,mAtoms,nMax,iTabBonds,
     &                        nBonds,iTabAtoms)
            nNeighbor_k = iTabAtoms(1,0,kAtom_)
            nCoBond_k=nCoBond(kAtom_,mAtoms,nMax,iTabBonds,
     &                        nBonds,iTabAtoms)
            nFgBond_k=nFgBond(kAtom_,mAtoms,nMax,iTabBonds,
     &                        nBonds,iTabAtoms)
            If ( nCoBond_j.lt.2.and.nFgBond_j.eq.0) Go To 250
            If ( nCoBond_k.lt.2.and.nFgBond_k.eq.0) Go To 250
*
#ifdef _DEBUGPRINT_
            Write (6,*) 'nNeighbor_j,nNeighbor_k=',
     &                   nNeighbor_j,nNeighbor_k
            Write (6,*)
#endif
*
         Do iNeighbor = 1, nNeighbor_j
            iAtom_ = iTabAtoms(1,iNeighbor,jAtom_)
            nCoBond_i=nCoBond(iAtom_,mAtoms,nMax,iTabBonds,
     &                        nBonds,iTabAtoms)
            nFgBond_i=nFgBond(iAtom_,mAtoms,nMax,iTabBonds,
     &                        nBonds,iTabAtoms)
            jBond  = iTabAtoms(2,iNeighbor,jAtom_)
            If (jBond.eq.iBond) Go To 301
            jBondType=iTabBonds(3,jBond)
#ifdef _DEBUGPRINT_
            Write (6,*) 'jBond,jBondType=',jBond,
     &                  BondType(Min(3,jBondType))
#endif
            If (jBondType.eq.vdW_Bond   .or.
     &          jBondType.eq.Magic_Bond) Go To 301
*           If (nCoBond_j.gt.2.and.
*    &          (nCoBond_i.ge.4.and.nFgBond_i.eq.0)) Go To 301
*           If (nCoBond_i.ge.8 .and.
*    &          nCoBond_j.ge.8 .and.
*    &          nCoBond_k.ge.8      ) Go To 301
            iAtom = iTabAI(1,iAtom_)
            If (iBondType.gt.Magic_Bond .and. iAtom.eq.iMagic) Go To 301
            ir = iTabRow(iANr(iAtom))
            Ind(1) = iAtom
            iDCR(1)= iTabAI(2,iAtom_)
            If (R_Stab_A(iDCR(1),jStab(0,iAtom),nStab(iAtom)))
     &          iDCR(1)=iOper(0)
*
*---------- Torsion should be A-..., eliminate P(A)-...
*
            If (iDCR(1).ne.iOper(0)) Go To 301
*
*----------- Eliminate A-R(B)-T(C)-TS(C) over A-B-RT(C)-RTS(D)
*            Proceed if A-R-T(C)-TS(C)
*
            If (R_Stab_A(iDCR(2),jStab(0,iAtom),nStab(iAtom)) .and.
     &          iDCR(2).ne.iOper(0)) Go To 301
*
*----------- Eliminate A-R(B)-T(C)-TS(D) over A-R(B)-C-S(D)
*            Proceed if A-R(B)-C-S(D)
*
            If (R_Stab_A(iDCR(3),jStab(0,iAtom),nStab(iAtom)) .and.
     &          R_Stab_A(iDCR(3),jStab(0,jAtom),nStab(jAtom)) .and.
     &          iDCR(3).ne.iOper(0)) Go To 301
#ifdef _DEBUGPRINT_
            Write (6,*)
            Write (6,*) 'E=',Name(iAtom),ChOp(iDCR(1))
#endif
            call dcopy_(3,Cx(1,iAtom,iIter),1,A,  1)
            call dcopy_(3,Cx(1,iAtom,iRef), 1,Ref,1)
            call dcopy_(3,Cx(1,iAtom,iPrv), 1,Prv,1)
            Help = ir.gt.3.or.jr.gt.3
*
*---------- Form double coset representatives for (iAtom,jAtom)
*
            Call DCR(Lambda,
     &               jStab(0,iAtom),nStab(iAtom),
     &               jStab(0,jAtom),nStab(jAtom),
     &               iDCRR,nDCRR)
            kDCRR=iDCR(2)
#ifdef _DEBUGPRINT_
            Write (6,'(10A)') 'R={',(ChOp(iDCRR(i)),
     &                            i=0,nDCRR-1),'}  '
            Write (6,'(2A)') 'R=',ChOp(kDCRR)
#endif
            Call OA(kDCRR,Cx(1:3,jAtom,iIter),  A(1:3,2))
            Call OA(kDCRR,Cx(1:3,jAtom,iRef ),Ref(1:3,2))
            Call OA(kDCRR,Cx(1:3,jAtom,iPrv ),Prv(1:3,2))
#ifdef _DEBUGPRINT_
            Write (6,'(10A)') 'U={',(ChOp(jStab(i,iAtom)),
     &                         i=0,nStab(iAtom)-1),'}  '
            Write (6,'(10A)') 'V={',(ChOp(jStab(i,jAtom)),
     &                         i=0,nStab(jAtom)-1),'}  '
#endif
*
*---------- Form stabilizer for (iAtom,jAtom)
*
            Call Inter(jStab(0,iAtom),nStab(iAtom),
     &                 jStab(0,jAtom),nStab(jAtom),
     &                 iStabM,nStabM)
#ifdef _DEBUGPRINT_
            Write (6,'(10A)') 'M={',
     &            (ChOp(iStabM(i)),i=0,nStabM-1),'}  '
#endif
            If (Help) Then
               rij2=(Ref(1,1)-Ref(1,2))**2
     &             +(Ref(2,1)-Ref(2,2))**2
     &             +(Ref(3,1)-Ref(3,2))**2
               Rab=Sqrt(rij2)
               RabCov=(CovRadT(iANr(iAtom))
     &               +CovRadT(iANr(jAtom)))/bohr
               r0=Zero
               Alpha=Zero
               f_Const_ij=Zero
               f_Const_ij_Ref=f_Const_ij
            Else
               r0=rAV(ir,jr)
               Alpha=aAv(ir,jr)
               rij2_Ref=(Ref(1,1)-Ref(1,2))**2
     &                 +(Ref(2,1)-Ref(2,2))**2
     &                 +(Ref(3,1)-Ref(3,2))**2
               f_Const_ij_Ref=rkt*Exp(Alpha*(r0**2-rij2_Ref))
               rij2=(A(1,1)-A(1,2))**2
     &             +(A(2,1)-A(2,2))**2
     &             +(A(3,1)-A(3,2))**2
               f_Const_ij=rkt*Exp(Alpha*(r0**2-rij2))
            End If
*
            Do lNeighbor = 1, nNeighbor_k
               lAtom_ = iTabAtoms(1,lNeighbor,kAtom_)
               nCoBond_l=nCoBond(lAtom_,mAtoms,nMax,iTabBonds,
     &                           nBonds,iTabAtoms)
               nFgBond_l=nFgBond(lAtom_,mAtoms,nMax,iTabBonds,
     &                           nBonds,iTabAtoms)
               lAtom = iTabAI(1,lAtom_)
               Ind(4) = lAtom
               iDCR(4) = iTabAI(2,lAtom_)
               If (R_Stab_A(iDCR(4),jStab(0,lAtom),nStab(lAtom)))
     &             iDCR(4)=iOper(0)
               kBond  = iTabAtoms(2,lNeighbor,kAtom_)
               If (kBond.eq.iBond) Go To 401
               If (lAtom_.eq.iAtom_) Go To 401
               kBondType=iTabBonds(3,kBond)
#ifdef _DEBUGPRINT_
         Write (6,*) 'kBond,kBondType=',kBond,Bondtype(Min(3,kBondType))
#endif
               If (kBondType.eq.vdW_Bond    .or.
     &             kBondType.eq.Magic_Bond) Go To 401
*              If (nCoBond_k.gt.2.and.
*    &             (nCoBond_l.ge.4.and.nFgBond_l.eq.0)) Go To 401
*              If (nCoBond_j.ge.8 .and.
*    &             nCoBond_k.ge.8 .and.
*    &             nCoBond_l.ge.8      ) Go To 401
*
               If (iBondType.gt.Magic_Bond .and. lAtom.eq.iMagic)
     &            Go To 401
               lr = iTabRow(iANr(lAtom))
               kDCRT= iDCR(3)
               kDCRTS=iDCR(4)
               kDCRS=iEor(kDCRTS,kDCRT)
#ifdef _DEBUGPRINT_
               Write (6,*) 'i,j,k,l=',
     &               Name(iAtom),ChOp(iDCR(1)),
     &               Name(jAtom),ChOp(iDCR(2)),
     &               Name(kAtom),ChOp(iDCR(3)),
     &               Name(lAtom),ChOp(iDCR(4))
#endif
*
*------------- Eliminate A-R(B)-T(C)-TS(D) over A-TSR(B)-S(C)-D
*
               If (R_Stab_A(iDCR(4),jStab(0,iAtom),nStab(iAtom)) .and.
     &             R_Stab_A(iDCR(4),jStab(0,jAtom),nStab(jAtom)) .and.
     &             R_Stab_A(iDCR(4),jStab(0,kAtom),nStab(kAtom)) .and.
     &             iDCR(4).ne.iOper(0)) Go To 401
*
               nE=1
               If (iDCR(2).eq.iOper(0)) nE=nE+1
               If (iDCR(3).eq.iOper(0)) nE=nE+1
               If (iDCR(4).eq.iOper(0)) nE=nE+1
               mE=1
               If (R_Stab_A(iDCR(4),
     &             jStab(0,iAtom),nStab(iAtom))) mE=mE+1
               If (R_Stab_A(iEOr(iDCR(4),iDCR(3)),
     &             jStab(0,kAtom),nStab(kAtom))) mE=mE+1
               If (R_Stab_A(iEOr(iDCR(4),iDCR(2)),
     &             jStab(0,jAtom),nStab(jAtom))) mE=mE+1
               If (nE.lt.mE) Go To 401
               If (nE.eq.mE .and. iAtom.gt.lAtom) Go To 401

#ifdef _DEBUGPRINT_
               Write (6,*)
               Write (6,*) 'TS=',Name(lAtom),ChOp(iDCR(4))
#endif
*
               Help = ir.gt.3.or.jr.gt.3.or.kr.gt.3.or.lr.gt.3
*
               Write (Label,'(A,I2,A,I2,A,I2,A,I2,A)')
     &                'D(',iAtom,',',jAtom,',',kAtom,',',lAtom,')'
#ifdef _DEBUGPRINT_
               Write (6,'(A,I2,A,I2,A,I2,A,I2,A)')
     &                'D(',iAtom,',',jAtom,',',kAtom,',',lAtom,')'
#endif
*
*------------- Form double coset representatives for (kAtom,lAtom)
*
               Call DCR(Lambda,
     &                  jStab(0,kAtom),nStab(kAtom),
     &                  jStab(0,lAtom),nStab(lAtom),
     &                  iDCRS,nDCRS)
#ifdef _DEBUGPRINT_
               Write (6,'(10A)') 'S={',(ChOp(iDCRS(i)),
     &                            i=0,nDCRS-1),'}  '
               Write (6,'(10A)') 'W={',(ChOp(jStab(i,kAtom)),
     &                            i=0,nStab(kAtom)-1),'}  '
               Write (6,'(10A)') 'X={',(ChOp(jStab(i,lAtom)),
     &                            i=0,nStab(lAtom)-1),'}  '
               Write (6,'(2A)') 'S=',ChOp(kDCRS)
#endif
*
               Ref(1:3,3) = Cx(1:3,kAtom,iRef)
               Call OA(kDCRS,Cx(1:3,lAtom,iRef),Ref(1:3,4))
               Prv(1:3,3) = Cx(1:3,kAtom,iPrv)
               Call OA(kDCRS,Cx(1:3,lAtom,iPrv),Prv(1:3,4))
*
               If (Help) Then
                  rkl2=(Ref(1,3)-Ref(1,4))**2
     &                +(Ref(2,3)-Ref(2,4))**2
     &                +(Ref(3,3)-Ref(3,4))**2
                  Rcd=Sqrt(rkl2)
                  RcdCov=(CovRadT(iANr(kAtom))
     &                  +CovRadT(iANr(lAtom)))/bohr
               Else
                  r0=rAV(kr,lr)
                  rkl2=(Ref(1,3)-Ref(1,4))**2
     &                +(Ref(2,3)-Ref(2,4))**2
     &                +(Ref(3,3)-Ref(3,4))**2
               End If
*
*------------- Form stabilizer for (kAtom,lAtom)
*
               Call Inter(jStab(0,kAtom),nStab(kAtom),
     &                    jStab(0,lAtom),nStab(lAtom),
     &                    iStabN,nStabN)
*
#ifdef _DEBUGPRINT_
               Write (6,'(10A)') 'N={',
     &               (ChOp(iStabN(i)),i=0,nStabN-1),'}  '
#endif
*
*------------- Form double coset representatives for
*              ((iAtom,jAtom),(kAtom,lAtom))
*
               Call DCR(Lambda,
     &                  iSTabM,nStabM,
     &                  iStabN,nStabN,
     &                  iDCRT,nDCRT)
*
*------------- Take care of some special cases which normally
*              are not included. If A=B we will normally exclude
*              the pairs R(A)-A and TS(C)-T(C).
*
               Call iCopy(nDCRT,iDCRT,1,iDCRX,1)
               Call iCopy(nDCRT,iDCRT,1,iDCRY,1)
               nDCRX=nDCRT
               nDCRY=nDCRT
               If (iAtom.eq.jAtom) Then
*                 Write (6,*) ' Special fix'
                  Call Union(iDCRX,nDCRX,iDCRY,nDCRY,
     &                       kDCRR,iDCRT,nDCRT)
               Else If (kAtom.eq.lAtom) Then
*                 Write (6,*) ' Special fix'
                  Call Union(iDCRX,nDCRX,iDCRY,nDCRY,
     &                       kDCRS,iDCRT,nDCRT)
               End If
*
#ifdef _DEBUGPRINT_
               Write (6,'(10A)') 'T={',
     &               (ChOp(iDCRT(i)),i=0,nDCRT-1),'}  '
               Write (6,'(2A)') 'kDCRT=',ChOp(kDCRT)
               Write (6,'(2A)') 'T=',ChOp(kDCRT)
#endif
*
               kDCRTS=iEor(kDCRT,kDCRS)
*
               Call OA(kDCRT ,Cx(1:3,kAtom,iIter),  A(1:3,3))
               Call OA(kDCRT ,Cx(1:3,kAtom,iRef ),Ref(1:3,3))
               Call OA(kDCRT ,Cx(1:3,kAtom,iPrv ),Prv(1:3,3))
               Call OA(kDCRTS,Cx(1:3,lAtom,iIter),  A(1:3,4))
               Call OA(kDCRTS,Cx(1:3,lAtom,iRef ),Ref(1:3,4))
               Call OA(kDCRTS,Cx(1:3,lAtom,iPrv ),Prv(1:3,4))
*
*------------- Form the stabilizer for the torsion
*
               If (iAtom.eq.lAtom .and.
     &             jAtom.eq.kAtom.and.kDCRR.eq.kDCRS) Then
                  Call Union(iStabM,nStabM,
     &                       iStabN,nStabN,
     &                       kDCRTS,iStabO,nStabO)
               Else
                  Call Inter(iStabM,nStabM,
     &                       iStabN,nStabN,
     &                       iStabO,nStabO)
               End If
*
#ifdef _DEBUGPRINT_
               Write (6,'(10A)') 'M={',
     &               (ChOp(iStabM(i)),i=0,nStabM-1),'}  '
               Write (6,'(10A)') 'N={',
     &               (ChOp(iStabN(i)),i=0,nStabN-1),'}  '
               Write (6,'(10A)') 'O={',
     &               (ChOp(iStabO(i)),i=0,nStabO-1),'}  '
#endif
*
*------------- Compute the degeneracy of the torsion
*
               iDeg=nIrrep/nStabO
               Deg=Sqrt(DBLE(iDeg))
*
*------------- Test if coordinate should be included
*
               If (Help) Then
                  rjk2=(Ref(1,2)-Ref(1,3))**2
     &                +(Ref(2,2)-Ref(2,3))**2
     &                +(Ref(3,2)-Ref(3,3))**2
                  Rbc=Sqrt(rjk2)
                  RbcCov=(CovRadT(iANr(jAtom))
     &                  +CovRadT(iANr(kAtom)))/bohr
                  Diff=RbcCov-Rbc
                  If (Diff.lt.Zero) Diff = Zero
                  f_Const=A_Trsn(1)+A_Trsn(2)*Diff
                  f_Const=f_Const*Fact
                  r0=Zero
                  Alpha=Zero
                  f_Const_Ref=f_Const
               Else
                  r0=rAV(jr,kr)
                  Alpha=aAv(jr,kr)
                  rjk2_Ref=(Ref(1,2)-Ref(1,3))**2
     &                    +(Ref(2,2)-Ref(2,3))**2
     &                    +(Ref(3,2)-Ref(3,3))**2
                  f_Const_ijk_Ref=f_Const_ij_Ref
     &                       *Exp(Alpha*(r0**2-rjk2_Ref))
                  rjk2=(A(1,2)-A(1,3))**2
     &                +(A(2,2)-A(2,3))**2
     &                +(A(3,2)-A(3,3))**2
                  f_Const_ijk=f_Const_ij
     &                       *Exp(Alpha*(r0**2-rjk2))
*
                  r0=rAV(kr,lr)
                  Alpha=aAv(kr,lr)
                  rkl2_Ref=(Ref(1,3)-Ref(1,4))**2
     &                    +(Ref(2,3)-Ref(2,4))**2
     &                    +(Ref(3,3)-Ref(3,4))**2
                  f_Const_Ref=f_Const_ijk_Ref
     &                   *Exp(Alpha*(r0**2-rkl2_Ref))
                  rkl2=(A(1,3)-A(1,4))**2
     &                +(A(2,3)-A(2,4))**2
     &                +(A(3,3)-A(3,4))**2
                  f_Const=f_Const_ijk
     &                   *Exp(Alpha*(r0**2-rkl2))
               End If
               If (Torsion_Check(iAtom,jAtom,kAtom,lAtom,
     &                           Ref,iTabAtoms,
     &                           nMax,mAtoms)) Then
                  f_Const_Ref=Max(f_Const_Ref,10.D0*F_Const_Min)
               End If

               If (f_Const_Ref.lt.f_Const_Min .and.
     &             iBondType.ne.Fragments_Bond .and.
     &             iBondType.le.Fragments_Bond .and.
     &             jBondType.ne.Fragments_Bond .and.
     &             kBondType.ne.Fragments_Bond
     &            )  Go To 401
*
*
*------------- Check that valence angles are above threshold
*
               mCent=3
               delta = (15.0D0/180.D0)*Pi
               If (nAtoms.eq.4) delta = -Ten
               Call Bend(Ref(1,1),mCent,Fi2,Grad_ref,.False.,
     &                   .False.,'        ',Hess,.False.)
               If (Fi2.gt.Pi-delta) Go To 401
               If (Fi2.lt.delta) Go To 401
               Call Bend(Ref(1,2),mCent,Fi3,Grad_ref,.False.,
     &                   .False.,'        ',Hess,.False.)
               If (Fi3.gt.Pi-delta) Go To 401
               If (Fi3.lt.delta) Go To 401
*              Write (6,*) ' T Force Constant:',f_Const
*
               nq = nq + 1
               If (.Not.Process) mB_Tot = mB_Tot + mB
               If (.Not.Proc_dB) mdB_Tot = mdB_Tot + mB**2
*              Write (6,*) 'nq=',nq
*
               nqT = nqT + 1
               iF1=1
               Call NxtWrd(Name(iAtom),iF1,iE1)
               Lbls(1)=Name(iAtom)(iF1:iE1)
               iF2=1
               Call NxtWrd(Name(jAtom),iF2,iE2)
               Lbls(2)=Name(jAtom)(iF2:iE2)
               If (kDCRR.ne.0) Then
                  Lbls(2)(iE2+1:iE2+1)='('
                  Lbls(2)(iE2+2:iE2+1+iChOp(kDCRR))=
     &                  ChOp(kDCRR)(1:iChOp(kDCRR))
                  Lbls(2)(iE2+2+iChOp(kDCRR):
     &                    iE2+2+iChOp(kDCRR))=')'
                  Call NxtWrd(Lbls(2),iF2,iE2)
               End If
               iF3=1
               Call NxtWrd(Name(kAtom),iF3,iE3)
               Lbls(3)=Name(kAtom)(iF3:iE3)
               If (kDCRT.ne.0) Then
                  Lbls(3)(iE3+1:iE3+1)='('
                  Lbls(3)(iE3+2:iE3+1+iChOp(kDCRT))=
     &                  ChOp(kDCRT)(1:iChOp(kDCRT))
                  Lbls(3)(iE3+2+iChOp(kDCRT):
     &                    iE3+2+iChOp(kDCRT))=')'
                  Call NxtWrd(Lbls(3),iF3,iE3)
               End If
               iF4=1
               Call NxtWrd(Name(lAtom),iF4,iE4)
               Lbls(4)=Name(lAtom)(iF4:iE4)
               If (kDCRTS.ne.0) Then
                  Lbls(4)(iE4+1:iE4+1)='('
                  Lbls(4)(iE4+2:iE4+1+iChOp(kDCRTS))=
     &                  ChOp(kDCRTS)(1:iChOp(kDCRTS))
                  Lbls(4)(iE4+2+iChOp(kDCRTS):
     &                    iE4+2+iChOp(kDCRTS))=')'
                  Call NxtWrd(Lbls(4),iF4,iE4)
               End If
               Write (LuIC,'(A,I3.3,8A)')
     &                't',nqT,' = Dihedral ',
     &                    Lbls(1)(iF1:iE1),
     &                ' ',Lbls(2)(iF2:iE2),
     &                ' ',Lbls(3)(iF3:iE3),
     &                ' ',Lbls(4)(iF4:iE4)
#ifdef _DEBUGPRINT_
               Write (6,'(A,I3.3,8A)')
     &                't',nqT,' = Dihedral ',
     &                    Lbls(1)(iF1:iE1),
     &                ' ',Lbls(2)(iF2:iE2),
     &                ' ',Lbls(3)(iF3:iE3),
     &                ' ',Lbls(4)(iF4:iE4)
               Write (6,*) 'iDeg=',iDeg
#endif
               Label=' '
               Write (Label,'(A,I3.3)') 't',nqT
*
               Call Trsn(A,nCent,Val,Grad,.False.,.False.,
     &                   '        ',Hess,Proc_dB)
               If (iIter.eq.iPrv) Then
                  Val_Prv=Val
               Else
                  Val_Prv=Value(nq,iPrv)
               End If
*
*------------- correct for 2Pi flip relative to reference.
*
               Range1=Two*Pi*0.8D0
               Range2=Pi*0.80D0
               Range3=Pi*1.20D0
               If (Abs(Val_Prv-Val).gt.Range1) Then
                  If (Sign(One,Val_Prv).eq.One) Then
                     Val=Val+Two*Pi
                  Else
                     Val=Val-Two*Pi
                  End If
               Else If (Abs(Val_Prv-Val).gt.Range2 .and.
     &                  Abs(Val_Prv-Val).lt.Range3) Then
                  If (Val_Prv-Val.gt.Zero) Then
                     Val=Val+Pi
                  Else
                     Val=Val-Pi
                  End If
               End If
#ifdef _DEBUGPRINT_
               Call RecPrt('Trsns:  B',' ',Grad,3,4)
               Call RecPrt('Trsns: dB',' ',Hess,12,12)
#endif
*
               If (Process) Then
*
                  Indq(1,nq)=5
                  ij = (jAtom-1)*nAtoms + iAtom
                  kl = (lAtom-1)*nAtoms + kAtom
                  Indq(2,nq) = (kl-1)*nAtoms**2 + ij
                  ijDCR = kDCRT*8 + kDCRR+1
                  Indq(3,nq) = kDCRS*8**2 + ijDCR
*
                  If (iMagic.ne.0) Then
                     f_Const=Max(f_Const,f_Const_Min2)
                  Else If (.Not. Help .and.
     &                    (iBondType.eq.Fragments_Bond .or.
     &                     jBondType.eq.Fragments_Bond .or.
     &                     kBondType.eq.Fragments_Bond )) Then
                     f_Const=Max(f_Const,f_Const_Min*1.D3)
                  Else
                     f_Const=Max(f_Const,f_Const_Min)
                  End If
                  fconst(nq)=Sqrt(f_Const)
                  rMult(nq)=Deg
*                 Scale down fconst if angles are close to linear
                  CosFi=Max(Abs(Cos(Fi2)),Abs(Cos(Fi3)))
                  CosThr=0.97D0
                  If (CosFi.gt.CosThr) Then
                     CosFact=(CosFi-CosThr)/(One-CosThr)
                     CosFact=One-(One-1.0D-3)*CosFact**2
                     fconst(nq)=CosFact*fconst(nq)
                  End If
*
                  Value(nq,iIter)=Val
                  qLbl(nq) = Label
*
*---------------- Project the gradient vector
*
                  Call ProjSym(nAtoms,nCent,Ind,nStab,
     &                         jStab,A,iDCR,Grad,
     &                         Smmtrc,Hess,
     &                         mB_Tot,mdB_Tot,
     &                         BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,
     &                         Proc_dB,nqB,nB,nq,rMult(nq))
*
               End If
*
  401          Continue
            End Do               ! iNeighbor_k
  301       Continue
         End Do                  ! iNeighbor_j
 250     Continue
         End Do                  ! iCase
  201    Continue
      End Do                     ! iBonds
#ifdef _DEBUGPRINT_
      Write (6,*) 'nqT=',nqT
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
      Function Torsion_Check(iAtom,jAtom,kAtom,lAtom,Ref,iTabAtoms,
     &                       nMax,mAtoms)
      Implicit Real*8 (a-h,o-z)
      Logical Torsion_Check
      Real*8 Ref(3,4), Coor(3,4)
      Integer iTabAtoms(2,0:nMax,mAtoms), iCase(4,12), iSet(4), jSet(4)
      Data iCase/1,2,3,4, 1,2,4,3, 2,1,3,4, 2,1,4,3,
     &           1,3,2,4, 1,3,4,2, 3,1,2,4, 3,1,4,2,
     &           1,4,2,3, 1,4,3,2, 4,1,2,3, 4,1,3,2/
*
      Torsion_Check=.False.
      iSet(1)=iAtom
      iSet(2)=jAtom
      iSet(3)=kAtom
      iSet(4)=lAtom
C     Write (*,*) 'Here we go!'
C     Write (*,*) 'iSet=',iSet
*
*     Check all twelve possible torsions to see which one is the one
*     which should at least be included.
*
      FC=0.0D0
      Test=0.0D0
      Do ii = 1, 12
         jSet(1) = iSet(iCase(1,ii))
         jSet(2) = iSet(iCase(2,ii))
         jSet(3) = iSet(iCase(3,ii))
         jSet(4) = iSet(iCase(4,ii))
C        Write (*,*) 'ii,jSet=',ii,jSet
         call dcopy_(3,Ref(1,iCase(1,ii)),1,Coor(1,1),1)
         call dcopy_(3,Ref(1,iCase(2,ii)),1,Coor(1,2),1)
         call dcopy_(3,Ref(1,iCase(3,ii)),1,Coor(1,3),1)
         call dcopy_(3,Ref(1,iCase(4,ii)),1,Coor(1,4),1)
C        Call RecPrt('Coor',' ',Coor,3,4)
         Test=FC_Torsion(jSet,Coor,iTabAtoms,nMax,mAtoms)
         If (Test.gt.FC) Then
            If (ii.gt.1) Return
            FC=Test
         End If
      End Do
      Torsion_Check=.True.
C     Write (*,*) 'Torsion_Check=',Torsion_Check
*
      Return
      End
      Function FC_Torsion(jSet,Coor,iTabAtoms,nMax,mAtoms)
      Implicit Real*8 (a-h,o-z)
      Real*8 FC_Torsion
      Integer iTabAtoms(2,0:nMax,mAtoms), iCase(2,3), jSet(4)
      Real*8  Coor(3,4), Gij(3)
#include "bondtypes.fh"
      Data iCase/1,2, 2,3, 3,4/
*
      FC_Torsion=0.0D0
      Do ii = 1, 3
         iHit=0
         iAtom=jSet(iCase(1,ii))
         jAtom=jSet(iCase(2,ii))
*
*        Loop over the neighbor atoms of iAtom to see if jAtom is a neighbor. If
*        it is check that the bond is classified as a covalent bond.
*
         Do i = 1, iTabAtoms(1,0,iAtom)
            If (jAtom.eq.iTabAtoms(1,i,iAtom) .and.
     &          iTabAtoms(2,i,iAtom).eq.Covalent_Bond) Then
                iHit=1
                iAtom_=iCase(1,ii)
                jAtom_=iCase(2,ii)
                Gij(ii)=1.0D0/Sqrt( (Coor(1,iAtom_)-Coor(1,jAtom_))**2
     &                             +(Coor(2,iAtom_)-Coor(2,jAtom_))**2
     &                             +(Coor(3,iAtom_)-Coor(3,jAtom_))**2 )
            End If
         End Do
         If (iHit.eq.0) Return
      End Do
      FC_Torsion=Gij(1)*Gij(2)*Gij(3)
*
      Return
      End
