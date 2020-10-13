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
      Subroutine Bond_List(
     &                 nq,nAtoms,iIter,nIter,Cx,jStab,
     &                 nStab,nDim,Smmtrc,Process,Value,
     &                 nB,iANr,qLbl,fconst,
     &                 rMult,iOptC,LuIC,Name,Indq,
     &                 Proc_dB,iTabBonds,nBonds,
     &                 iTabAI,mAtoms,mB_Tot,mdB_Tot,
     &                 BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,mqB)
      use Symmetry_Info, only: nIrrep, iOper
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
      Parameter (mB=2*3)
      Real*8 Cx(3,nAtoms,nIter), A(3,2), Grad(mB), Hess(mB**2),
     &       fconst(nB), Value(nB,nIter), rMult(nB),
     &       BM(nB_Tot), dBM(ndB_Tot)
      Integer   nStab(nAtoms), iDCRR(0:7), jStab(0:7,nAtoms),
     &          iStabM(0:7), Ind(2), iDCR(2), iANr(nAtoms), iChOp(0:7),
     &          Indq(3,nB), iTabBonds(3,nBonds), iTabAI(2,mAtoms),
     &          iBM(nB_Tot), idBM(2,ndB_Tot), mqB(nB)
      Logical Smmtrc(3,nAtoms), Process, Proc_dB,Help, R_Stab_A
      Character*14 Label, qLbl(nB)
      Character*3 ChOp(0:7)
      Character*(LENIN) Name(nAtoms)
      Character*(LENIN4) Lbls(2)
#include "bondtypes.fh"
#define _FMIN_
#define _VDW_
#include "ddvdt.fh"
#define _SCHLEGEL_
#include "ddvdt_bond.fh"
      Data ChOp/'E  ','X ','Y ','XY ','Z  ','XZ ','YZ ','XYZ'/
      Data iChOp/1,1,1,2,1,2,2,3/
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
      If (nBonds.lt.1) Return
*
*
      nqB=0
      Write (6,*)
      Write (6,*) ' ---> Enter Bonds.'
      Write (6,*)
      Write (6,*) 'Process=',Process
      Call RecPrt('CX',' ',CX,3*nAtoms,nIter)
      Write (6,'(20(1X,A))') (Name(i),i=1,nAtoms)
      Write (6,*)
      Write (6,*) ' iTabAI'
      Write (6,*)
      Do iAtom = 1, mAtoms
         Write (6,*) iTabAI(1,iAtom),iTabAI(2,iAtom)
      End Do
#endif
*
*---- Loop over bonds
*
      nCent=2
      Do iBond = 1, nBonds
         iBondType=iTabBonds(3,iBond)
*
*        We will only incorpotate covalent and fragment bonds
*
         If (iBondType.eq.vdW_Bond  ) Go To 1   ! vdW bonds
         If (iBondType.gt.Magic_Bond) Go To 1   ! magic bonds
*
         Do iCase = 1, 2
*
            If (iCase.eq.1) Then
               iAtom_ = iTabBonds(1,iBond)
               jAtom_ = iTabBonds(2,iBond)
            Else
               iAtom_ = iTabBonds(2,iBond)
               jAtom_ = iTabBonds(1,iBond)
            End If
            iAtom = iTabAI(1,iAtom_)
            jAtom = iTabAI(1,jAtom_)

            iDCR(1)=iTabAI(2,iAtom_)
            iDCR(2)=iTabAI(2,jAtom_)
            If (jAtom.gt.iAtom) Go To 2
            If (iDCR(1).ne.iOper(0)) Go To 2
            If (R_Stab_A(iDCR(2),jStab(0,iAtom),nStab(iAtom)) .and.
     &          iDCR(2).ne.iOper(0)) Go To 2
            iRow = iANr(iAtom)
            jRow = iANr(jAtom)
#ifdef _DEBUGPRINT_
            Write (6,*) 'iAtom,jAtom=',iAtom,jAtom
#endif
            Help = iRow.gt.3 .or. jRow.gt.3
            Ind(1)=iAtom
            Ind(2)=jAtom
            call dcopy_(3,Cx(1,iAtom,iIter),1,A,1)
            Write (Label,'(A,I2,A,I2,A)') 'B(',iAtom,',',jAtom,')'
*
#ifdef _DEBUGPRINT_
            Call RecPrt('A',' ',Cx(1,iAtom,iIter),1,3)
            Call RecPrt('B',' ',Cx(1,jAtom,iIter),1,3)
#endif

*
*------------- Form double coset representatives
*
            Call DCR(Lambda,
     &               jStab(0,iAtom),nStab(iAtom),
     &               jStab(0,jAtom),nStab(jAtom),iDCRR,nDCRR)
             kDCRR = iDCR(2)
*
#ifdef _DEBUGPRINT_
            Write (6,'(10A)') 'U={',
     &            (ChOp(jStab(i,iAtom)),i=0,nStab(iAtom)-1),'}  '
            Write (6,'(10A)') 'V={',
     &            (ChOp(jStab(i,jAtom)),i=0,nStab(jAtom)-1),'}  '
            Write (6,'(10A)') 'R={',
     &            (ChOp(iDCRR(i)),i=0,nDCRR-1),'}  '
            Write (6,'(2A)') 'R=',ChOp(iDCR(2))
#endif
*
            Call OA(iDCR(2),Cx(1:3,jAtom,iIter),A(1:3,2))
*
*---------- Compute the stabilizer of A & R(B), this is done in two ways.
*
*           A=/=B, the stabilizer is formed as the intersection of
*                  the stabilizers of A and B.
*
*           A=B, the stabilizer is formed as union of U and R(U)
*
            If (iAtom.ne.jAtom) Then
               Call Inter(jStab(0,iAtom),nStab(iAtom),
     &                    jStab(0,jAtom),nStab(jAtom),iStabM,nStabM)
            Else
               Call Union(jStab(0,iAtom),nStab(iAtom),
     &                    jStab(0,jAtom),nStab(jAtom),
     &                    kDCRR,iStabM,nStabM)
            End If
#ifdef _DEBUGPRINT_
            Write (6,'(10A)') 'M={',
     &            (ChOp(iStabM(i)),i=0,nStabM-1),'}  '
#endif
*
*---------- Now evaluate the degeneracy of the bond.
*
            iDeg=nIrrep/nStabM
            Deg=Sqrt(DBLE(iDeg))
#ifdef _DEBUGPRINT_
            Write (6,*)' nIrrep,nStabM=',nIrrep,nStabM
#endif
*
            nq = nq + 1
            If (.Not.Process) mB_Tot = mB_Tot + mB
            If (.Not.Proc_dB) mdB_Tot = mdB_Tot + mB**2
*
            nqB = nqB + 1
            iF1=1
            Call NxtWrd(Name(iAtom),iF1,iE1)
            Lbls(1)=Name(iAtom)(iF1:iE1)
            iF2=1
            Call NxtWrd(Name(jAtom),iF2,iE2)
            Lbls(2)=Name(jAtom)(iF2:iE2)
            If (kDCRR.ne.0) Then
               Lbls(2)(iE2+1:iE2+1)='('
               Lbls(2)(iE2+2:iE2+1+iChOp(kDCRR))=
     &               ChOp(kDCRR)(1:iChOp(kDCRR))
               Lbls(2)(iE2+2+iChOp(kDCRR):iE2+2+iChOp(kDCRR))=')'
               Call NxtWrd(Lbls(2),iF2,iE2)
            End If
            Write (LuIC,'(A,I3.3,4A)')
     &             'b',nqB,' = Bond ',
     &             Lbls(1)(iF1:iE1),' ',
     &             Lbls(2)(iF2:iE2)
#ifdef _DEBUGPRINT_
            Write (6,'(A,I3.3,4A)')
     &             'b',nqB,' = Bond ',
     &             Lbls(1)(iF1:iE1),' ',
     &             Lbls(2)(iF2:iE2)
#endif
            Label=' '
            Write (Label,'(A,I3.3)') 'b',nqB
            If (.Not.Proc_dB) Call FZero(Hess,36)
            Call Strtch(A,nCent,Val,Grad,.False.,'        ',Hess,
     &                  Proc_dB)
*
            If (Process) Then
*
               Indq(1,nq) = 1
               Indq(2,nq) = (jAtom-1)*nAtoms + iAtom
               Indq(3,nq) = kDCRR+1
*
               Rij2=(A(1,1)-A(1,2))**2
     &             +(A(2,1)-A(2,2))**2
     &             +(A(3,1)-A(3,2))**2
               Rab=Sqrt(Rij2)
               If (Help) Then
                  RabCov=CovRad(iANr(iAtom))+CovRad(iANr(jAtom))
                  If ((iRow.eq.1.and.jRow.eq.1).or.Help) Then
*                    Bond a la Fischer & Almlof
                     f_Const=A_StrH(1)*EXP(-A_StrH(2)*(Rab-RabCov))
                  Else
                     ij=Max(iRow,jRow)*(Max(iRow,jRow)+1)/2
     &                 +Min(iRow,jRow)
                     f_Const=A_Str/(Rab-B_Str(ij))**3
                  End If
               Else
                  If (iAnd(iOptC,2048).eq.2048.and.iBondType.eq.1) Then
                     r0 = r_ref_vdW(iRow,jRow)
                     f_Const=rkr_vdW*EXP(-Alpha_vdW*(Rab-r0)**2)
                  Else
                     r0 = rAV(iRow,jRow)
                     Alpha = aAV(iRow,jRow)
                     f_Const=rkr*EXP(Alpha*(r0**2-rij2))
                  End If
               End If
*
               f_Const=Max(f_Const,f_Const_Min)
               If (iBondType.eq.Fragments_Bond) f_Const=f_Const*1.D3
               fconst(nq)=Sqrt(f_Const)
               rMult(nq)=Deg
*
               Value(nq,iIter)=Val
               qLbl(nq)=Label
*
*-----------   Project the gradient vector
*
               Call ProjSym(nAtoms,nCent,Ind,nStab,jStab,A,
     &                      iDCR,Grad,Smmtrc,nDim,
     &                      Hess,mB_Tot,mdB_Tot,
     &                      BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,
     &                      Proc_dB,mqB,nB,nq,rMult(nq))
*
            End If
*
 2          Continue
         End Do     ! iCase
*
 1       Continue
      End Do        ! iBond
*
      Return
      End
