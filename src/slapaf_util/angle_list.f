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
      Subroutine Angle_List(
     &                 nq,
     &                 nAtoms,iIter,nIter,Cx,jStab,
     &                 nStab,nDim,Smmtrc,Process,Value,
     &                 nB,iANr,qLbl,iRef,
     &                 fconst,rMult,LuIC,Name,Indq,
     &                 Grad_all,iGlow,iGhi,iPrv,Proc_dB,
     &                 iTabBonds,nBonds,iTabAI,mAtoms,iTabAtoms,nMax,
     &                 mB_Tot,mdB_Tot,
     &                 BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,nqB,Thr_small)
      use Symmetry_Info, only: nIrrep, iOper
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "print.fh"
#include "Molcas.fh"
      Parameter (mB = 3*3)
      Real*8 Cx(3,nAtoms,nIter), A(3,3), Hess(mB**2),
     &       fconst(nB),Value(nB,nIter), rMult(nB),
     &       Ref(3,3), Prv(3,3),
     &       Grad_Ref(9), Axis(3), Perp_Axis(3,2), Grad(mB),
     &       Grad_all(9,iGlow:iGhi,nIter),
     &       BM(nB_Tot), dBM(ndB_Tot)
      Integer   nStab(nAtoms), iANr(nAtoms),
     &          iDCRR(0:7), jStab(0:7,nAtoms),
     &          iStabM(0:7), Ind(3), iDCR(3), iDCRT(0:7),
     &          iStabN(0:7), iChOp(0:7), Indq(3,nB), nqB(nB),
     &          iTabBonds(3,nBonds), iTabAI(2,mAtoms),
     &          iTabAtoms(2,0:nMax,mAtoms), iBM(nB_Tot),idBM(2,ndB_Tot)
      Logical Smmtrc(3,nAtoms), Process, PSPrint,
     &        MinBas, Help, Proc_dB, R_Stab_A
      Character*14 Label, qLbl(nB)
      Character*3 ChOp(0:7)
      Character*(LENIN) Name(nAtoms)
      Character*(LENIN4) Lbls(3)
#include "bondtypes.fh"
#define _FMIN_
#include "ddvdt.fh"
#include "ddvdt_bend.fh"
      Data ChOp/'E  ','X ','Y ','XY ','Z  ','XZ ','YZ ','XYZ'/
      Data iChOp/1,1,1,2,1,2,2,3/
#include "constants.fh"
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
      If (nBonds.lt.2) Return
      iRout=150
      iPrint=nPrint(iRout)
#ifdef _DEBUGPRINT_
      iPrint=99
#endif
*
      nqA=0
      PSPrint=.False.
#ifdef _DEBUGPRINT_
      If (iPrint.ge.99) PSPrint=.True.
      If (PSPrint) Write (6,*) ' Enter Bends.'
#endif
      Call FZero(Hess,81)
*
*---- Loop over bends
*
      bohr=CONST_BOHR_RADIUS_IN_SI_ * 1.0D+10
      MinBas=.False.
      If (MinBas) Then
         Fact=1.3d0
      Else
         Fact=One
      End If
      nCent=3
*     Write (6,*)
*
      Do mAtom_ = 1, mAtoms
         mAtom = iTabAI(1,mAtom_)
         mr = iTabRow(iANr(mAtom))
         Ind(2) = mAtom
         iDCR(2) = iTabAI(2,mAtom_)

         nNeighbor_m = iTabAtoms(1,0,mAtom_)
         nCoBond_m=nCoBond(mAtom_,mAtoms,nMax,iTabBonds,nBonds,
     &                     nBonds,iTabAtoms)
         If (nNeighbor_m.lt.2) Go To 100
*
         Do iNeighbor = 1, nNeighbor_m
            iAtom_ = iTabAtoms(1,iNeighbor,mAtom_)
            iAtom = iTabAI(1,iAtom_)
            nNeighbor_i = iTabAtoms(1,0,iAtom_)
            nCoBond_i=nCoBond(iAtom_,mAtoms,nMax,iTabBonds,nBonds,
     &                        nBonds,iTabAtoms)
            ir = iTabRow(iANr(iAtom))
            Ind(1) = iAtom
            iDCR(1) = iTabAI(2,iAtom_)
*
            If (iDCR(1).ne.iOper(0)) Go To 200
            If (R_Stab_A(iDCR(2),jStab(0,iAtom),nStab(iAtom)) .and.
     &          iDCR(2).ne.iOper(0)) Go To 200
*
            iBond = iTabAtoms(2,iNeighbor,mAtom_)
            iBondType=iTabBonds(3,iBond)
            If (iBondType.eq.vdW_Bond.or.
     &          iBondType.gt.Magic_Bond) Go To 200
#ifdef _DEBUGPRINT_
            Write (6,*)
            Write (6,*) 'iAtom,mAtom=',iAtom,mAtom
            Write (6,*) 'iBond,iBondType=',iBond,iBondType
            Write (6,*) 'E,R=',ChOp(iDCR(1)),ChOp(iDCR(2))
#endif
*
            call dcopy_(3,Cx(1,iAtom,iIter),1,A,1)
            call dcopy_(3,Cx(1,iAtom,iRef),1,Ref,1)
            call dcopy_(3,Cx(1,iAtom,iPrv),1,Prv,1)
*
            Do jNeighbor = 1, nNeighbor_m
               jAtom_ = iTabAtoms(1,jNeighbor,mAtom_)
               jAtom = iTabAI(1,jAtom_)
               nNeighbor_j = iTabAtoms(1,0,jAtom_)
               nCoBond_j=nCoBond(jAtom_,mAtoms,nMax,iTabBonds,nBonds,
     &                           nBonds,iTabAtoms)
               If (nCoBond_i.ge.8 .and.
     &             nCoBond_j.ge.8 .and.
     &             nCoBond_m.ge.8       ) Go To 300
*
               jr = iTabRow(iANr(jAtom))
               Ind(3) = jAtom
               iDCR(3) = iTabAI(2,jAtom_)
               If (R_Stab_A(iDCR(3),jStab(0,iAtom),nStab(iAtom)) .and.
     &             iDCR(3).ne.iOper(0)) Go To 300
               If (iDCR(3).eq.iOper(0) .and. iAtom.ge.jAtom) Go To 300
               kDCR=iEor(iDCR(2),iDCR(3))
               If (R_Stab_A(kDCR,jStab(0,mAtom),nStab(mAtom)) .and.
     &             iDCR(2).ne.iOper(0)) Go To 300
*
               Help=mr.gt.3 .or. ir.gt.3 .or. jr.gt.3
*
               jBond = iTabAtoms(2,jNeighbor,mAtom_)
               jBondType=iTabBonds(3,jBond)
               If (jBondType.eq.vdW_Bond.or.
     &             jBondType.gt.Magic_Bond) Go To 300
#ifdef _DEBUGPRINT_
               Write (6,*)
               Write (6,*) 'jAtom,mAtom=',jAtom,mAtom
               Write (6,*) 'jBond,jBondType=',jBond,jBondType
               Write (6,*) 'T=',ChOp(iDCR(3))
#endif
*
               Write (Label,'(A,I2,A,I2,A,I2,A)')
     &                'A(',iAtom,',',mAtom,',',jAtom,')'
*
#ifdef _DEBUGPRINT_
               If (PSPrint) Then
                  Call RecPrt('A',' ',Cx(1,iAtom,iIter),1,3)
                  Call RecPrt('B',' ',Cx(1,mAtom,iIter),1,3)
                  Call RecPrt('C',' ',Cx(1,jAtom,iIter),1,3)
               End If
#endif
*
*------------- Form double coset representatives for (iAtom,jAtom)
*
               Call DCR(Lambda,
     &                  jStab(0,iAtom),nStab(iAtom),
     &                  jStab(0,jAtom),nStab(jAtom),iDCRT,nDCRT)
#ifdef _DEBUGPRINT_
               Write (6,'(10A)') 'T={',
     &               (ChOp(iDCRT(i)),i=0,nDCRT-1),'}  '
#endif
               kDCRT=iDCR(3)
*
#ifdef _DEBUGPRINT_
               If (PSPrint) Then
                  Write (6,'(10A)') 'U={',
     &                  (ChOp(jStab(i,iAtom)),i=0,nStab(iAtom)-1),'}  '
                  Write (6,'(10A)') 'V={',
     &                  (ChOp(jStab(i,mAtom)),i=0,nStab(mAtom)-1),'}  '
                  Write (6,'(10A)') 'X={',
     &                  (ChOp(jStab(i,jAtom)),i=0,nStab(jAtom)-1),'}  '
                  Write (6,'(2A)') 'T=',ChOp(kDCRT)
               End If
#endif
*
               Call OA(kDCRT,Cx(1:3,jAtom,iIter),  A(1:3,3))
               Call OA(kDCRT,Cx(1:3,jAtom,iRef ),Ref(1:3,3))
               Call OA(kDCRT,Cx(1:3,jAtom,iPrv ),Prv(1:3,3))
*
*------------- Form the stabilizer for (iAtom,jAtom)
*
               If (iAtom.eq.jAtom) Then
                     Call Union(jStab(0,iAtom),nStab(iAtom),
     &                          jStab(0,jAtom),nStab(jAtom),
     &                          iDCR(3),iStabN,nStabN)
               Else
                     Call Inter(jStab(0,iAtom),nStab(iAtom),
     &                          jStab(0,jAtom),nStab(jAtom),
     &                          iStabN,nStabN)
               End If
*
#ifdef _DEBUGPRINT_
               If (PSPrint) Then
                  Write (6,'(10A)') 'N={',
     &                  (ChOp(iStabN(i)),i=0,nStabN-1),'}  '
               End If
#endif
*
*------------- Form double coset representatives for
*              ((iAtom,mAtom),jAtom)
*
               Call DCR(Lambda,
     &                  jStab(0,mAtom),nStab(mAtom),
     &                  iStabN,nStabN,iDCRR,nDCRR)
               kDCRR = iDCR(2)
*
#ifdef _DEBUGPRINT_
               If (PSPrint) Then
                  Write (6,'(10A)') 'R={',
     &                  (ChOp(iDCRR(i)),i=0,nDCRR-1),'}  '
                  Write (6,'(2A)') 'R=',ChOp(kDCRR)
               End If
#endif
*
               Call OA(kDCRR,Cx(1:3,mAtom,iIter),  A(1:3,2))
               Call OA(kDCRR,Cx(1:3,mAtom,iRef ),Ref(1:3,2))
               Call OA(kDCRR,Cx(1:3,mAtom,iPrv ),Prv(1:3,2))
*
*------------- Form the stabilizer for ((iAtom,mAtom),jAtom)
*
               Call Inter(jStab(0,mAtom),nStab(mAtom),
     &                    iStabN,nStabN,
     &                    iStabM,nStabM)
*
#ifdef _DEBUGPRINT_
               If (PSPrint) Then
                  Write (6,'(10A)') 'M={',
     &                  (ChOp(iStabM(i)),i=0,nStabM-1),'}  '
               End If
#endif
*
*------------- Compute the degeneracy of the angle
*
               ideg=nIrrep/nStabM
               Deg=Sqrt(DBLE(iDeg))
#ifdef _DEBUGPRINT_
               If (PSPrint) Write (6,*)' nIrrep,nStabM=',nIrrep,nStabM
#endif
*
*------------- Test if coordinate should be included
*              All angles which are lower than the explicit threshold
*              are rejected.
*
               If (Help) Then
                  rim2=(Ref(1,1)-Ref(1,2))**2
     &                +(Ref(2,1)-Ref(2,2))**2
     &                +(Ref(3,1)-Ref(3,2))**2
                  Rab=Sqrt(rim2)
                  RabCov=CovRad(iANr(iAtom))+CovRad(iANr(mAtom))
                  rmj2=(Ref(1,2)-Ref(1,3))**2
     &                +(Ref(2,2)-Ref(2,3))**2
     &                +(Ref(3,2)-Ref(3,3))**2
                  Rbc=Sqrt(rmj2)
                  RbcCov=CovRad(iANr(jAtom))+CovRad(iANr(mAtom))
                  If (ir.eq.1.or.jr.eq.1) Then
                     f_Const=A_Bend(1)
                  Else
                     f_Const=A_Bend(2)
                  End If
                  f_Const=f_Const*Fact
                  f_Const_Ref=f_Const
               Else
                  r0=rAV(ir,mr)
                  Alpha=aAv(ir,mr)
                  rim2_Ref=(Ref(1,1)-Ref(1,2))**2
     &                    +(Ref(2,1)-Ref(2,2))**2
     &                    +(Ref(3,1)-Ref(3,2))**2
                  f_Const_Ref=rkf*Exp(Alpha*(r0**2-rim2_Ref))
                  rim2=(A(1,1)-A(1,2))**2
     &                +(A(2,1)-A(2,2))**2
     &                +(A(3,1)-A(3,2))**2
                  f_Const=rkf*Exp(Alpha*(r0**2-rim2))
*
                  r0=rAV(mr,jr)
                  Alpha=aAv(mr,jr)
                  rmj2_Ref=(Ref(1,2)-Ref(1,3))**2
     &                    +(Ref(2,2)-Ref(2,3))**2
     &                    +(Ref(3,2)-Ref(3,3))**2
                  f_Const_Ref=f_Const_Ref*Exp(Alpha*(r0**2-rmj2_Ref))
                  rmj2=(A(1,2)-A(1,3))**2
     &                +(A(2,2)-A(2,3))**2
     &                +(A(3,2)-A(3,3))**2
                  f_Const=f_Const*Exp(Alpha*(r0**2-rmj2))
*
               End If
               If (f_Const_Ref.lt.f_Const_Min .and.
     &             iBondType.ne.Fragments_Bond .and.
     &             jBondType.ne.Fragments_Bond ) Go To 300
#ifdef _DEBUGPRINT_
               Write (6,*) ' A Force Constant:',f_Const
               Write (6,*) iAtom,mAtom,jAtom, f_Const
#endif
*
               Call Bend(Ref,nCent,Val_Ref,Grad_Ref,.False.,
     &                   .False.,'        ',Hess,Proc_dB)
               Call Bend(A,nCent,Val,Grad,.False.,
     &                   .False.,'        ',Hess,Proc_dB)
*
*------------- Skip cases with a too small angle.
*
               If (Abs(Val_Ref).lt.Thr_Small .and.
     &             iBondType.ne.Fragments_Bond .and.
     &             jBondType.ne.Fragments_Bond) Go To 300
*
*                                                                      *
************************************************************************
*                                                                      *
*              We'd like to avoid the problem with angles that have the
*              value Pi. We introduce "linear" angles under two
*              conditions.
*              (1) the reference angle is within delta of Pi.
*              (2) the actual angle is within 1.0D-11 of Pi.(?)
*              Delta is set to 1.0D-11 if there are 3 atoms.(?)
*
               delta=(45.0D0/180.0D0)*Pi
*              If (mAtoms.eq.3) delta=1.0D-11 ! I do not understand
*              although it is probably me who introduced it!
*
               If ( (Abs(Val_Ref-Pi).lt.Delta .or.
     &               Abs(Val-Pi).lt.1.0D-11) .and.
     &            .NOT.(
     &                  (iBondType.eq.Fragments_Bond .or.
     &                   jBondType.eq.Fragments_Bond)
     &                   .and. mAtoms.le.4
     &                  )
     &            ) Then
*
*---------------- Reference is linear(a) or
*---------------- reference is NOT linear but the new structure is(b).
*
                  If (Abs(Val-Pi).lt.1.0D-11 .and.
     &                .Not.(Abs(Val_Ref-Pi).lt.Delta)) Then
*                    Case b
                     nk=1
                  Else
*                    Case a
                     nk=2
                  End If
                  Call CoSys(Prv,Axis,Perp_Axis)
*
                  Label=' '
                  Label(1:1)='L'
C                 Do k = 1, 2
                  Do k = 1, nk
                     nq = nq + 1
                     If (.Not.Process) mB_Tot = mB_Tot + mB
                     If (.Not.Proc_dB) mdB_Tot = mdB_Tot + mB**2
*
                     nqA = nqA + 1
                     iF1=1
                     Call NxtWrd(Name(iAtom),iF1,iE1)
                     Lbls(1)=Name(iAtom)(iF1:iE1)
                     iF2=1
                     Call NxtWrd(Name(mAtom),iF2,iE2)
                     Lbls(2)=Name(mAtom)(iF2:iE2)
                     If (kDCRR.ne.0) Then
                        Lbls(2)(iE2+1:iE2+1)='('
                        Lbls(2)(iE2+2:iE2+1+iChOp(kDCRR))=
     &                        ChOp(kDCRR)(1:iChOp(kDCRR))
                        Lbls(2)(iE2+2+iChOp(kDCRR):
     &                          iE2+2+iChOp(kDCRR))=')'
                        Call NxtWrd(Lbls(2),iF2,iE2)
                     End If
                     iF3=1
                     Call NxtWrd(Name(jAtom),iF3,iE3)
                     Lbls(3)=Name(jAtom)(iF3:iE3)
                     If (kDCRT.ne.0) Then
                        Lbls(3)(iE3+1:iE3+1)='('
                        Lbls(3)(iE3+2:iE3+1+iChOp(kDCRT))=
     &                        ChOp(kDCRT)(1:iChOp(kDCRT))
                        Lbls(3)(iE3+2+iChOp(kDCRT):
     &                          iE3+2+iChOp(kDCRT))=')'
                        Call NxtWrd(Lbls(3),iF3,iE3)
                     End If
                     Write (LuIC,'(A,I3.3,A,I1.1,6A)')
     &                      'a',nqA,' = LAngle(',k,') ',
     &                          Lbls(1)(iF1:iE1),
     &                      ' ',Lbls(2)(iF2:iE2),
     &                      ' ',Lbls(3)(iF3:iE3)
#ifdef _DEBUGPRINT_
                     If (iPrint.ge.49)
     &               Write (6,'(A,I3.3,A,I1.1,6A)')
     &                      'a',nqA,' = LAngle(',k,') ',
     &                          Lbls(1)(iF1:iE1),
     &                      ' ',Lbls(2)(iF2:iE2),
     &                      ' ',Lbls(3)(iF3:iE3)
                     Write (6,*) 'iDeg=',iDeg
#endif
                     Label=' '
                     Write (Label,'(A,I3.3)') 'a',nqA
*
                     If (Process) Then
*
                        Call LBend(A,nCent,Val,
     &                             Grad_all(1,nq,iIter),
     &                             .False.,.False.,
     &                             '        ',Hess,Proc_dB,Axis,
     &                             Perp_Axis(1,k),(k.eq.2))
*
*---------------------- Flip Angle value and gradient if needed!
*
                        BB=DDot_(9,Grad_all(1,nq, iPrv),1,
     &                             Grad_all(1,nq,iIter),1)
                        If (BB.lt.Zero) Then
*                          Write (6,*) ' Angle flips, corrected!'
                           Val=Two*Pi-Val
                           Call DScal_(9,-One,
     &                                Grad_all(1,nq,iIter),1)
                           Call DScal_(81,-One,Hess,1)
                        End If
*
                        Indq(1,nq)=2+k
                        mi = (iAtom-1)*nAtoms + mAtom
                        Indq(2,nq)= (jAtom-1)*nAtoms**2 + mi
                        Indq(3,nq)= kDCRT*8 + kDCRR+1
*
                        f_Const=Max(f_Const,f_Const_Min)
                        If (.Not.Help .and.
     &                     (iBondType.eq.Fragments_Bond .or.
     &                      jBondType.eq.Fragments_Bond ))
     &                      f_Const=f_Const*1.D3
                        fconst(nq)=Sqrt(f_Const)
                        rMult(nq)=Deg
*
                        Value(nq,iIter)=Val
                        qLbl(nq)=Label
*
*---------------------- Project the gradient vector
*
                        Call ProjSym(nAtoms,nCent,Ind,nStab,
     &                               jStab,A,iDCR,
     &                               Grad_all(1,nq,iIter),
     &                               Smmtrc,nDim,PSPrint,Hess,
     &                               mB_Tot,mdB_Tot,
     &                               BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,
     &                               Proc_dB,nqB,nB,nq,rMult(nq))
*
                     End If
*
                  End Do
*
               Else      ! Non-linear case
*
                  nq = nq + 1
                  If (.Not.Process) mB_Tot = mB_Tot + mB
                  If (.Not.Proc_dB) mdB_Tot = mdB_Tot + mB**2
*
                  nqA = nqA + 1
                  iF1=1
                  Call NxtWrd(Name(iAtom),iF1,iE1)
                  Lbls(1)=Name(iAtom)(iF1:iE1)
                  iF2=1
                  Call NxtWrd(Name(mAtom),iF2,iE2)
                  Lbls(2)=Name(mAtom)(iF2:iE2)
                  If (kDCRR.ne.0) Then
                     Lbls(2)(iE2+1:iE2+1)='('
                     Lbls(2)(iE2+2:iE2+1+iChOp(kDCRR))=
     &                     ChOp(kDCRR)(1:iChOp(kDCRR))
                     Lbls(2)(iE2+2+iChOp(kDCRR):
     &                       iE2+2+iChOp(kDCRR))=')'
                     Call NxtWrd(Lbls(2),iF2,iE2)
                  End If
                  iF3=1
                  Call NxtWrd(Name(jAtom),iF3,iE3)
                  Lbls(3)=Name(jAtom)(iF3:iE3)
                  If (kDCRT.ne.0) Then
                     Lbls(3)(iE3+1:iE3+1)='('
                     Lbls(3)(iE3+2:iE3+1+iChOp(kDCRT))=
     &                     ChOp(kDCRT)(1:iChOp(kDCRT))
                     Lbls(3)(iE3+2+iChOp(kDCRT):
     &                       iE3+2+iChOp(kDCRT))=')'
                     Call NxtWrd(Lbls(3),iF3,iE3)
                  End If
                  Write (LuIC,'(A,I3.3,6A)')
     &                   'a',nqA,' = Angle ',
     &                       Lbls(1)(iF1:iE1),
     &                   ' ',Lbls(2)(iF2:iE2),
     &                   ' ',Lbls(3)(iF3:iE3)
#ifdef _DEBUGPRINT_
                  If (iPrint.ge.49)
     &            Write (6,'(A,I3.3,6A)')
     &                   'a',nqA,' = Angle ',
     &                       Lbls(1)(iF1:iE1),
     &                   ' ',Lbls(2)(iF2:iE2),
     &                   ' ',Lbls(3)(iF3:iE3)
#endif
                  Label=' '
                  Write (Label,'(A,I3.3)') 'a',nqA
*
                  If (Process) Then
*
                     Call Bend(A,nCent,Val,Grad_all(1,nq,iIter),
     &                         .False.,.False.,'        ',Hess,
     &                         Proc_dB)
*
*---------------------- Flip Angle value and gradient if needed!
*
                     BB=DDot_(9,Grad_all(1,nq,iPrv),1,
     &                       Grad_all(1,nq,iIter),1)
                     If (BB.lt.Zero) Then
*                       Write (*,*) ' Angle flips, corrected!'
*                       Write (*,*) ' iRef,iIter=', iRef,iIter
                        Val=Two*Pi-Val
                        Call DScal_(9,-One,Grad_all(1,nq,iIter),1)
                        Call DScal_(81,-One,Hess,1)
                     End If
*
                     Indq(1,nq)=2
                     mi = (iAtom-1)*nAtoms + mAtom
                     Indq(2,nq)= (jAtom-1)*nAtoms**2 + mi
                     Indq(3,nq)= kDCRT*8 + kDCRR+1
*
                     If (.Not.Help .and.
     &                  (iBondType.eq.Fragments_Bond .or.
     &                   jBondType.eq.Fragments_Bond ))
     &                   f_Const=f_Const*100.D0
                     fconst(nq)=Sqrt(f_Const)
                     rMult(nq)=Deg
*
                     Value(nq,iIter)=Val
                     qLbl(nq)=Label
*
*------------------- Project the gradient vector
*
                     Call ProjSym(nAtoms,nCent,Ind,nStab,
     &                            jStab,A,iDCR,
     &                            Grad_all(1,nq,iIter),
     &                            Smmtrc,nDim,PSPrint,Hess,
     &                            mB_Tot,mdB_Tot,
     &                            BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,
     &                            Proc_dB,nqB,nB,nq,rMult(nq))
*
                  End If
*
               End If
*                                                                      *
************************************************************************
*                                                                      *
*
  300          Continue
            End Do             ! End loop over jNeighbor
  200       Continue
         End Do                ! End loop over iCase
  100    Continue
      End Do                   ! End loop over iBond
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
