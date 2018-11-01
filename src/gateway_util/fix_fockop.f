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
* Copyright (C) 2017, Roland Lindh                                     *
************************************************************************
      Subroutine Fix_FockOp(Info,nInfo,LuRd,DInf,nDInf)
************************************************************************
*                                                                      *
*    Objective: To compute the fock operator for basis sets which do   *
*               not carry this information explicitly in the basis set *
*               file or are inline basis sets.                         *
*                                                                      *
*    The Fock operator is defined as F=\sum_{k,m} |k>e_{k,m}<m|.       *
*                                                                      *
* Called from: Input                                                   *
*                                                                      *
* Calling    : GetBS                                                   *
*                                                                      *
*     Author:  Roland Lindh (thanks to Ben Swerts)                     *
*                                                                      *
*     Modified for muons by R. Lindh January 2017                      *
************************************************************************
      use Her_RW
      use Real_Spherical
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "stdalloc.fh"
#include "real.fh"
#include "print.fh"
#include "status.fh"
#include "periodic_table.fh"
      External MltPrm, KnEPrm, NAPrm
      Real*8 DInf(nDInf)
      Real*8, Dimension(:), Allocatable :: FockOp_t
      Character*13 DefNm
      Character*80 Ref(2), Bsl_, BSLbl
      Character *256 Basis_lib, Fname
      Character*180 STDINP(mxAtom*2) ! CGGn
      Integer BasisTypes(4), nDel(MxAng)
      Integer List_AE(0:iTabMx), List(0:iTabMx), List_Add(0:iTabMx)
      Logical Try_Again
      Real*8 A(4)
      Data DefNm/'basis_library'/
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUG_
#ifdef _DEBUG_
C     nPrint(113)=99
C     nPrint(114)=99
C     nPrint(116)=99
C     nPrint(122)=99
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Call qEnter('Fix_FockOp')
*                                                                      *
************************************************************************
*                                                                      *
*     Generate a dummy center. This is fine since we will only do
*     1-center overlap integrals here.
*
      call dcopy_(3,Zero,0,A,1)
*
      nOrdOp=2
      iComp = 1
*
      nPrp=Max(4,nMltpl)
      nDiff = 0
*
      Call ICopy(1+iTabMx,0,0,List   ,1)
      Call ICopy(1+iTabMx,0,0,List_AE,1)
      BasisTypes(1)=0
      BasisTypes(2)=0
      BasisTypes(3)=0
      BasisTypes(4)=0
      lSTDINP=0
*
*     Loop over all valence shell with a non-funtional FockOp
*
      mCnttp = nCnttp   ! to be restored at the end
*
      Do 1000 iCnttp = 1, mCnttp
*
         iFerm=1
         If (fMass(iCnttp).ne.1.0D0) iFerm=2
*
         If (FockOp(iCnttp).and.Charge(iCnttp).eq.0.0D0) Then
            Do iAng = 0, nVal_Shells(iCnttp)-1
               iShll_a    = ipVal(iCnttp) + iAng
               ipFockOp_a = ipFockOp(iShll_a)
               nCntrc_a = nBasis_Cntrct(iShll_a)
               Call FZero(Dinf(ipFockOp_a),nCntrc_a**2)
            End Do
         End If
*
         If(AuxCnttp(iCnttp) .or.
     &      FragCnttp(iCnttp) .or.
     &      nFragType(iCnttp).gt.0 .or.
     &      FockOp(iCnttp)) Then
           Goto 1000
         End If
*
*        Special treatment for muonic basis sets
*
         If (iFerm.eq.2) Then
*
            iShll = Mx_Shll-1
            jShll = iShll
*
*           The Fock operator will simply be the one-particle
*           Hamiltonian (kinetic + nuclear-attraction operator)
*
            xFactor=1.0D0/fMass(iCnttp)
            If (FNMC) Then
               iAtom=iAtmNr(iCnttp)
*              Get the atom mass in au (me=1)
               xMass=CntMass(iCnttp)
*              Substract the electron mass to get the nuclear mass.
               xMass=xMass-DBLE(iAtom)
               xfactor=xfactor+One/xMass
            End If
*
            Do iAng = 0, nVal_Shells(iCnttp)-1
*
               iShll_a    = ipVal(iCnttp) + iAng
               ipCff_a    = ipCff_Cntrct(iShll_a)
               ipExp_a    = ipExp(iShll_a)
               nPrim_a  = nExp(iShll_a)
               If (nPrim_a.eq.0) Cycle
               nCntrc_a = nBasis_Cntrct(iShll_a)
               iCmp_a = (iAng+1)*(iAng+2)/2
               If (Prjct(iShll_a)) iCmp_a = 2*iAng+1
               naa = nElem(iAng)*nElem(iAng)
               nScr1 = Max(nPrim_a,nPrim_a)*Max(nCntrc_a,nCntrc_a)*naa
               nScr2 = Max(nCntrc_a,nCntrc_a)**2*naa
*
               ipFockOp_a = ipFockOp(iShll_a)
*                                                                      *
************************************************************************
*                                                                      *
*              Compute the kinetic integrals
*
               nOrdOp=2
               ip = ipExp(iShll+1)
               Call One_Int(KnEPrm,DInf,nDInf,A,ip,Info,nInfo,jShll,
     &                      iAng,iComp,nOrdOp,nScr1,nScr2,naa,ipKnE,
     &                      nSAA,
     &                      iShll_a,nPrim_a,ipExp_a,nCntrc_a,ipCff_a,
     &                      iCmp_a,
     &                      iShll_a,nPrim_a,ipExp_a,nCntrc_a,ipCff_a,
     &                      iCmp_a)
*define _DEBUG_
#ifdef _DEBUG_
               Call DScal_(nCntrc_a**2*iCmp_a**2,
     &                     xFactor,DInf(ipKnE),1)
               Call RecPrt('Kinetric Energy Integrals',' ',
     &                     DInf(ipKnE),nCntrc_a**2,iCmp_a**2)
               Call DScal_(nCntrc_a**2*iCmp_a**2,
     &                     1.0D0/xFactor,DInf(ipKnE),1)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*              Compute the nuclear-attraction integrals
*
               nOrdOp=0
               A(4) = DBLE(iCnttp) ! Dirty tweak
               Call One_Int(NAPrm,DInf,nDInf,A,ip,Info,nInfo,jShll,
     &                      iAng,iComp,nOrdOp,nScr1,nScr2,naa,ipNAE,
     &                      nSBB,
     &                      iShll_a,nPrim_a,ipExp_a,nCntrc_a,ipCff_a,
     &                      iCmp_a,
     &                      iShll_a,nPrim_a,ipExp_a,nCntrc_a,ipCff_a,
     &                      iCmp_a)
#ifdef _DEBUG_
               Call RecPrt('Nuclear-attraction Integrals',' ',
     &                     DInf(ipNAE),nCntrc_a**2,iCmp_a**2)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*              Add together the kinetic and nuclear-attraction
*
               Call DaXpY_(nCntrc_a**2*iCmp_a**2,
     &                     xFactor,
     &                     DInf(ipKnE),1,
     &                     DInf(ipNAE),1)
*
*              Change to proper order (nCntrc_a * iCmp_a)
*
               jp1Hm = ip
               ip = ip + nCntrc_a**2 * iCmp_a**2
               Call Reorder(DInf(ipNAE),DInf(jp1Hm),
     &                      nCntrc_a,nCntrc_a,iCmp_a,iCmp_a)
*                                                                      *
************************************************************************
*                                                                      *
*              Compute the overlap integrals
*
               nOrdOp=0
               Call One_Int(MltPrm,DInf,nDInf,A,ip,Info,nInfo,jShll,
     &                      iAng,iComp,nOrdOp,nScr1,nScr2,naa,ipOvr,
     &                      nSCC,
     &                      iShll_a,nPrim_a,ipExp_a,nCntrc_a,ipCff_a,
     &                      iCmp_a,
     &                      iShll_a,nPrim_a,ipExp_a,nCntrc_a,ipCff_a,
     &                      iCmp_a)
#ifdef _DEBUG_
               Call RecPrt('Overlap Integrals',' ',
     &                     DInf(ipOvr),nCntrc_a**2,iCmp_a**2)
#endif
*
*              Change to proper order (nCntrc_a * iCmp_a)
*
               jpOvr = ip
               ip = ip + nCntrc_a**2 * iCmp_a**2
               Call Reorder(DInf(ipOvr),DInf(jpOvr),
     &                      nCntrc_a,nCntrc_a,iCmp_a,iCmp_a)
*                                                                      *
************************************************************************
*                                                                      *
*              Now we need to convert it to the Fock operator!
*
*              Solve F C = e S C
*
*              S^(-1/2) F S^(-1/2) S^(1/2) C = e S^(1/2) C
*              Set F' = S^(-1/2) F S^(-1/2)
*                  C' = S^(1/2) C
*
*              Solve F' C' = e C' , generate C = S^-(1/2) C'
*
               nBF = nCntrc_a*iCmp_a
               ipS = ip
               ip = ip + nBF**2
               ipS12i= ip
               ip = ip + nBF**2
*
*              1) Compute the eigenvectors and eigenvalues of the
*                 overlap matrix
*
               ipEVal = ip
               ip = ip + nBF*(nBF+1)/2
               ipEVec = ip
               ip = ip + nBF**2
               Call FZero(DInf(ipEVec),nBF**2)
               Call DCopy_(nBF,1.0D0,0,DInf(ipEVec),nBF+1)
               Do iBF = 1, nBF
                  Do jBF = 1, iBF
                     ij    =  (jBF-1)*nBF + iBF
                     ijTri = (iBF-1)*iBF/2 + jBF
                     DInf(ipEVal-1 + ijTri) = DInf(jpOvr-1 + ij)
                  End Do
               End Do
               Call NIDiag_new(DInf(ipEVal),DInf(ipEVec),nBF,nBF,0)
*
*              2) Construct S^(1/2) and S^(-1/2)
*
               Call FZero(DInf(ipS12i),nBF**2)
               Do kEval = 1, nBF
                  e   = DInf(ipEVal-1 + kEval*(kEval+1)/2)
                  e12i= 1.0D0/Sqrt(e)
                  Do iBF = 1, nBF
                     C_ik = DInf(ipEVec-1 + (kEVal-1)*nBF + iBF)
                     Do jBF = 1, nBF
                        C_jk = DInf(ipEVec-1 + (kEVal-1)*nBF + jBF)
                        ij = (jBF-1)*nBF + iBF
                        DInf(ipS12i-1 + ij) = DInf(ipS12i -1 + ij)
     &                                      + C_ik * e12i * C_jk
                     End Do
                  End Do
               End Do
*
*              3) Form F' =  S^(-1/2) F S^(-1/2)
*
               ipFPrim = ip
               ip = ip + nBF**2
               ipTemp = ip
               ip = ip + nBF**2
               Call FZero(Dinf(ipFPrim),nBF**2)
               Call DGEMM_('N','N',
     &                     nBF,nBF,nBF,
     &                     1.0d0,DInf(ipS12i),nBF,
     &                     DInf(jp1Hm),nBF,
     &                     0.0d0,DInf(ipTemp),nBF)
               Call DGEMM_('N','N',
     &                     nBF,nBF,nBF,
     &                     1.0d0,DInf(ipTemp),nBF,
     &                     DInf(ipS12i),nBF,
     &                     0.0d0,DInf(ipFPrim),nBF)
*
*              4) Compute C' and the eigenvalues
*
               Call FZero(DInf(ipEVec),nBF**2)
               Call DCopy_(nBF,1.0D0,0,DInf(ipEVec),nBF+1)
               Do iBF = 1, nBF
                  Do jBF = 1, iBF
                     ij    =  (jBF-1)*nBF + iBF
                     ijTri = (iBF-1)*iBF/2 + jBF
                     DInf(ipEVal-1 + ijTri) = DInf(ipFPrim-1 + ij)
                  End Do
               End Do
               Call NIDiag_new(DInf(ipEVal),DInf(ipEVec),nBF,nBF,0)
*
*              5) Form C = S^(-1/2) C'
*
               ipC = ip
               ip = ip + nBF**2
               Call DGEMM_('N','N',
     &                     nBF,nBF,nBF,
     &                     1.0d0,DInf(ipS12i),nBF,
     &                     DInf(ipEVec),nBF,
     &                     0.0d0,DInf(ipC),nBF)
#ifdef _DEBUG_
      Call RecPrt('Cs for F',' ',DInf(ipC),nBF,nBF)
#endif
*
*              6) Form the matrix representation of the Fock operator
*
               Call FZero(DInf(jp1Hm),nBF**2)
               Do kEval = 1, nBF
                  e   = DInf(ipEVal-1 + kEval*(kEval+1)/2)
                  Do iBF = 1, nBF
                     C_ik = DInf(ipC-1 + (kEVal-1)*nBF + iBF)
                     Do jBF = 1, nBF
                        C_jk = DInf(ipC-1 + (kEVal-1)*nBF + jBF)
                        ij = (jBF-1)*nBF + iBF
                        DInf(jp1Hm-1 + ij) = DInf(jp1Hm-1 + ij)
     &                                      + C_ik * e * C_jk
                     End Do
                  End Do
               End Do
*
               Call Reorder(DInf(jp1Hm),DInf(ipNAE),
     &                      nCntrc_a,iCmp_a,nCntrc_a,iCmp_a)
*
*              Make result isotropic and distribute
*
               Do iB = 1, nCntrc_a
                  Do jB = 1, nCntrc_a
                     ijB=(jB-1)*nCntrc_a+iB
                     iTo   = ipFockOp_a-1 + (jB-1)*nCntrc_a+iB
                     DInf(iTo) = Zero
                     Tmp = Zero
                     Do iC = 1, iCmp_a
                        ijC=(iC-1)*iCmp_a+iC
                        iFrom = ipNAE-1 + (ijC-1)*nCntrc_a**2+ijB
                        Tmp = Tmp + DInf(iFrom)
                     End Do
                     DInf(iTo) = DInf(iTo) + Tmp/DBLE(iCmp_a)
                  End Do
               End Do
#ifdef _DEBUG_
               Call RecPrt('Actual Fock operator',' ',DInf(ipFockOp_a),
     &                     nCntrc_a,nCntrc_a)
#endif
            End Do
            FockOp(iCnttp)=.TRUE.
            Go To 1000
         End If
*
*
*        create a new basis set index (temporary)
*
         nCnttp = mCnttp + 1
         If (nCnttp.gt.Mxdbsc) Then
            Call WarningMessage(2,'Fix_FockOp: Increase Mxdbsc')
            Call Abend()
         End If
*
*        create the temporary basis set label for this element to
*        read the corresponding ANO-RCC basis set.
*
         BSLbl=' '
         BSLbl=PTab(iAtmNr(iCnttp))
*
         If (BSLbl(1:1).eq.' ') Then
            BSLbl=BSLbl(2:2)//'.ANO-RCC.....'
         Else
            BSLbl=BSLbl(1:2)//'.ANO-RCC.....'
         End If
*
         LenBSL = Len(BSLbl)
         Last = iCLast(BSLbl,LenBSL)
         Indx = Index(BSLbl,'/')
*
         Bsl_=' '
         If (Indx.eq.0) Then
            Call WhichMolcas(Basis_lib)
            If (Basis_lib(1:1).ne.' ') then
               ib = index(Basis_lib,' ')-1
               If(ib.lt.1) Call SysAbendMsg('Read_ANO_RCC',
     &                     'Too long PATH to MOLCAS',' ')
               Fname=Basis_lib(1:ib)//'/basis_library'
            Else
               Fname=DefNm
            Endif
            Indx = Last+1
            Bsl_ = BSLbl
         Else
            Fname = BSLbl(Indx+2:Last)
            If (Fname.eq.' ') Then
               Call WarningMessage(2,
     &                     ' No basis set library specified for'
     &                         //',BSLbl='//BSLbl
     &                         //',Fname='//Fname)
               Call Quit_OnUserError()
            End If
 1001       If (Fname(1:1).eq.' ') Then
               Fname(1:79)=Fname(2:80)
               Fname(80:80) = ' '
               Go To 1001
            End If
            Bsl_=BSLbl(1:Indx-1)
         End If
*
#ifdef _DEBUG_
         Write (6,*)
         Write (6,*)
         Write(6,'(1X,A,I5,A,A)')
     &         'Basis Set ',nCnttp,' Label: ', BSLbl(1:Indx-1)
         Write(6,'(1X,A,A)') 'Basis set is read from library:',Fname
#endif
*
*        Let's get the reference basis set (ANO-RCC).
*
         iShll = Mx_Shll-1
         jShll = iShll
         SODK(nCnttp)=.False.
         Call GetBS(Fname,Bsl_,Indx-1,lAng,ipExp,
     &              ipCff,ipCff_Cntrct,ipCff_Prim,ipFockOp,nExp,
     &              nBasis,nBasis_Cntrct,MxShll,iShll,MxAng,
     &              Charge(nCnttp),iAtmNr(nCnttp),BLine,Ref,
     &              PAM2(nCnttp),
     &              ipPAM2xp(nCnttp),ipPAM2cf(nCnttp),nPAM2(nCnttp),
     &              FockOp(nCnttp),
     &              ECP(nCnttp),NoPairL(nCnttp),SODK(nCnttp),
     &              ipM1xp(nCnttp),ipM1cf(nCnttp),nM1(nCnttp),
     &              ipM2xp(nCnttp),ipM2cf(nCnttp),nM2(nCnttp),ipBk,
     &              CrRep(nCnttp),nProj,nAIMP,ipAkl,ip_Occ,iOptn,
     &              UnNorm,nDel,
     &              nVal,   nPrj,   nSRO,   nSOC,  nPP,
     &              ipVal_, ipPrj_, ipSRO_, ipSOC_,ipPP_,
     &              LuRd,BasisTypes,AuxCnttp(nCnttp),
     &              idummy,idummy,idummy,idummy,
     &              idummy,idummy,idummy,idummy,idummy,
     &              STDINP,lSTDINP,.False.,.true.,' ',
     &              DInf,nDInf)
*
         If (.Not.FockOp(nCnttp)) Then
            Write (6,*) 'Fix_FockOp: reference basis doesn''t contain a'
     &                //' proper Fock operator'
            Go To 1000
         End If
         Transf(jShll+1)=.False.
         Prjct(jShll+1)=.False.
         Transf(jShll+2)=.False.
         Prjct(jShll+2)=.False.
         ipVal(nCnttp) = ipVal_
         ipPrj(nCnttp) = ipPrj_
         ipSRO(nCnttp) = ipSRO_
         ipSOC(nCnttp) = ipSOC_
         ipPP(nCnttp)  = ipPP_
         nVal_Shells(nCnttp) = nVal
         nPrj_Shells(nCnttp) = nPrj
         nSRO_Shells(nCnttp) = nSRO
         nSOC_Shells(nCnttp) = nSOC
         nPP_Shells(nCnttp)  = nPP
         nTot_Shells(nCnttp) = nVal+nPrj+nSRO+nSOC+nPP
*                                                                      *
************************************************************************
*                                                                      *
*        Start processing shells of iCnttp and mCnttp. Loop only over
*        the shells of iCnttp (mCnttp might be larger!)
*
*
         Try_Again=.True.
         Call ICopy(1+iTabMx,0,0,List_Add,1)
 777     Continue
         Test_Charge=Zero
         Do iAng = 0, nVal_Shells(iCnttp)-1
*
*           Pointers to the actuall shell
*
            iShll_a    = ipVal(iCnttp) + iAng
            ipCff_a    = ipCff_Cntrct(iShll_a)
            ipExp_a    = ipExp(iShll_a)
            ipFockOp_a = ipFockOp(iShll_a)
            nPrim_a  = nExp(iShll_a)
ccjd
*           if (nPrim_a==0) cycle
            If (nPrim_a.eq.0) Go To 999
ccjd
            nCntrc_a = nBasis_Cntrct(iShll_a)
            iCmp_a = (iAng+1)*(iAng+2)/2
            If (Prjct(iShll_a)) iCmp_a = 2*iAng+1
*
*           Pointers to the reference shell
*
            iShll_r = ipVal(nCnttp) + iAng
            ipCff_r    = ipCff_Cntrct(iShll_r)
            ipExp_r    = ipExp(iShll_r)
            ipFockOp_r = ipFockOp(iShll_r)
            nPrim_r  = nExp(iShll_r)
            nCntrc_r = nBasis_Cntrct(iShll_r)
            iCmp_r = (iAng+1)*(iAng+2)/2
            If (Prjct(iShll_r)) iCmp_r = 2*iAng+1
*
*                                                                      *
************************************************************************
*                                                                      *
            If (ECP(iCnttp)) Then
#ifdef _DEBUG_
               If (lPP) Then
                  Write (6,*) 'Reference is ECP (Pseudo Potential)'
               Else
                  Write (6,*) 'Reference is ECP (Huzinaga type)'
               End If
               Call RecPrt('Reference Exponents',' ',
     &                     DInf(ipExp_r),1,nPrim_r)
               Call RecPrt('Reference Coefficients',' ',
     &                     DInf(ipCff_r),nPrim_r,nCntrc_r)
               Call RecPrt('Reference Fock operator',' ',
     &                     DInf(ipFockOp_r),nCntrc_r,nCntrc_r)
#endif
               Call OrbType(iAtmNr(nCnttp),List_AE,31)
               Call ECP_Shells(iAtmNr(iCnttp),List)
               If (lPP.or.nM1(iCnttp).eq.0) Then
*
*                 Pseud potential case
*
                  nRemove = List_AE(iAng)-List(iAng)
*
               Else
*
*                 Huzinaga type, remove according to the number of
*                 projected shells.
*
                  iAngMax_Proj=nPrj_Shells(iCnttp)
                  If (iAng.le.iAngMax_Proj) Then
                     iShll_Proj_r = ipPrj(iCnttp) + iAng
                     nCntrc_Proj = nBasis(iShll_Proj_r)
                     nRemove = nCntrc_Proj
                  Else
                     nRemove=0
                  End If
*
*                 If too many try the default
*
                  If (nRemove.gt.nCntrc_r) Then
                     nRemove = List_AE(iAng)-List(iAng)
                  End If
*
               End If ! lPP
#ifdef _DEBUG_
               Write(6,*) 'nRemove=',nRemove
               Write(6,*) 'List_Add(iAng)=',List_Add(iAng)
#endif
               nRemove = nRemove - List_Add(iAng)
#ifdef _DEBUG_
               Write(6,*) 'nRemove=',nRemove
#endif
               Test_Charge=Test_Charge + DBLE(2*(2*iAng+1)*nRemove)
*
*              Update pointers in case of ECP
*
*              Update the number of contracted functions of ref.
               nCntrc_t = nCntrc_r - nRemove
*              Update pointer to contraction coeffs of ref
               ipCff_r = ipCff_r + nRemove*nPrim_r
*              Pick up relevant parts of the FockOp matrix of ref.
               Call mma_allocate(FockOp_t,nCntrc_t**2)
               ipFockOp_t=1
               iOff_t = ipFockOp_t
               iOff_r = ipFockOp_r + nRemove*nCntrc_r + nRemove
               Do i = 1, nCntrc_t
                  call dcopy_(nCntrc_t,DInf(iOff_r),1,
     &                                FockOp_t(iOff_t),1)
                  iOff_r = iOff_r + nCntrc_r
                  iOff_t = iOff_t + nCntrc_t
               End Do
               nCntrc_r = nCntrc_t
            End If
*                                                                      *
************************************************************************
*                                                                      *
*
#ifdef _DEBUG_
            Call RecPrt('Actual Exponents',' ',
     &                  DInf(ipExp_a),1,nPrim_a)
            Call RecPrt('Actual Coefficients',' ',
     &                  DInf(ipCff_a),nPrim_a,nCntrc_a)
            Call RecPrt('Reference Exponents',' ',
     &                  DInf(ipExp_r),1,nPrim_r)
            Call RecPrt('Reference Coefficients',' ',
     &                  DInf(ipCff_r),nPrim_r,nCntrc_r)
            If (Allocated(FockOp_t)) Then
               Call RecPrt('Reference Fock operator',' ',
     &                     FockOp_t,nCntrc_r,nCntrc_r)
            Else
               Call RecPrt('Reference Fock operator',' ',
     &                     DInf(ipFockOp_r),nCntrc_r,nCntrc_r)
          End If
#endif
            If (Allocated(FockOp_t)) Then
               Check=DDot_(nCntrc_r**2,FockOp_t,1,
     &                                FockOp_t,1)
            Else
               Check=DDot_(nCntrc_r**2,DInf(ipFockOp_r),1,
     &                                DInf(ipFockOp_r),1)
            End If
            If (Check.eq.Zero) Go To 999
            If (Charge(iCnttp).eq.Zero) Go To 999
*                                                                      *
************************************************************************
*                                                                      *
            naa = nElem(iAng)*nElem(iAng)
            nScr1 = Max(nPrim_a,nPrim_r)*Max(nCntrc_a,nCntrc_r)*naa
            nScr2 = Max(nCntrc_a,nCntrc_r)**2*naa
*                                                                      *
************************************************************************
*                                                                      *
            ip = ipExp(iShll+1)
*                                                                      *
************************************************************************
*                                                                      *
*           Compute S_AA
*
            nOrdOp=0
            Call One_Int(MltPrm,DInf,nDInf,A,ip,Info,nInfo,jShll,iAng,
     &                   iComp,nOrdOp,nScr1,nScr2,naa,ipSAA,nSAA,
     &                   iShll_a,nPrim_a,ipExp_a,nCntrc_a,ipCff_a,
     &                   iCmp_a,
     &                   iShll_a,nPrim_a,ipExp_a,nCntrc_a,ipCff_a,
     &                   iCmp_a)
*                                                                      *
************************************************************************
*                                                                      *
*           Compute S_AR
*
            nOrdOp=0
            Call One_Int(MltPrm,DInf,nDInf,A,ip,Info,nInfo,jShll,iAng,
     &                   iComp,nOrdOp,nScr1,nScr2,naa,ipSAR,nSAR,
     &                   iShll_a,nPrim_a,ipExp_a,nCntrc_a,ipCff_a,
     &                   iCmp_a,
     &                   iShll_r,nPrim_r,ipExp_r,nCntrc_r,ipCff_r,
     &                   iCmp_r)
*
            nSRR = nCntrc_r**2 * naa
*                                                                      *
************************************************************************
*                                                                      *
*           Reorder and compute the inverse of SAA
*
            ipS_AA = ip
            ip = ip + nSAA
            Call Reorder(DInf(ipSAA),DInf(ipS_AA),
     &                   nCntrc_a,nCntrc_a,iCmp_a,iCmp_a)
#ifdef _DEBUG_
            Call RecPrt('Reordered SAA',' ',DInf(ipS_AA),
     &                  nCntrc_a*iCmp_a,nCntrc_a*iCmp_a)
#endif
            Call MInv(DInf(ipS_AA),DInf(ipSAA),iSing,D,nCntrc_a*iCmp_a)
            ip = ip -nSAA
#ifdef _DEBUG_
            Write (6,*) 'iSing=',iSing
            Write (6,*) 'Det=',D
            Call RecPrt('Inverse of SAA',' ',DInf(ipSAA),
     &                  nCntrc_a*iCmp_a,nCntrc_a*iCmp_a)
#endif
*
*           Reorder SAR
            ipS_AR = ip
            ip = ip + nSAR
            Call Reorder(DInf(ipSAR),DInf(ipS_AR),
     &                   nCntrc_a,nCntrc_r,iCmp_a,iCmp_r)
#ifdef _DEBUG_
            Call RecPrt('Reordered SAR',' ',DInf(ipS_AR),
     &                  nCntrc_a*iCmp_a,nCntrc_r*iCmp_r)
#endif
*
*           Expand and reorder the reference fock operator
*
            ipE_R=ip
            ip = ip + nSRR
            ipTmp = ip
            ip = ip + nSRR
            Call FZero(DInf(ipTmp),nSRR)
            If (Allocated(FockOp_t)) Then
               Do iB = 1, nCntrc_r
                  Do jB = 1, nCntrc_r
                     ijB=(jB-1)*nCntrc_r+iB
                     iFrom = ipFockOp_t-1 + (jB-1)*nCntrc_r+iB
                     Temp = FockOp_t(iFrom)
                     Do iC = 1, iCmp_r
                        ijC=(iC-1)*iCmp_r+iC
                        iTo = ipTmp-1 + (ijC-1)*nCntrc_r**2+ijB
                        DInf(iTo) = Temp
                     End Do
                  End Do
               End Do
            Else
               Do iB = 1, nCntrc_r
                  Do jB = 1, nCntrc_r
                     ijB=(jB-1)*nCntrc_r+iB
                     iFrom = ipFockOp_r-1 + (jB-1)*nCntrc_r+iB
                     Temp = DInf(iFrom)
                     Do iC = 1, iCmp_r
                        ijC=(iC-1)*iCmp_r+iC
                        iTo = ipTmp-1 + (ijC-1)*nCntrc_r**2+ijB
                        DInf(iTo) = Temp
                     End Do
                  End Do
               End Do
            End If
#ifdef _DEBUG_
            Call RecPrt('Expanded ER',' ',DInf(ipTmp),
     &                  nCntrc_r*nCntrc_r,iCmp_r*iCmp_r)
#endif
            Call Reorder(DInf(ipTmp),DInf(ipE_R),
     &                   nCntrc_r,nCntrc_r,iCmp_r,iCmp_r)
            ip = ip - nSRR ! Release ipTmp
#ifdef _DEBUG_
            Call RecPrt('Reordered ER',' ',DInf(ipE_R),
     &                  nCntrc_r*iCmp_r,nCntrc_r*iCmp_r)
#endif
*
*           Form (SAA)-1 SAR
*
            ipTmp1 = ip
            ip = ip + nSAR
            Call DGEMM_('N','N',
     &                  nCntrc_a*iCmp_a,nCntrc_r*iCmp_r,nCntrc_a*iCmp_a,
     &                  1.0d0,DInf(ipSAA),nCntrc_a*iCmp_a,
     &                  DInf(ipS_AR),nCntrc_a*iCmp_a,
     &                  0.0d0,DInf(ipTmp1),nCntrc_a*iCmp_a)
#ifdef _DEBUG_
            Call RecPrt('(SAA)^-1 SAR',' ',DInf(ipTmp1),
     &                  nCntrc_a*iCmp_a,nCntrc_r*iCmp_r)
#endif
*
*           Form (SAA)-1 SAR ER
*
            ipTmp2 = ip
            ip = ip + nSAR
            Call DGEMM_('N','N',
     &                  nCntrc_a*iCmp_a,nCntrc_r*iCmp_r,nCntrc_r*iCmp_r,
     &                  1.0d0,DInf(ipTmp1),nCntrc_a*iCmp_a,
     &                  DInf(ipE_R),nCntrc_r*iCmp_r,
     &                  0.0d0,DInf(ipTmp2),nCntrc_a*iCmp_a)
#ifdef _DEBUG_
            Call RecPrt('(SAA)^-1 SAR ER',' ',DInf(ipTmp2),
     &                  nCntrc_a*iCmp_a,nCntrc_r*iCmp_r)
#endif
*
*           Form (SAA)-1 SAR ER (SAR)^T (SAA)-1
*
            Call DGEMM_('N','T',
     &                  nCntrc_a*iCmp_a,nCntrc_a*iCmp_a,nCntrc_r*iCmp_r,
     &                  1.0d0,DInf(ipTmp2),nCntrc_a*iCmp_a,
     &                  DInf(ipTmp1),nCntrc_a*iCmp_a,
     &                  0.0d0,DInf(ipSAA),nCntrc_a*iCmp_a)
#ifdef _DEBUG_
            Call RecPrt('EA',' ',DInf(ipSAA),
     &                  nCntrc_a*iCmp_a,nCntrc_a*iCmp_a)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*           Now we just need to reorder and put it into place!
*
            ipTmp3 = ip
            ip = ip + nSRR
            Call Reorder(DInf(ipSAA),DInf(ipTmp3),
     &                   nCntrc_a,iCmp_a,nCntrc_a,iCmp_a)
#ifdef _DEBUG_
            Call RecPrt('Reordered EA',' ',DInf(ipTmp3),
     &                  nCntrc_a*nCntrc_a,iCmp_a*iCmp_a)
#endif
*
            Do iB = 1, nCntrc_a
               Do jB = 1, nCntrc_a
                  ijB=(jB-1)*nCntrc_a+iB
                  iTo   = ipFockOp_a-1 + (jB-1)*nCntrc_a+iB
                  DInf(iTo) = Zero
                  Do iC = 1, iCmp_a
                     ijC=(iC-1)*iCmp_a+iC
                     iFrom = ipTmp3-1 + (ijC-1)*nCntrc_a**2+ijB
                     DInf(iTo) = DInf(iTo) + DInf(iFrom)
                  End Do
                  DInf(iTo) = DInf(iTo)/DBLE(iCmp_a)
               End Do
            End Do
 999        Continue
            If (Allocated(FockOp_t)) Call mma_deallocate(FockOp_t)
#ifdef _DEBUG_
            Call RecPrt('Actual Fock operator',' ',DInf(ipFockOp_a),
     &                  nCntrc_a,nCntrc_a)
#endif
*                                                                      *
************************************************************************
*                                                                      *
         End Do  ! iAng
*
         Charge_Actual=DBLE(iAtmNr(iCnttp))
         Charge_Effective=Charge(iCnttp)
         qTest=Test_Charge -
     &         (Charge_Actual-Charge_Effective)
c         write(6,*)'qtest, Test_Charge = ',qtest, Test_Charge
c         write(6,*)'Charge_Actual,Charge_Effective = ',
c     &               Charge_Actual,Charge_Effective
         If (qTest.eq.Zero.or.Charge(iCnttp).eq.Zero) Then
            FockOp(iCnttp)=.TRUE.
         Else If (Try_Again) Then
            If (qTest.eq.2.0D0) Then
*              s
               List_Add(0)=1
            Else If (qTest.eq.6.0D0) Then
*              p
               List_Add(1)=1
            Else If (qTest.eq.10.0D0) Then
*              d
               List_Add(2)=1
            Else If (qTest.eq.8.0D0) Then
*              s,p
               List_Add(0)=1
               List_Add(1)=1
            Else If (qTest.eq.12.0D0) Then
*              s,d
               List_Add(0)=1
               List_Add(2)=1
            Else If (qTest.eq.16.0D0) Then
*              p,d
               List_Add(1)=1
               List_Add(2)=1
            Else If (qTest.eq.18.0D0) Then
*              s,p,d
               List_Add(0)=1
               List_Add(1)=1
               List_Add(2)=1
            Else If (qTest.eq.26.0D0) Then
*              2s,2p,d
               List_Add(0)=2
               List_Add(1)=2
               List_Add(2)=1
            End If
            Try_Again=.False.
            Go To 777
         Else
            Write (6,*) 'GuessOrb option turned off!'
            FockOp(iCnttp)=.FALSE.
         End If
*                                                                      *
************************************************************************
*                                                                      *
 1000 Continue
*
*                                                                      *
************************************************************************
*                                                                      *
*     Restore the correct nCnttp value
*
      nCnttp=mCnttp
*
#ifdef _INSANE_DEBUG_
      nPrint(113)=5
      nPrint(114)=5
      nPrint(116)=5
      nPrint(122)=5
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Check if we can activate the computation of FckInt!
*
      Do_FckInt=.True.
      Do 2000 iCnttp = 1, nCnttp
         If(AuxCnttp(iCnttp) .or.
     &      FragCnttp(iCnttp) .or.
     &      nFragType(iCnttp).gt.0 .or.
     &      FockOp(iCnttp)) Then
           Goto 2000
         End If
*
         Do_FckInt = Do_FckInt .and. FockOp(iCnttp) ! To be activated!
*
 2000 Continue
*                                                                      *
************************************************************************
*                                                                      *
      Call qExit('Fix_FockOp')
      Return
      End
