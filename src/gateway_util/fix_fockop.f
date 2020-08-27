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
* Copyright (C) 2017,2020, Roland Lindh                                *
************************************************************************
      Subroutine Fix_FockOp(LuRd)
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
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "stdalloc.fh"
#include "real.fh"
#include "print.fh"
#include "status.fh"
#include "periodic_table.fh"
      External MltPrm, KnEPrm, NAPrm
      Real*8, Allocatable :: FockOp_t(:)
      Real*8, Allocatable :: Scr1(:), Scr2(:), Scr3(:)
      Real*8, Allocatable :: S12i(:,:), EVec(:,:), EVal(:)
      Real*8, Allocatable :: FPrim(:,:), Temp(:,:), C(:,:)
      Real*8, Allocatable :: Hm1(:,:), Ovr(:,:)
      Real*8, Allocatable :: S_AA(:), S_AR(:), E_R(:)
      Real*8, Allocatable :: Tmp1(:), Tmp2(:), Tmp3(:)
      Real*8, Allocatable :: KnE(:), NAE(:), Ovrlp(:)
      Real*8, Allocatable :: SAA(:), SAR(:)
      Character*13 DefNm
      Character*80 Ref(2), Bsl_, BSLbl
      Character *256 Basis_lib, Fname
      Character*180, Allocatable :: STDINP(:) ! CGGn
      Integer BasisTypes(4), nDel(MxAng)
      Integer List_AE(0:iTabMx), List(0:iTabMx), List_Add(0:iTabMx)
      Logical Try_Again
      Real*8 A(4)
      Data DefNm/'basis_library'/
*                                                                      *
************************************************************************
*                                                                      *
#include "getbs_interface.fh"
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*                                                                      *
************************************************************************
*                                                                      *
!#define _DEBUG_
#ifdef _DEBUG_
      nPrint(114)=99
      nPrint(116)=99
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Generate a dummy center. This is fine since we will only do
*     1-center overlap integrals here.
*
      A(:)=Zero
      Call mma_allocate(STDINP,mxAtom*2,label='STDINP')
*
      nOrdOp=2
      iComp = 1
*
      nPrp=Max(4,nMltpl)
      nDiff = 0
*
      List   (:)=0
      List_AE(:)=0
      BasisTypes(:)=0
      lSTDINP=0
*
*     Loop over all valence shell with a non-funtional FockOp
*
      mCnttp = nCnttp   ! to be restored at the end
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Do iCnttp = 1, mCnttp
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
         iFerm=1
         If (fMass(iCnttp).ne.1.0D0) iFerm=2
*
         If (dbsc(iCnttp)%FOp.and.Charge(iCnttp).eq.0.0D0) Then
            Do iAng = 0, dbsc(iCnttp)%nVal-1
               iShll_a    = dbsc(iCnttp)%iVal + iAng
               Shells(iShll_a)%FockOp(:,:)=Zero
            End Do
         End If
*
         If(dbsc(iCnttp)%Aux .or.
     &      dbsc(iCnttp)%Frag .or.
     &      dbsc(iCnttp)%nFragType.gt.0 .or.
     &      dbsc(iCnttp)%FOp) Then
           Cycle
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
            Do iAng = 0, dbsc(iCnttp)%nVal-1
*
               iShll_a    = dbsc(iCnttp)%iVal + iAng
               nPrim_a  = Shells(iShll_a)%nExp
               If (nPrim_a.eq.0) Cycle
               nCntrc_a = Shells(iShll_a)%nBasis_C
               iCmp_a = (iAng+1)*(iAng+2)/2
               If (Shells(iShll_a)%Prjct) iCmp_a = 2*iAng+1
               naa = nElem(iAng)*nElem(iAng)
               nScr1 = Max(nPrim_a,nPrim_a)*Max(nCntrc_a,nCntrc_a)*naa
               nScr2 = Max(nCntrc_a,nCntrc_a)**2*naa
               Call mma_allocate(Scr1,nScr1,Label='Scr1')
               Call mma_allocate(Scr2,nScr2,Label='Scr2')
*                                                                      *
************************************************************************
*                                                                      *
*              Compute the kinetic integrals
*
               nOrdOp=2
               nSAA=nCntrc_a**2 * naa
*
               Call KnEMmP(nHer,MmKnEP,iAng,iAng,nOrdOp)
               nScr3=nPrim_a**2 * MmKnEP
               Call mma_allocate(Scr3,nScr3,Label='Scr1')
*
               Call mma_Allocate(KnE,NSAA,Label='KnE')
               Call One_Int(KnEPrm,Scr3,nScr3,A,iAng,iComp,nOrdOp,
     &                      Scr1,nScr1,Scr2,nScr2,naa,KnE,nSAA,
     &                      iShll_a,nPrim_a,Shells(iShll_a)%Exp,
     &                     nCntrc_a,Shells(iShll_a)%Cff_c(1,1,1),iCmp_a,
     &                      iShll_a,nPrim_a,Shells(iShll_a)%Exp,
     &                     nCntrc_a,Shells(iShll_a)%Cff_c(1,1,1),iCmp_a)
               Call mma_deallocate(Scr3)
*define _DEBUG_
#ifdef _DEBUG_
               Call DScal_(nCntrc_a**2*iCmp_a**2,
     &                     xFactor,KnE,1)
               Call RecPrt('Kinetric Energy Integrals',' ',
     &                     KnE,nCntrc_a**2,iCmp_a**2)
               Call DScal_(nCntrc_a**2*iCmp_a**2,
     &                     1.0D0/xFactor,KnE,1)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*              Compute the nuclear-attraction integrals
*
               nOrdOp=0
               A(4) = DBLE(iCnttp) ! Dirty tweak
               nSBB=nCntrc_a**2 * naa
               Call mma_Allocate(NAE,nSBB,Label='NAE')
*
               Call NAMem(nHer,MemNA ,iAng,iAng,nOrdOp)
               nScr3=nPrim_a**2 * MemNA
               Call mma_allocate(Scr3,nScr3,Label='Scr3')
*
               Call One_Int(NAPrm,Scr3,nScr3,A,iAng,iComp,nOrdOp,
     &                      Scr1,nScr1,Scr2,nScr2,naa,NAE,nSBB,
     &                      iShll_a,nPrim_a,Shells(iShll_a)%Exp,
     &                     nCntrc_a,Shells(iShll_a)%Cff_c(1,1,1),iCmp_a,
     &                      iShll_a,nPrim_a,Shells(iShll_a)%Exp,
     &                     nCntrc_a,Shells(iShll_a)%Cff_c(1,1,1),iCmp_a)
               Call mma_deallocate(Scr3)
#ifdef _DEBUG_
               Call RecPrt('Nuclear-attraction Integrals',' ',
     &                     NAE,nCntrc_a**2,iCmp_a**2)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*              Add together the kinetic and nuclear-attraction
*
               Call DaXpY_(nCntrc_a**2*iCmp_a**2,
     &                     xFactor,
     &                     KnE,1,
     &                     NAE,1)
               Call mma_deallocate(KnE)
*
*              Change to proper order (nCntrc_a * iCmp_a)
*
               Call mma_allocate(Hm1,nCntrc_a**2,iCmp_a**2,Label='Hm1')
               Call Reorder_GW(NAE,Hm1,nCntrc_a,nCntrc_a,iCmp_a,iCmp_a)
*                                                                      *
************************************************************************
*                                                                      *
*              Compute the overlap integrals
*
               nOrdOp=0
               nSCC=nCntrc_a**2 * naa
               Call mma_allocate(Ovrlp,nSCC,Label='Ovrlp')
*
               Call MltMmP(nHer,MmMltp,iAng,iAng,nOrdOp)
               nScr3=nPrim_a**2 * MmMltp
               Call mma_allocate(Scr3,nScr3,Label='Scr3')
*
               Call One_Int(MltPrm,Scr3,nScr3,A,iAng,iComp,nOrdOp,
     &                      Scr1,nScr1,Scr2,nScr2,naa,Ovrlp,nSCC,
     &                      iShll_a,nPrim_a,Shells(iShll_a)%Exp,
     &                     nCntrc_a,Shells(iShll_a)%Cff_c(1,1,1),iCmp_a,
     &                      iShll_a,nPrim_a,Shells(iShll_a)%Exp,
     &                     nCntrc_a,Shells(iShll_a)%Cff_c(1,1,1),iCmp_a)
               Call mma_deallocate(Scr3)
#ifdef _DEBUG_
               Call RecPrt('Overlap Integrals',' ',
     &                     Ovrlp,nCntrc_a**2,iCmp_a**2)
#endif
*
*              Change to proper order (nCntrc_a * iCmp_a)
*
               nBF = nCntrc_a*iCmp_a
               Call mma_allocate(Ovr,nBF,nBF,Label='Ovr')
               Call Reorder_GW(Ovrlp,Ovr,
     &                      nCntrc_a,nCntrc_a,iCmp_a,iCmp_a)
               Call mma_deallocate(Ovrlp)
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
               Call mma_Allocate(S12i,nBF,nBF,Label='S')
               S12i(:,:)=Zero
*
*              1) Compute the eigenvectors and eigenvalues of the
*                 overlap matrix
*
               Call mma_allocate(EVal,nBF*(nBF+1)/2,Label='EVal')
               Call mma_allocate(EVec,nBF,nBF,Label='EVec')
               EVec(:,:)=Zero
               Do iBF = 1, nBF
                  EVec(iBF,iBF)=One
                  Do jBF = 1, iBF
                     ijTri = (iBF-1)*iBF/2 + jBF
                     EVal(ijTri) = Ovr(iBF,jBF)
                  End Do
               End Do
               Call mma_deallocate(Ovr)
               Call NIDiag_new(EVal,EVec,nBF,nBF,0)
*
*              2) Construct S^(1/2) and S^(-1/2)
*
               Do kEval = 1, nBF
                  e   = EVal(kEval*(kEval+1)/2)
                  e12i= 1.0D0/Sqrt(e)
                  Do iBF = 1, nBF
                     C_ik = EVec(iBF,kEVal)
                     Do jBF = 1, nBF
                        C_jk = EVec(jBF,kEVal)
                        S12i(iBF,jBF) = S12i(iBF,jBF)
     &                                + C_ik * e12i * C_jk
                     End Do
                  End Do
               End Do
*
*              3) Form F' =  S^(-1/2) F S^(-1/2)
*
               Call mma_allocate(FPrim,nBF,nBF,Label='FPrim')
               FPrim(:,:)=Zero
               Call mma_allocate(Temp,nBF,nBF,Label='Temp')
               Call DGEMM_('N','N',
     &                     nBF,nBF,nBF,
     &                     1.0d0,S12i,nBF,
     &                           Hm1,nBF,
     &                     0.0d0,Temp,nBF)
               Call DGEMM_('N','N',
     &                     nBF,nBF,nBF,
     &                     1.0d0,Temp,nBF,
     &                           S12i,nBF,
     &                     0.0d0,FPrim,nBF)
*
*              4) Compute C' and the eigenvalues
*
               EVec(:,:)=Zero
               Do iBF = 1, nBF
                  EVec(iBF,iBF)=One
                  Do jBF = 1, iBF
                     ij    = (jBF-1)*nBF   + iBF
                     ijTri = (iBF-1)*iBF/2 + jBF
                     EVal(ijTri) = FPrim(iBF,jBF)
                  End Do
               End Do
               Call mma_deallocate(Temp)
               Call mma_deallocate(FPrim)
               Call NIDiag_new(EVal,EVec,nBF,nBF,0)
*
*              5) Form C = S^(-1/2) C'
*
               Call mma_allocate(C,nBF,nBF,Label='C')
               C(:,:)=Zero
               Call DGEMM_('N','N',
     &                     nBF,nBF,nBF,
     &                     1.0d0,S12i,nBF,
     &                           EVec,nBF,
     &                     0.0d0,C,nBF)
#ifdef _DEBUG_
      Call RecPrt('Cs for F',' ',C,nBF,nBF)
#endif
*
*              6) Form the matrix representation of the Fock operator
*
               Call mma_deallocate(Hm1)
               Call mma_allocate(Hm1,nBF,nBF,Label='Hm1')
               Hm1(:,:)=Zero
               Do kEval = 1, nBF
                  e   = EVal(kEval*(kEval+1)/2)
                  Do iBF = 1, nBF
                     C_ik = C(iBF,kEVal)
                     Do jBF = 1, nBF
                        C_jk = C(jBF,kEVal)
                        ij = (jBF-1)*nBF + iBF
                        Hm1(iBF,jBF) = Hm1(iBF,jBF)
     &                                      + C_ik * e * C_jk
                     End Do
                  End Do
               End Do
               Call mma_deallocate(C)
*
               Call Reorder_GW(Hm1,NAE,
     &                      nCntrc_a,iCmp_a,nCntrc_a,iCmp_a)
               Call mma_deallocate(Hm1)
*
*              Make result isotropic and distribute
*
               Do iB = 1, nCntrc_a
                  Do jB = 1, nCntrc_a
                     ijB=(jB-1)*nCntrc_a+iB
                     Tmp = Zero
                     Do iC = 1, iCmp_a
                        ijC=(iC-1)*iCmp_a+iC
                        iFrom = (ijC-1)*nCntrc_a**2+ijB
                        Tmp = Tmp + NAE(iFrom)
                     End Do
                     Shells(iShll_a)%FockOp(iB,jB) = Tmp/DBLE(iCmp_a)
                  End Do
               End Do
               Call mma_deallocate(NAE)
#ifdef _DEBUG_
               Call RecPrt('Actual Fock operator',' ',
     &                     Shells(iShll_a)%FockOp,nCntrc_a,nCntrc_a)
#endif
               Call mma_deallocate(EVal)
               Call mma_deallocate(EVec)
               Call mma_deallocate(S12i)
               Call mma_deallocate(Scr1)
               Call mma_deallocate(Scr2)
            End Do
*
            dbsc(iCnttp)%FOp=.TRUE.
            Cycle
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
         Call GetBS(Fname,Bsl_,Indx-1,lAng,iShll,MxAng,
     &              Charge(nCnttp),iAtmNr(nCnttp),BLine,Ref,
     &              PAM2(nCnttp),NoPairL(nCnttp),SODK(nCnttp),
     &              CrRep(nCnttp),nProj,nAIMP,UnNorm,nDel,
     &              nVal,   nPrj,   nSRO,   nSOC,  nPP,
     &              ipVal_, ipPrj_, ipSRO_, ipSOC_,ipPP_,
     &              LuRd,BasisTypes,
     &              STDINP,lSTDINP,.False.,.true.,' ')
*
         If (.Not.dbsc(nCnttp)%FOp) Then
            Write (6,*) 'Fix_FockOp: reference basis doesn''t contain a'
     &                //' proper Fock operator'
            Cycle
         End If
         Shells(jShll+1)%Transf=.False.
         Shells(jShll+1)%Prjct =.False.
         Shells(jShll+2)%Transf=.False.
         Shells(jShll+2)%Prjct =.False.
         dbsc(nCnttp)%iVal = ipVal_
         dbsc(nCnttp)%iPrj = ipPrj_
         dbsc(nCnttp)%iSRO = ipSRO_
         dbsc(nCnttp)%iSOC = ipSOC_
         ipPP(nCnttp)  = ipPP_
         dbsc(nCnttp)%nVal = nVal
         dbsc(nCnttp)%nPrj = nPrj
         dbsc(nCnttp)%nSRO = nSRO
         dbsc(nCnttp)%nSOC = nSOC
         nPP_Shells(nCnttp)  = nPP
         nTot_Shells(nCnttp) = nVal+nPrj+nSRO+nSOC+nPP
*                                                                      *
************************************************************************
*                                                                      *
*        Start processing shells of iCnttp and mCnttp. Loop only over
*        the shells of iCnttp (mCnttp might be larger!)
*
         Try_Again=.True.
         Call ICopy(1+iTabMx,[0],0,List_Add,1)
 777     Continue
*
         Test_Charge=Zero
         Do iAng = 0, dbsc(iCnttp)%nVal-1
*
*           Pointers to the actuall shell
*
            iShll_a    = dbsc(iCnttp)%iVal + iAng
            nPrim_a  = Shells(iShll_a)%nExp
            If (nPrim_a.eq.0) Cycle
            nCntrc_a = Shells(iShll_a)%nBasis_C
            iCmp_a = (iAng+1)*(iAng+2)/2
            If (Shells(iShll_a)%Prjct) iCmp_a = 2*iAng+1
*
*           Pointers to the reference shell
*
            iShll_r = dbsc(nCnttp)%iVal + iAng
            nPrim_r  = Shells(iShll_r)%nExp
            If (nPrim_r.eq.0) Then
               Write (6,*) 'GuessOrb option turned off!'
               dbsc(iCnttp)%FOp=.FALSE.
               Exit
            End If
            nCntrc_r = Shells(iShll_r)%nBasis_C
            iCmp_r = (iAng+1)*(iAng+2)/2
            If (Shells(iShll_r)%Prjct) iCmp_r = 2*iAng+1
*
*                                                                      *
************************************************************************
*                                                                      *
            If (dbsc(iCnttp)%ECP) Then
#ifdef _DEBUG_
               If (lPP) Then
                  Write (6,*) 'Reference is ECP (Pseudo Potential)'
               Else
                  Write (6,*) 'Reference is ECP (Huzinaga type)'
               End If
               Call RecPrt('Reference Exponents',' ',
     &                    Shells(iShll_r)%Exp,1,nPrim_r)
               Call RecPrt('Reference Coefficients',' ',
     &                    Shells(iShll_r)%Cff_c(1,1,1),nPrim_r,nCntrc_r)
               Call RecPrt('Reference Fock operator',' ',
     &                    Shells(iShll_r)%FockOp,nCntrc_r,nCntrc_r)
#endif
               Call OrbType(iAtmNr(nCnttp),List_AE,31)
               Call ECP_Shells(iAtmNr(iCnttp),List)
               If (lPP.or.dbsc(iCnttp)%nM1.eq.0) Then
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
                  iAngMax_Proj=dbsc(iCnttp)%nPrj
                  If (iAng.le.iAngMax_Proj) Then
                     iShll_Proj_r = dbsc(iCnttp)%iPrj + iAng
                     nCntrc_Proj = Shells(iShll_Proj_r)%nBasis
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
*              Pick up relevant parts of the FockOp matrix of ref.
               Call mma_allocate(FockOp_t,nCntrc_t**2)
               ipFockOp_t=1
               iOff_t = ipFockOp_t
               Do i = 1, nCntrc_t
                  call dcopy_(nCntrc_t,
     &                    Shells(iShll_r)%FockOp(nRemove+1,nRemove+i),1,
     &                    FockOp_t(iOff_t),1)
                  iOff_t = iOff_t + nCntrc_t
               End Do
               nCntrc_r = nCntrc_t
            Else
               nRemove = 0
            End If
*                                                                      *
************************************************************************
*                                                                      *
*
#ifdef _DEBUG_
            Call RecPrt('Actual Exponents',' ',
     &                 Shells(iShll_a)%Exp,1,nPrim_a)
            Call RecPrt('Actual Coefficients',' ',
     &                 Shells(iShll_a)%Cff_c(1,1,1),nPrim_a,nCntrc_a)
            Call RecPrt('Reference Exponents',' ',
     &                 Shells(iShll_r)%Exp,1,nPrim_r)
            Call RecPrt('Reference Coefficients',' ',
     &                 Shells(iShll_r)%Cff_c(1,nRemove+1,1),
     &                                              nPrim_r,nCntrc_r)
            If (Allocated(FockOp_t)) Then
               Call RecPrt('Reference Fock operator',' ',
     &                     FockOp_t,nCntrc_r,nCntrc_r)
            Else
               Call RecPrt('Reference Fock operator',' ',
     &                     Shells(iShll_r)%FockOp,nCntrc_r,nCntrc_r)
            End If
#endif
            If (Allocated(FockOp_t)) Then
               Check=DDot_(nCntrc_r**2,FockOp_t,1,
     &                                FockOp_t,1)
            Else
               Check=DDot_(nCntrc_r**2,Shells(iShll_r)%FockOp,1,
     &                                 Shells(iShll_r)%FockOp,1)
            End If
            If (Check.eq.Zero .or.  Charge(iCnttp).eq.Zero) Then
               If (Allocated(FockOp_t)) Call mma_deallocate(FockOp_t)
               Cycle
            End If
*                                                                      *
************************************************************************
*                                                                      *
            naa = nElem(iAng)*nElem(iAng)
            nScr1 = Max(nPrim_a,nPrim_r)*Max(nCntrc_a,nCntrc_r)*naa
            nScr2 = Max(nCntrc_a,nCntrc_r)**2*naa
            Call mma_allocate(Scr1,nScr1,Label='Scr1')
            Call mma_allocate(Scr2,nScr2,Label='Scr2')
*                                                                      *
************************************************************************
*                                                                      *
*           Compute S_AA
*
            nOrdOp=0
            nSAA= nCntrc_a**2 * naa
            Call mma_allocate(SAA,nSAA,Label='SAA')
*
            Call MltMmP(nHer,MmMltp,iAng,iAng,nOrdOp)
            nScr3=nPrim_a**2 * MmMltp
            Call mma_allocate(Scr3,nScr3,Label='Scr3')
*
            Call One_Int(MltPrm,Scr3,nScr3,A,iAng,iComp,nOrdOp,
     &                   Scr1,nScr1,Scr2,nScr2,naa,SAA,nSAA,
     &                   iShll_a,nPrim_a,Shells(iShll_a)%Exp,
     &                   nCntrc_a,Shells(iShll_a)%Cff_c(1,1,1),iCmp_a,
     &                   iShll_a,nPrim_a,Shells(iShll_a)%Exp,
     &                   nCntrc_a,Shells(iShll_a)%Cff_c(1,1,1),iCmp_a)
            Call mma_deallocate(Scr3)
*                                                                      *
************************************************************************
*                                                                      *
*           Compute S_AR
*
            nOrdOp=0
            nSAR=nCntrc_a*nCntrc_r * naa
            Call mma_allocate(SAR,nSAR,Label='SAR')
*
            Call MltMmP(nHer,MmMltp,iAng,iAng,nOrdOp)
            nScr3=nPrim_a*nPrim_r * MmMltp
            Call mma_allocate(Scr3,nScr3,Label='Scr3')
*
            Call One_Int(MltPrm,Scr3,nScr3,A,iAng,iComp,nOrdOp,
     &                   Scr1,nScr1,SCr2,nScr2,naa,SAR,nSAR,
     &                   iShll_a,nPrim_a,Shells(iShll_a)%Exp,
     &                   nCntrc_a,Shells(iShll_a)%Cff_c(1,1,1),iCmp_a,
     &                   iShll_r,nPrim_r,Shells(iShll_r)%Exp,
     &                   nCntrc_r,Shells(iShll_r)%Cff_c(1,1+nRemove,1),
     &                                                         iCmp_a)
            Call mma_deallocate(Scr3)
*
            nSRR = nCntrc_r**2 * naa
*                                                                      *
************************************************************************
*                                                                      *
*           Reorder and compute the inverse of SAA
*
            Call mma_allocate(S_AA,nSAA,Label='S_AA')
            Call Reorder_GW(SAA,S_AA,
     &                   nCntrc_a,nCntrc_a,iCmp_a,iCmp_a)
#ifdef _DEBUG_
            Call RecPrt('Reordered SAA',' ',S_AA,
     &                  nCntrc_a*iCmp_a,nCntrc_a*iCmp_a)
#endif
            Call MInv(S_AA,SAA,iSing,D,nCntrc_a*iCmp_a)
            Call mma_deallocate(S_AA)
#ifdef _DEBUG_
            Write (6,*) 'iSing=',iSing
            Write (6,*) 'Det=',D
            Call RecPrt('Inverse of SAA',' ',SAA,
     &                  nCntrc_a*iCmp_a,nCntrc_a*iCmp_a)
#endif
*
*           Reorder SAR
            Call mma_allocate(S_AR,nSAR,Label='S_AR')
            Call Reorder_GW(SAR,S_AR,
     &                   nCntrc_a,nCntrc_r,iCmp_a,iCmp_r)
            Call mma_deallocate(SAR)
#ifdef _DEBUG_
            Call RecPrt('Reordered SAR',' ',S_AR,
     &                  nCntrc_a*iCmp_a,nCntrc_r*iCmp_r)
#endif
*
*           Expand and reorder the reference fock operator
*
            Call mma_allocate(E_R,nSRR,Label='E_R')
            Call mma_allocate(Tmp1,nSRR,Label='Tmp1')
            Tmp1(:)=Zero
            If (Allocated(FockOp_t)) Then
               Do iB = 1, nCntrc_r
                  Do jB = 1, nCntrc_r
                     ijB=(jB-1)*nCntrc_r+iB
                     iFrom = ipFockOp_t-1 + (jB-1)*nCntrc_r+iB
                     Tmp = FockOp_t(iFrom)
                     Do iC = 1, iCmp_r
                        ijC=(iC-1)*iCmp_r+iC
                        iTo = (ijC-1)*nCntrc_r**2+ijB
                        Tmp1(iTo) = Tmp
                     End Do
                  End Do
               End Do
            Else
               Do iB = 1, nCntrc_r
                  Do jB = 1, nCntrc_r
                     ijB=(jB-1)*nCntrc_r+iB
                     Tmp = Shells(iShll_r)%FockOp(iB,jB)
                     Do iC = 1, iCmp_r
                        ijC=(iC-1)*iCmp_r+iC
                        iTo = (ijC-1)*nCntrc_r**2+ijB
                        Tmp1(iTo) = Tmp
                     End Do
                  End Do
               End Do
            End If
#ifdef _DEBUG_
            Call RecPrt('Expanded ER',' ',Tmp1,
     &                  nCntrc_r*nCntrc_r,iCmp_r*iCmp_r)
#endif
            Call Reorder_GW(Tmp1,E_R,
     &                   nCntrc_r,nCntrc_r,iCmp_r,iCmp_r)
#ifdef _DEBUG_
            Call RecPrt('Reordered ER',' ',E_R,
     &                  nCntrc_r*iCmp_r,nCntrc_r*iCmp_r)
#endif
            Call mma_deallocate(Tmp1)
*
*           Form (SAA)-1 SAR
*
            Call mma_allocate(Tmp1,nSAR,Label='Tmp1')
            Call DGEMM_('N','N',
     &                  nCntrc_a*iCmp_a,nCntrc_r*iCmp_r,nCntrc_a*iCmp_a,
     &                  1.0d0,SAA,nCntrc_a*iCmp_a,
     &                        S_AR,nCntrc_a*iCmp_a,
     &                  0.0d0,Tmp1,nCntrc_a*iCmp_a)
#ifdef _DEBUG_
            Call RecPrt('(SAA)^-1 SAR',' ',Tmp1,
     &                  nCntrc_a*iCmp_a,nCntrc_r*iCmp_r)
#endif
            Call mma_deallocate(S_AR)
*
*           Form (SAA)-1 SAR ER
*
            Call mma_allocate(Tmp2,nSAR,Label='Tmp2')
            Call DGEMM_('N','N',
     &                  nCntrc_a*iCmp_a,nCntrc_r*iCmp_r,nCntrc_r*iCmp_r,
     &                  1.0d0,Tmp1,nCntrc_a*iCmp_a,
     &                        E_R,nCntrc_r*iCmp_r,
     &                  0.0d0,Tmp2,nCntrc_a*iCmp_a)
#ifdef _DEBUG_
            Call RecPrt('(SAA)^-1 SAR ER',' ',Tmp2,
     &                  nCntrc_a*iCmp_a,nCntrc_r*iCmp_r)
#endif
            Call mma_deallocate(E_R)
*
*           Form (SAA)-1 SAR ER (SAR)^T (SAA)-1
*
            Call DGEMM_('N','T',
     &                  nCntrc_a*iCmp_a,nCntrc_a*iCmp_a,nCntrc_r*iCmp_r,
     &                  1.0d0,Tmp2,nCntrc_a*iCmp_a,
     &                        Tmp1,nCntrc_a*iCmp_a,
     &                  0.0d0,SAA,nCntrc_a*iCmp_a)
#ifdef _DEBUG_
            Call RecPrt('EA',' ',SAA,
     &                  nCntrc_a*iCmp_a,nCntrc_a*iCmp_a)
#endif
            Call mma_deallocate(Tmp2)
            Call mma_deallocate(Tmp1)
*                                                                      *
************************************************************************
*                                                                      *
*           Now we just need to reorder and put it into place!
*
            Call mma_allocate(Tmp3,nSAA,Label='Tmp3')
            Call Reorder_GW(SAA,Tmp3,
     &                   nCntrc_a,iCmp_a,nCntrc_a,iCmp_a)
            Call mma_deallocate(SAA)
#ifdef _DEBUG_
            Call RecPrt('Reordered EA',' ',Tmp3,
     &                  nCntrc_a*nCntrc_a,iCmp_a*iCmp_a)
#endif
*
            Do iB = 1, nCntrc_a
               Do jB = 1, nCntrc_a
                  ijB = iB + (jB-1)*nCntrc_a
                  Tmp=Zero
                  Do iC = 1, iCmp_a
                     ijC = iC + (iC-1)*iCmp_a
                     iFrom = ijB + (ijC-1)*nCntrc_a**2
                     Tmp = Tmp + Tmp3(iFrom)
                  End Do
                  Shells(iShll_a)%FockOp(iB,jB) = Tmp/DBLE(iCmp_a)
               End Do
            End Do
            If (Allocated(FockOp_t)) Call mma_deallocate(FockOp_t)
#ifdef _DEBUG_
            Call RecPrt('Actual Fock operator',' ',
     &                  Shells(iShll_a)%FockOp,nCntrc_a,nCntrc_a)
#endif
            Call mma_deallocate(Tmp3)
            Call mma_deallocate(Scr1)
            Call mma_deallocate(Scr2)
*                                                                      *
************************************************************************
*                                                                      *
         End Do  ! iAng
*                                                                      *
************************************************************************
*                                                                      *
*        Deallocate the memory for the reference Fock operator
*
         Do iShll_r = jShll+1, iShll
            If (Allocated(Shells(iShll_r)%Exp))
     &          Call mma_deallocate(Shells(iShll_r)%Exp)
            Shells(iShll_r)%nExp=0
            If (Allocated(Shells(iShll_r)%FockOp))
     &          Call mma_deallocate(Shells(iShll_r)%FockOp)
            Shells(iShll_r)%nFockOp=0
            If (Allocated(Shells(iShll_r)%pCff))
     &          Call mma_deallocate(Shells(iShll_r)%pCff)
            If (Allocated(Shells(iShll_r)%Cff_c))
     &          Call mma_deallocate(Shells(iShll_r)%Cff_c)
            If (Allocated(Shells(iShll_r)%Cff_p))
     &          Call mma_deallocate(Shells(iShll_r)%Cff_p)
            Shells(iShll_r)%nExp=0
            Shells(iShll_r)%nBasis=0
         End Do
*                                                                      *
************************************************************************
*                                                                      *
*
         Charge_Actual=DBLE(iAtmNr(iCnttp))
         Charge_Effective=Charge(iCnttp)
         qTest=Test_Charge -
     &         (Charge_Actual-Charge_Effective)
c         write(6,*)'qtest, Test_Charge = ',qtest, Test_Charge
c         write(6,*)'Charge_Actual,Charge_Effective = ',
c     &               Charge_Actual,Charge_Effective
         If (qTest.eq.Zero.or.Charge(iCnttp).eq.Zero) Then
            dbsc(iCnttp)%FOp=.TRUE.
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
            dbsc(iCnttp)%FOp=.FALSE.
         End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      End Do ! iCnttp
*                                                                      *
************************************************************************
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
      Do iCnttp = 1, nCnttp
         If(dbsc(iCnttp)%Aux .or.
     &      dbsc(iCnttp)%Frag .or.
     &      dbsc(iCnttp)%nFragType.gt.0 .or.
     &      dbsc(iCnttp)%FOp) Cycle
*
         Do_FckInt = Do_FckInt .and. dbsc(iCnttp)%FOp ! To be activated!
*
      End Do
      Call mma_deallocate(STDINP)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
