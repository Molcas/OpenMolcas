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
* Copyright (C) 2012, Roland Lindh                                     *
************************************************************************
      Subroutine Mk_aCD_acCD_Shells(Info,nInfo,iCnttp,W2L)
************************************************************************
*                                                                      *
*    Objective: To generate aCD auxiliary basis sets on-the-fly.       *
*                                                                      *
* Called from: Mk_RICD_Shells                                          *
*                                                                      *
*     Author: Roland Lindh, Dept. of Chemistry - Angstrom              *
*                                                                      *
************************************************************************
      Use Basis_Info
      Implicit Real*8 (A-H,O-Z)
      External Integral_RICD
#include "itmax.fh"
#include "info.fh"
#include "SysDef.fh"
#include "real.fh"
#include "print.fh"
#include "status.fh"
#include "stdalloc.fh"
      Integer, Allocatable :: iList2_c(:,:), iList2_p(:,:), iD_c(:),
     &                        Con(:), ConR(:,:), Prm(:), Indkl_p(:),
     &                        AL(:), LTP(:,:), iD_p(:), Indkl(:)
      Real*8, Allocatable :: Wg(:), Vec(:), Scr(:), TP(:), tVt(:),
     &                       Q(:), A(:), Z(:), tVp(:), tVtF(:), C(:),
     &                       Temp(:), QTmp(:), Tmp(:)
      Real*8, Allocatable :: TInt_c(:), TInt_p(:), ADiag(:)
      Real*8 :: Dummy(1)=[0.0D0]
*                                                                      *
************************************************************************
*                                                                      *
!#define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
      Real*8, Allocatable :: H(:), U(:), tVtInv(:)
#endif
      Logical Hit, Found, Diagonal, Keep_Basis, In_Core, W2L
      Character*80 BSLbl, Label
      Character*80 atom,type,author,basis,CGTO, Aux
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Interface
         Subroutine  Drv2El_Atomic_NoSym(Integral_RICD,ThrAO,
     &                                   iCnttp,jCnttp,
     &                                   TInt,nTInt,
     &                                   In_Core,ADiag,Lu_A,ijS_req,
     &                                   Keep_Shell)
         External Integral_RICD
         Real*8 ThrAO
         Integer iCnttp, jCnttp, nTInc, Lu_A, ijS_req, Keep_Shell
         Logical In_Core
         Real*8, Allocatable :: TInt(:), ADiag(:)
         End Subroutine
      End Interface
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*---- Statement Function
*
      iTri(i,j)=Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
       iPrint=99
#else
       iPrint=5
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Max_Cnt=0
      ThrAO=Zero
      mData=4
      nCnttp_Start = nCnttp
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Loop now over all unique valence basis sets and generate the
*     corresponding aCD auxiliary basis sets. Note that there are two
*     different types of aCD auxiliary basis sets, aCD and acCD.
*
      nSO_p=0
      nTheta_All=0
*                                                                      *
************************************************************************
*                                                                      *
*
*     Pick up the threshold for the CD procedure. Note that basis
*     sets might have individual accuracy!
*
      mdc = mdciCnttp(iCnttp)
      Thr_aCD=aCD_Thr(iCnttp)*Thrshld_CD
*
      nTest= nVal_Shells(iCnttp)-1
*                                                                      *
************************************************************************
*                                                                      *
      If (Skip_High_AC) Then
*
*        Pick up the angular index of the highest valence shell
*
         If (iAtmNr(iCnttp).le.2) Then
            iVal=0
         Else If (iAtmNr(iCnttp).le.10) Then
            iVal=1
         Else If (iAtmNr(iCnttp).le.18) Then
            iVal=1
         Else If (iAtmNr(iCnttp).le.36) Then
            iVal=2
         Else If (iAtmNr(iCnttp).le.54) Then
            iVal=2
         Else If (iAtmNr(iCnttp).le.86) Then
            iVal=3
         Else
            iVal=3
         End If
         Keep_All = 2*nTest
*        Find the number of polarization shells
         iZ=Max(0,nTest-iVal)
*        Reduce the product basis from excessive shells
         Keep_Shell=Keep_All-iZ
      Else
         Keep_Shell=iTabmx
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Define some parameters to facilitate the atomic calculation
*
      iShell = nVal_Shells(iCnttp)
      nShlls=iShell
*                                                                      *
************************************************************************
*                                                                      *
*     Use the name of the old valence basis
*
      Label=Bsl_Old(iCnttp)
*
      Hit=.True.
      Call Decode(Label,atom,1,Hit)
      Hit=.True.
      Call Decode(Label,type,2,Hit)
      Hit=.True.
      Call Decode(Label,author,3,Hit)
      Hit=.True.
      Call Decode(Label,basis,4,Hit)
      Hit=.True.
      Call Decode(Label,CGTO,5,Hit)
      Hit=.False.
      Call Decode(Label,Aux,6,Hit)
      If (.Not.Hit) Aux = ' '
*
      n=Index(Atom,' ')-1
      Label=' '
      Label(1:n+1)=atom(1:n)//'.'
      nn = n + 1
*
      n=Index(Type,' ')-1
      If (Do_acCD_Basis) Then
         Label(nn+1:nn+n+23)=Type(1:n)//'....acCD-aux-basis.'
      Else
         Label(nn+1:nn+n+22)=Type(1:n)//'....aCD-aux-basis.'
      End If
*
      Indx=Index(Label,' ')
      BSLbl=' '
      BSLbl(1:Indx-1)=Label(1:Indx-1)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     C O N T R A C T E D    S E C T I O N
*
*     Run in contracted mode to generate the auxiliary basis for the
*     aCD primitive product basis.
*
      Call Flip_Flop(.False.)
*                                                                      *
************************************************************************
*                                                                      *
*     Define AOtSO
*
      iAO = 0
      iSO = 0
      nSO=0
      Do iAng = 0, nTest
         iShll_ = ipVal(iCnttp) + iAng
         nCmp = (iAng+1)*(iAng+2)/2
         If (Prjct(iShll_)) nCmp = 2*iAng+1
         iSO = 0
         If (nBasis_Cntrct(iShll_).ne.0 .and.
     &       nExp(iShll_).ne.0) Then
            Do iCmp = 1, nCmp
               iAO = iAO + 1
               iAOtSO(iAO,0) = iSO + 1
               nCont = nBasis(iShll_)
               Do iCont = 1, nCont
                   iSO = iSO + 1
               End Do
            End Do
         End If
         nSO=nSO+iSO
      End Do
*
*     Generate list
*
      nPhi_All=nSO*(nSO+1)/2
      Call mma_allocate(iList2_c,mData*2,nPhi_All,label='iList2_c')
      Call Mk_List2(iList2_c,nPhi_All,mData,nSO,iCnttp,
     &              nTest,ipVal,Mxdbsc,Prjct,MxShll,nBasis,0)
*                                                                      *
************************************************************************
*                                                                      *
*     If the full product basis is used no need for decomposition!
*
      If (Thr_aCD.eq.0.0D0) Then
         nTInt_c=nPhi_All
         Call mma_allocate(iD_c,nTInt_c,label='iD_c')
         Do i = 1, nTInt_c
            iD_c(i) = i
         End Do
         NumCho_c=nTInt_c
         Go To 1881
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Generate atomic two-electron integrals to decompose.
*
      ijS_req=0
      Call Drv2El_Atomic_NoSym(Integral_RICD,ThrAO,iCnttp,iCnttp,
     &                         TInt_c,nTInt_c,
     &                         In_Core,ADiag,Lu_A,ijS_req,Keep_Shell)
*                                                                      *
************************************************************************
*                                                                      *
*     Let us now decompose and retrieve the most important
*     contracted products, indicies stored in iD_c
*
      Call mma_allocate(iD_c,nTInt_c,label='iD_c')
*
*     Temporary code for weights to be used in a MS-aCD/acCD
*     scheme. Currently set to unit giving the convential
*     all purpose aCD/acCD auxiliary basis sets.
*
      Call mma_allocate(Wg,nTInt_c,label='Wg')
      call dcopy_(nTInt_c,[1.0D0],0,Wg,1)
*
      If (In_Core) Then
#ifdef _DEBUG_
         If (iPrint.ge.99)
     &   Call RecPrt('TInt_c',' ',TInt_c,nTInt_c,nTInt_c)
#endif
         Call mma_allocate(Vec,nTInt_c**2,label='Vec')
*
         Call CD_InCore_p_w(TInt_c,nTInt_c,
     &                      Wg,Vec,nTInt_c,iD_c,NumCho_c,
     &                      Thr_aCD,iRC)
*
         If (iRC.ne.0) Then
            Call WarningMessage(2,'Error in Mk_RICD_Shells')
            Write (6,*) 'Mk_aCD_Shells: CD_InCore_p(c) failed!'
            Call Abend()
         End If
#ifdef _DEBUG_
         If (iPrint.ge.99)
     &   Call RecPrt('Vec',' ',Vec,nTInt_c,NumCho_c)
#endif
         Call mma_deallocate(TInt_c)
         Call mma_deallocate(Vec)
*
      Else    ! out-of-core part
*
         Call GetMem('Scr','Max','Real',iDum,lScr)
         lScr=Min(lScr-2*nTInt_c,nTInt_c**2+3*nTInt_c)
         Call mma_Allocate(Scr,lScr,label='Scr')
*
         iSeed=Lu_A+1
         Lu_B=IsFreeUnit(iSeed)
         Call DaName_MF_WA(Lu_B,'AVEC1')
*
         Call Get_Pivot_idx_w(ADiag,Wg,nTInt_c,
     &                        NumCho_c,Lu_A,Lu_B,iD_c,
     &                        Scr,lScr,Thr_aCD)
*
         Call mma_deallocate(Scr)
         Call mma_deallocate(ADiag)
         Call DaEras(Lu_B)
         Call DaEras(Lu_A)
*
      End If
*
      Call mma_deallocate(Wg)
*
 1881 Continue
*
      If (NumCho_c.lt.1) Then
         Call WarningMessage(2,'Error in Mk_RICD_Shells')
         Write (6,*) 'Mk_aCD_Shells: NumCho_c.lt.1 is illegal!'
         Call Abend()
      End If
*
#ifdef _DEBUG_
      If (iPrint.ge.49) Then
         Write (6,*) ' Thr_aCD:',Thr_aCD
         Write (6,*) 'NumCho_c:',NumCho_c
         Call iVcPrt('iD_c',' ',iD_c,NumCho_c)
      End If
#endif
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Define AOtSO for primitive integral calculations.
*
      If (Do_acCD_Basis) Then
         iAO = 0
         iSO = 0
         nSO_p=0
         Do iAng = 0, nTest
            iShll_ = ipVal(iCnttp) + iAng
            nCmp = (iAng+1)*(iAng+2)/2
            If (Prjct(iShll_)) nCmp = 2*iAng+1
            iSO = 0
            Do iCmp = 1, nCmp
               iAO = iAO + 1
               iAOtSO(iAO,0) = iSO + 1
               nCont = nExp(iShll_)
               Do iCont = 1, nCont
                   iSO = iSO + 1
               End Do
            End Do
            nSO_p=nSO_p+iSO
         End Do
      End If
*                                                                      *
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
*        Loop through angular products. Note that all the products
*        of an atom require multiple basis sets since Seward is not
*        structured to handle more than one shell of a specific
*        angular at the time, i.e. a basis set contains only, for
*        example, one d-shell. For an atomic basis spd we will have
*        the p*p and d*s resulting in two independent shells with
*        the same total angular momentum, d.
*
         iShll=Mx_Shll - 1
*
*        Start now looping over the products and analys the result
*        of the CD. Note the very peculiar loop structure over
*        iBS, iAng, and jAng. This to reduce the number of
*        created basis sets.
*
         nBS = (nTest+2)/2
         Do iBS = 0, nBS-1
            iAngMin=iBS
            iAngMax=nTest - iBS
*
            nCnttp=nCnttp+1
            Keep_Basis = .False.
*
            If (nCnttp.gt.Mxdbsc) Then
               Call WarningMessage(2,'Error in Mk_RICD_Shells')
               Write (6,*) 'Mk_RICD_Shells: Increase Mxdbsc'
               Call Abend()
            End If
*
*           Some generic setting of information
*
            SODK(nCnttp)=.False.
            Bsl(nCnttp)=Label
            Bsl_Old(nCnttp)=Bsl(nCnttp)
            AuxCnttp(nCnttp)=.True.
            Charge(nCnttp)=Zero
            pChrg(nCnttp)=pChrg(iCnttp)
            Fixed(nCnttp)=Fixed(iCnttp)
            Parent_iCnttp(nCnttp)=iCnttp
            nOpt(nCnttp) = 0
            ipVal(nCnttp) = iShll+1
            ipPrj(nCnttp) = -1
            ipSRO(nCnttp) = -1
            ipSOC(nCnttp) = -1
            ipPP(nCnttp)  = -1
            nPrj_Shells(nCnttp) = 0
            nSRO_Shells(nCnttp) = 0
            nSOC_Shells(nCnttp) = 0
            nPP_Shells(nCnttp)  = 0
            AuxCnttp(nCnttp) =.True.
            lAux =.True.
            ECP(nCnttp)=.False.
            aCD_Thr(nCnttp)=aCD_Thr(iCnttp)
            fmass(nCnttp)=fmass(iCnttp)
*                                                                      *
************************************************************************
*                                                                      *
*           Loop over shell pairs
*
            jShll=iShll
            Do iAng = 0, iAngMax
               jAngMax=Min(iAng,iAngMin)
               iShll_=ipVal(iCnttp)+iAng
               If (iAng.eq.iAngMax) jAngMax=iAngMax
               If (iAng.lt.iAngMin) jAngMax=0
               jAngMin=iAngMin
               If (iAng.le.iAngMin) jAngMin=0
               Do jAng = jAngMin, jAngMax
                  jShll_=ipVal(iCnttp)+jAng
*
                  iShll = iShll + 1
#ifdef _DEBUG_
                  If (iAng.eq.3.and.jAng.eq.1) Then
                     iPrint=49
C                    iPrint=99
                  Else
                     iPrint= 5
                  End If
                  If (iPrint.ge.49) Then
                     Write (6,*)
                     Write (6,*) 'iAng,jAng=',iAng,jAng
                     Write (6,*) 'iAngMax=',iAngMax
                  End If
#endif
                  If (iShll.gt.MxShll) Then
                     Call WarningMessage(2,'Error in Mk_RICD_Shells')
                     Write (6,*) 'Mk_RICD_Shells: iShll.gt.MxShll'
                     Write (6,*) 'iShll,MxShll=',iShll,MxShll
                     Call Abend()
                  End If
                  Diagonal=iAng.eq.jAng
*
*                 Examine if any contracted products of these two shells
*                 survived the CD procedure, or that it is an empty shell.
*
                  Found=.False.
                  kShll=-1
                  lShll=-1
                  Do iCho_c = 1, NumCho_c
                     ijSO = iD_c(iCho_c)
                     kAng=iList2_c(1,ijSO)
                     lAng=iList2_c(2,ijSO)
                     Found = Found .or.
     &                   iAng.eq.kAng .and. jAng.eq.lAng
                     If ( iAng.eq.kAng .and. jAng.eq.lAng ) Then
                        kShll=iList2_c(7,ijSO)
                        lShll=iList2_c(8,ijSO)
                     End If
                  End Do
*
*                 Fake not found for shells which should explicitly be
*                 empty.
*
                  Found = Found .and. jAng.ge.iAngMin
     &                          .and. iAng.ge.iAngMin
     &                          .and. iAng+jAng.le.Keep_Shell
                  Keep_Basis = Found .or. Keep_Basis
#ifdef _DEBUG_
                  If (iPrint.ge.49)
     &            Write (6,*) 'Found,kShll,lShll=',Found,kShll,lShll
#endif
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*                    P R I M I T I V E   S E C T I O N
*
*                    Run in uncontracted mode to produce a SLIM
*                    primitive product  basis.
*
                  If (Do_acCD_Basis.and.Found) Then
*
                     Call Flip_Flop(.True.)
*
*                    Generate list
*
                     npi=nExp(iShll_)
                     nCmpi=(iAng+1)*(iAng+2)/2
                     If (Prjct(iShll_)) nCmpi=2*iAng+1
                     npj=nExp(jShll_)
                     nCmpj=(jAng+1)*(jAng+2)/2
                     If (Prjct(jShll_)) nCmpj=2*jAng+1
                     If (iAng.eq.jAng) Then
                        nTheta_All=npi*nCmpi*(npi*nCmpi+1)/2
                     Else
                        nTheta_All=npi*nCmpi*npj*nCmpj
                     End If
*
                     Call mma_allocate(iList2_p,2*mData,nTheta_all,
     &                                 label='iList2_p')
*
                     ijS_Req=(iAng+1)*iAng/2 + jAng + 1
*
                     Call Mk_List2(iList2_p,nTheta_All,mData,
     &                             nSO_p,iCnttp, nTest,ipVal,Mxdbsc,
     &                             Prjct,MxShll,nBasis,ijS_Req)
*                                                                      *
************************************************************************
*                                                                      *
*                    Generate atomic two-electron integrals
*
                     Call Drv2El_Atomic_NoSym(Integral_RICD,
     &                                        ThrAO,iCnttp,iCnttp,
     &                                        TInt_p,nTInt_p,
     &                                        In_Core,ADiag,Lu_A,
     &                                        ijS_Req,Keep_Shell)
*
                     If (.NOT.In_Core) Then
                        Call WarningMessage(2,'Error in Mk_RICD_Shells')
                        Write (6,*) 'Out-of-core acCD not implemented!'
                        Call Abend()
                     End If
#ifdef _DEBUG_
                     If (iPrint.ge.99)
     &               Call RecPrt('TInt_p','(5G20.11)',
     &                           TInt_p,nTInt_p,nTInt_p)
#endif
                     Call Flip_Flop(.False.)
*
                  End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*
*                 Now mimic the procedure of GetBS!
*
                  iStrt=ipExp(iShll)
*                                                                      *
************************************************************************
*                                                                      *
*                 Working on the CONTRACTED functions.
*
*                 This section is identical for acCD and aCD auxiliary
*                 basis sets!
*                                                                      *
************************************************************************
*                                                                      *
                  If (Found) Then
*
                     lAng = iAng + jAng
*
*                    Now figure out how many and which!
*
                     nk=nBasis_Cntrct(kShll)
                     nl=nBasis_Cntrct(lShll)
                     If (Diagonal) Then
                        nCntrc_Max=nk*(nk+1)/2
                     Else
                        nCntrc_Max=nk*nl
                     End If
#ifdef _DEBUG_
                     If (iPrint.ge.49) Write (6,*) 'nCntrc_Max=',
     &                                              nCntrc_Max
#endif
                     Call mma_allocate(Con,nCntrc_Max,label='Con')
                     Call mma_allocate(ConR,2,nCntrc_Max,label='ConR')
                     Call IZero(Con,nCntrc_Max)
                     Call IZero(ConR,2*nCntrc_Max)
                     nCntrc=0
                     Do iCho_c = 1, NumCho_c
                        ijSO = iD_c(iCho_c)
                        kAng=iList2_c(1,ijSO)
                        lAng=iList2_c(2,ijSO)
                        If (kAng.eq.iAng.and.lAng.eq.jAng) Then
*
*                          Pick up the radial index!
*
                           ik=iList2_c(5,ijSO)
                           il=iList2_c(6,ijSO)
*
                           If (Diagonal) Then
                              ikl=iTri(ik,il)
                           Else
                              ikl=(il-1)*nk + ik
                           End If
*
*                          Note that this migh be done several time
*                          since several angular pairs might have the same
*                          radial function!
*
                           If (Con(ikl).eq.0) Then
                               nCntrc=nCntrc+1
                               Con(ikl)=1
                               ConR(1,nCntrc)=ik
#ifdef _DEBUG_
                              If (iPrint.ge.49) Then
                                 Write (6,*) 'iCho_c,  ijSO=',
     &                                        iCho_c+1,ijSO
                              End If
#endif
                               ConR(2,nCntrc)=il
                           End If
                        End If
                     End Do    !  iCho_c
#ifdef _DEBUG_
                     If (iPrint.ge.49) Then
                        Write (6,*) 'nCntrc=',nCntrc
                        Call iVcPrt('Con',' ',Con,nCntrc_Max)
                        Call iVcPrt('ConR',' ',ConR,2*nCntrc)
                        Write (6,*)
                        Write (6,*) 'ConR'
                        Write (6,'(30I3)')
     &                        (ConR(1,i),i=1,nCntrc)
                        Write (6,'(30I3)')
     &                        (ConR(2,i),i=1,nCntrc)
                     End If
#endif
*
                  Else
*
*                    Let us put in an empty shell!
*
                     nk=0
                     nl=0
                     nCntrc=0
                  End If
*                                                                      *
************************************************************************
*                                                                      *
*                 Work on the PRIMITIVE products!
*
*                 Here the work is trivial in case of the aCD basis
*                                                                      *
************************************************************************
*                                                                      *
                  If (Found) Then
*                                                                      *
************************************************************************
*                                                                      *
*                    Produce the SLIM primitive products
*                                                                      *
************************************************************************
*                                                                      *
                     If (Do_acCD_Basis) Then
*
*                       Now figure out how many and which!
*
                        npk=nExp(kShll)
                        npl=nExp(lShll)
                        If (Diagonal) Then
                           nPrim_Max=npk*(npk+1)/2
                        Else
                           nPrim_Max=npk*npl
                        End If
#ifdef _DEBUG_
                        If (iPrint.ge.49)
     &                     Write (6,*) 'nPrim_Max:',nPrim_Max
#endif
                        Call mma_allocate(Prm,nPrim_Max,label='Prm')
                        Call IZero(Prm,nPrim_Max)
*
*                       Pick up the diagonal elements from TInt_p
*                       corresponding to this shell pair. We sum over
*                       the angular parts identical to those of the
*                       contracted.
*
                        nCompA=(iAng+1)*(iAng+2)/2
                        nCompB=(jAng+1)*(jAng+2)/2
                        nAB=nCompA*nCompB
                        Call mma_allocate(AL,nAB,label='AL')
*
*                       First make a list from the contracted which
*                       angular products to include.
*
                        Call Mk_AngList(AL,nCompA,nCompB,
     &                                  iD_c,NumCho_c,
     &                                  iList2_c,nPhi_All,
     &                                  2*mData,iAng,jAng)
*
                        Call mma_allocate(TP,nPrim_Max**2,label='TP')
                        Call mma_allocate(LTP,2,nPrim_Max,label='LTP')
                        Call Mk_TInt_P(TInt_p,nTheta_All,
     &                                 TP,nPrim_Max,
     &                                 AL,nCompA,nCompB,
     &                                 iList2_p,nTheta_All,
     &                                 2*mData,iAng,jAng,npk,npl,LTP)
*
#ifdef _DEBUG_
                        If (iPrint.ge.99) Then
                        Call RecPrt('TIntP','(5G20.10)',
     &                              TP,nPrim_Max,nPrim_Max)
                        Call iVcPrt('List_TP',' ',LTP,2*nPrim_Max)
                        End If
#endif
*                       Let us now decompose and retrieve the most
*                       important primitive products, indicies stored in
*                       iD_p
*
                        Call mma_allocate(iD_p,nPrim_Max,label='iD_p')
                        Call mma_allocate(Vec,nPrim_Max**2,label='Vec')
*
                        Thrshld_CD_p = Thr_aCD*2.0D-1
 3377                   Continue
                        Call CD_InCore_p(TP,nPrim_Max,Vec,nPrim_Max,
     &                                   iD_p,NumCho_p,Thrshld_CD_p,iRC)
                        If (NumCho_p.lt.1) Then
                           Call WarningMessage(2,
     &                           'Error in Mk_RICD_Shells')
                           Write (6,*) 'Mk_aCD_Shells: '
     &                               //'NumCho_p.lt.1 is illegal!'
                           Write (6,*) 'iAng,jAng=',iAng,jAng
                           Write (6,*) 'nPrim_Max=',nPrim_Max
                           Write (6,*) 'NumCho_p=',NumCho_p
                           Write (6,*) 'iRC=',iRC
                           Call Abend()
                        End If
*
#ifdef _DEBUG_
                        If (iPrint.ge.49) Then
                        Write (6,*) 'Thrshld_CD_p:',Thrshld_CD_p
                        Write (6,*) 'NumCho_p    :',NumCho_p
                        Call iVcPrt('iD_p',' ',iD_p,NumCho_p)
                        End If
                        If (iPrint.ge.99)
     &                  Call RecPrt('Vec',' ',Vec,nPrim_Max,NumCho_p)
#endif
                        If (NumCho_p.lt.nCntrc) Then
                           If (iPrint.ge.6) Then
                           Write (6,*) 'W a r n i n g!'
                           Write (6,*) 'Fewer primitive functions than'
     &                               //' contracted functions!'
                           Write (6,*) 'NumCho_p=',NumCho_p
                           Write (6,*) '  nCntrc=',nCntrc
                           End If
                           Thrshld_CD_p=Thrshld_CD_p*0.5
                           If (Thrshld_CD_p.le.1.0D-14) Then
                              Call WarningMessage(2,
     &                           'Error in Mk_RICD_Shells')
                              Write (6,*) 'Thrshld_CD_p is too low!'
                              Write (6,*) 'iAng, jAng:',iAng,jAng
                              Call Abend()
                           End If
                           Call Mk_TInt_P(TInt_p,nTheta_All,
     &                                    TP,nPrim_Max,
     &                                    AL,nCompA,nCompB,
     &                                    iList2_p,nTheta_All,
     &                                    2*mData,iAng,jAng,npk,npl,LTP)
                           Go To 3377
                        End If
                        Call mma_deallocate(TP)
                        Call mma_deallocate(Vec)
*
                        Do iCho_p = 1, NumCho_p
                           ijSO = iD_p(iCho_p)
                           Prm(ijSO)=1
                        End Do
                        nPrim=NumCho_p
                        Call mma_allocate(Shells(iShll)%Exp,nPrim,
     &                                    Label='ExpacCD')
                        Shells(iShll)%nExp=nPrim
*
                        iEnd = iStrt - 1
*
#ifdef _DEBUG_
                        If (iPrint.ge.49) Then
                        Write (6,*) 'nPrim=',nPrim
                        Call iVcPrt('Prm',' ',Prm,nPrim_Max)
                        End If
#endif
                        Call mma_allocate(Indkl_p,nPrim_Max,
     &                                    label='Indkl_p')
                        Call Mk_Indkl(Prm,Indkl_p,nPrim_Max)
*
*                       Observe that the exponents are ordered according
*                       to their importance as given by the CD.
*
                        Do iCho_p = 1, NumCho_p
                           iTheta = iD_p(iCho_p)
                           ik=LTP(1,iTheta)
                           il=LTP(2,iTheta)
                           Exp_i=Shells(kShll)%Exp(ik)
                           Exp_j=Shells(lShll)%Exp(il)
                           Shells(iShll)%Exp(iCho_p)=Exp_i+Exp_j
                        End Do
#ifdef _DEBUG_
                        If (iPrint.ge.49)
     &                    Call RecPrt('SLIM Exponents',' ',
     &                                Shells(iShll)%Exp,1,nPrim)
#endif
*                                                                      *
************************************************************************
*                                                                      *
                     Else ! Do_aCD_Basis
*                                                                      *
************************************************************************
*                                                                      *
*                       Put in the aCD set of exponents, i.e. all unique
*                       sums.
*
                        If (Diagonal) Then
                           nPrim=nExp(kShll)*(nExp(kShll)+1)/2
                        Else
                           nPrim=nExp(kShll)*nExp(lShll)
                        End If
                        Call mma_allocate(Shells(iShll)%Exp,nPrim,
     &                                    Label='ExpaCD')
                        Shells(iShll)%nExp=nPrim
                        iEnd = iStrt - 1
*
                        iOff = 0
                        Do ip_Exp = 1, nExp(kShll)
                           jp_Exp_Max = nExp(lShll)
                           If (Diagonal) jp_Exp_Max = ip_Exp
                           Do jp_Exp = 1, jp_Exp_Max
                              iOff = iOff + 1
                              Shells(iShll)%Exp(iOff)=
     &                                   Shells(kShll)%Exp(ip_Exp)
     &                                  +Shells(lShll)%Exp(jp_Exp)
                           End Do
                        End Do
*
                        If (iOff.ne.nPrim) Then
                           Call WarningMessage(2,
     &                           'Error in Mk_RICD_Shells')
                           Write (6,*) 'Mk_aCD_Shell: iOff.ne.iEnd'
                           Call Abend()
                        End If
*
#ifdef _DEBUG_
                        If (iPrint.ge.49) Then
                        If (Diagonal) Then
                           Call TriPrt('aCD Exponents',' ',
     &                                 Shells(iShll)%Exp,
     &                                 nExp(kShll))
                        Else
                           Call RecPrt('aCD Exponents',' ',
     &                                 Shells(iShll)%Exp,
     &                                             nExp(kShll),
     &                                             nExp(lShll))
                        End If
                        End If
#endif
                     End If
*
                  Else
*                                                                      *
************************************************************************
*                                                                      *
*                    An empty shell
*
                     nPrim=0
                     Shells(iShll)%nExp=nPrim
                     iEnd = iStrt - 1
*                                                                      *
************************************************************************
*                                                                      *
                  End If ! Found
*                                                                      *
************************************************************************
*                                                                      *
                  lAng=iAng+jAng
                  iAngMx=Max(iAngMx,lAng)
                  MaxPrm(lAng)=Max(MaxPrm(lAng),nPrim)
*
#ifdef _DEBUG_
                  If (iPrint.ge.49) Then
                  Write (6,*)
                  Write (6,*) 'iShll=',iShll
                  Write (6,*) 'nPrim,nCntrc=',nPrim,nCntrc
                  Write (6,*) 'lAng=',lAng
                  Write (6,*) 'MaxPrm(lAng)=',MaxPrm(lAng)
                  End If
#endif
*
                  nExp(iShll)=nPrim
                  nBasis_Cntrct(iShll)=nCntrc
                  iStrt = iEnd + 1
*                                                                      *
************************************************************************
*                                                                      *
                  Call mma_allocate(Shells(iShll)%Cff_c,
     &                              nPrim,nCntrc,2,Label='Cff_c')
                  Call mma_allocate(Shells(iShll)%pCff,
     &                              nPrim,nCntrc,Label='pCff')
                  Shells(iShll)%nBasis = nCntrc
                  Call mma_allocate(Shells(iShll)%Cff_p,
     &                              nPrim,nPrim ,2,Label='Cff_p')
*
                  iEnd  = iStrt -1
                  iEndc = iEnd
*                                                                      *
************************************************************************
*                                                                      *
*                 C O N T R A C T I O N    C O E F F I C I E N T S
*
                  If (Found) Then
*                                                                      *
************************************************************************
*                                                                      *
*                    For SLIM basis sets
*                                                                      *
************************************************************************
*                                                                      *
                     If (Do_acCD_Basis) Then
*
*                       Alright this is a bit more elaborate than for
*                       the aCD basis set. Surprise!
*
*                       Some care has to be taken here. There might be
*                       different angular products, for example, px*px
*                       and px*py, which carry the same radial part but
*                       have different angular part! To overcome this
*                       possible source of redundancy we use the sum of
*                       such terms in the fitting procedure!
*
                        nTheta=nPrim
                        If (iAng.eq.jAng) Then
                           nTheta_Full=nExp(kShll)*(nExp(kShll)+1)/2
                        Else
                           nTheta_Full=nExp(kShll)*nExp(lShll)
                        End If
                        nPhi=nCntrc
*
*                       Generate the (theta'|V|theta') matrix in the
*                       SLIM primitive product basis.
*
                        Call mma_allocate(tVt,nTheta**2,label='tVt')
                        Call Mk_tVt(TInt_p,nTInt_p,
     &                              tVt,nTheta,iList2_p,2*mData,
     &                              Prm,nPrim_Max,
     &                              iAng,jAng,nExp(kShll),nExp(lShll),
     &                              Indkl_p,nPrim_Max,
     &                              AL,nCompA,nCompB)
*
#ifdef _DEBUG_
                        If (iPrint.ge.49) Then
                           Write (6,*)
                           Write (6,*) 'tVt(Diag)'
                           Write (6,*)
     &                     (tVt(i),i=1,nTheta**2,nTheta+1)
                        End If
                        If (iPrint.ge.99) Then
                        Call RecPrt('tVt',' ',tVt,nTheta,nTheta)
                        Call iVcPrt('iD_p',' ',iD_p,NumCho_p)
                        End If
#endif
*
*                       Generate (theta'|V|theta')^{-1}
*
*                       Let's do a Cholesky decomposition with pivoting
*                       according to the previous CD.
*
                        nTri=nTheta*(nTheta+1)/2
                        Call mma_allocate(Q,nTri,label='Q')
                        Call mma_allocate(A,nTri,label='A')
                        Call mma_allocate(Z,nTheta,label='Z')
                        Do iCho_p = 1, NumCho_p
                           iTheta_full=iD_p(iCho_p)
                           iTheta     =Indkl_p(iTheta_full)
                           Do jCho_p = 1, iCho_p
                              jTheta_full=iD_p(jCho_p)
                              jTheta     =Indkl_p(jTheta_full)
                              ijT= iCho_p*(iCho_p-1)/2 + jCho_p
                              ijS= (jTheta-1)*nTheta + iTheta
                              A(ijT)=tVt(ijS)
                           End Do
                        End Do
#ifdef _DEBUG_
                        If (iPrint.ge.99) Then
                        Call TriPrt('A',' ',A,nTheta)
*
                        Call mma_allocate(H,nTri,label='H')
                        Call mma_allocate(U,nTri,label='U')
                        call dcopy_(nTri,A,1,H,1)
                        Call FZero(U,nTheta**2)
                        call dcopy_(nTheta,[One],0,U,nTheta+1)
                        Call Jacob(H,U,nTheta,nTheta)
                        Call TriPrt('H','(10G20.10)',H,nTheta)
                        Call RecPrt('U',' ',U,nTheta,nTheta)
                        Call mma_deallocate(H)
                        Call mma_deallocate(U)
                        End If
#endif
*
#ifdef _DEBUG_
                        Call mma_allocate(tVtInv,nTheta*2,
     &                                    label='tVtInv')
                        If (iPrint.ge.49) Then
                        iSing=0
                        Det=0.0D0
                        Call MInv(tVt,tVtInv,iSing,Det,nTheta)
                        Write (6,*) 'iSing,Det=',iSing,Det
                        End If
#endif
                        Call mma_deallocate(tVt)
*
                        Do iTheta = 1, nTheta
                           iOff_Ak=(iTheta-1)*iTheta/2 + 1
                           iOff_Qk=(iTheta-1)*iTheta/2 + 1
                           Thrs= Thrshld_CD_p
C                          Thrs= 1.0D-12
                           Call Inv_Cho_Factor(A(iOff_Ak),iTheta,A,Q,
     &                                         iTheta,iDum,iDum,
     &                                         Dummy,0,Z,
     &                                         Dummy,Thrs,
     &                                         Q(iOff_Qk),LinDep)
                           If (LinDep.eq.1) Then
                              Call WarningMessage(2,
     &                           'Error in Mk_RICD_Shells')
                              Write (6,*) 'Mk_aCD_Shells: linear '
     &                                  //'dependence'
     &                                  //' found in tVt!'
                              Write (6,*) 'Found for vector:',iTheta
                              Call Abend()
                           End If
                        End Do
                        Call mma_deallocate(Z)
                        Call mma_deallocate(A)
#ifdef _DEBUG_
                        If (iPrint.ge.49)
     &                  Call TriPrt('Q','(9G10.3)',Q,nTheta)
#endif
*
*                       Generate the (theta'|V|theta) matrix in the
*                       mixed product basis.
*
                        Call mma_allocate(tVp,nTheta*nPhi,label='tVp')
                        Call mma_allocate(tVtF,nTheta*nTheta_Full,
     &                                    label='tVtF')
                        Call Mk_tVtF(TInt_p,nTInt_p,
     &                               tVtF,nTheta,
     &                               iList2_p,2*mData,
     &                               Prm,nPrim_Max,
     &                               iAng,jAng,nExp(kShll),nExp(lShll),
     &                               Indkl_p,nPrim_Max,
     &                               nTheta_Full,
     &                               AL,nCompA,nCompB)
                        Call mma_deallocate(AL)
#ifdef _DEBUG_
                        If (iPrint.ge.49)
     &                  Call RecPrt('tVtF',' ',tVtF,nTheta,nTheta_Full)
#endif
*
*                       Pick up the contraction coefficients of the aCD
*                       basis set. Be careful what this means in the
*                       case that the shells are identical!
*
                        Call mma_allocate(Indkl,nCntrc_Max,
     &                                    label='Indkl')
                        Call Mk_Indkl(Con,Indkl,nCntrc_Max)
                        Call mma_allocate(C,nTheta_Full*nPhi,label='C')
                        Call Mk_Coeffs(Shells(kShll)%Cff_c(1,1,1),
     &                                 nExp(kShll),nBasis_Cntrct(kShll),
     &                                 Shells(lShll)%Cff_c(1,1,1),
     &                                 nExp(lShll),nBasis_Cntrct(lShll),
     &                                 C,nTheta_Full,nPhi,
     &                                 iD_c,NumCho_c,
     &                                 iList2_c,2*mData,
     &                                 nPhi_All,
     &                                 Indkl,nCntrc_Max,
     &                                 nBasis_Cntrct(kShll),
     &                                 nBasis_Cntrct(lShll),
     &                                 iAng,jAng,
     &                                 Shells(kShll)%Cff_p(1,1,1),
     &                                 Shells(lShll)%Cff_p(1,1,1))
                        Call mma_deallocate(Indkl)
#ifdef _DEBUG_
                        If (iPrint.ge.49)
     &                  Call RecPrt('C',' ',C,nTheta_Full,nPhi)
#endif
*
*                       Generate the (theta'|V|phi') matrix.
*
                        Call DGEMM_('N','N',
     &                              nTheta,nPhi,nTheta_Full,
     &                              1.0d0,tVtF,nTheta,
     &                              C,nTheta_Full,
     &                              0.0d0,tVp,nTheta)
                        Call mma_deallocate(tVtF)
#ifdef _DEBUG_
                        If (iPrint.ge.49)
     &                  Call RecPrt('tVp',' ',tVp,nTheta,nPhi)
#endif
                        Call mma_deallocate(C)
*
*                       Generate the contraction coefficients of the
*                       SLIM contracted product basis in terms of the
*                       SLIM primitive product basis as
*                       Sum(nu') (mu'|V|nu')^-1  (nu'|V|i')=C_{mu',i'}
*
*                       To simplify life I will put the Q matrix into
*                       square storage.
*
                        Call mma_Allocate(Temp,nTheta**2,label='Temp')
                        Call FZero(Temp,nTheta**2)
                        Do iTheta = 1, nTheta
                           Do jTheta = 1, iTheta
                              ijT = iTheta*(iTheta-1)/2 + jTheta
                              ijS = (iTheta-1)*nTheta + jTheta
                              Temp(ijS)=Q(ijT)
                           End Do
                        End Do
                        call mma_deallocate(Q)
#ifdef _DEBUG_
                        If (iPrint.ge.49)
     &                  Call RecPrt('Q',' ',Temp,nTheta,nTheta)
#endif
*
*                       Resort the external index back to original
*                       order. The column index is external.
*
                        Call mma_allocate(QTmp,nTheta**2,label='QTmp')
                        Do iCho_p = 1, NumCho_p
                           iTheta_Full = iD_p(iCho_p)
                           iTheta      = Indkl_p(iTheta_Full)
                           call dcopy_(nTheta,Temp(iCho_p),nTheta,
     &                                       QTmp(iTheta),nTheta)
                        End Do
                        Call mma_deallocate(Temp)
#ifdef _DEBUG_
                        If (iPrint.ge.99)
     &                  Call RecPrt('Q',' ',QTmp,nTheta,nTheta)
#endif
*                       Q(T) tVp
                        Call mma_allocate(Scr,nTheta*nPhi,label='Scr')
                        Call DGEMM_('T','N',
     &                              nTheta,nPhi,nTheta,
     &                              1.0d0,QTmp,nTheta,
     &                              tVp,nTheta,
     &                              0.0d0,Scr,nTheta)
*                       QQ(T) tVp
                        Call DGEMM_('N','N',
     &                              nTheta,nPhi,nTheta,
     &                              1.0d0,QTmp,nTheta,
     &                              Scr,nTheta,
     &                              0.0d0,
     &                              Shells(iShll)%Cff_c(1,1,1),nTheta)
#ifdef _DEBUG_
                        If (iPrint.ge.99)
     &                  Call RecPrt('SLIM coeffcients',' ',
     &                               Shells(iShll)%Cff_c(1,1,1),
     &                               nTheta,nPhi)
                        If (iPrint.ge.49) Then
                        Call DGEMM_('N','N',
     &                              nTheta,nPhi,nTheta,
     &                              1.0d0,tVtInv,nTheta,
     &                              tVp,nTheta,
     &                              0.0d0,Scr,nTheta)
                        Call RecPrt('SLIM coeffcients2',' ',Scr,
     &                               nTheta,nPhi)
                        End If
                        Call mma_deallocate(tVtInv)
#endif
                        Call mma_deallocate(tVp)
*
*                       Now reorder the coefficients to the CD order of
*                       the exponents.
*
                        Call mma_allocate(Tmp,nTheta*nPhi,label='Tmp')
                        call dcopy_(nTheta*nPhi,
     &                              Shells(iShll)%Cff_c(1,1,1),1,Tmp,1)
                        Do iCho_p = 1, NumCho_p
                           iTheta_full = iD_p(iCho_p)
                           iTheta      = Indkl_p(iTheta_full)
                           iTo   =  iStrt-1 + iCho_p
                           call dcopy_(nPhi,
     &                                 Tmp(iTheta),nTheta,
     &                                 Shells(iShll)%Cff_c(iCho_p,1,1),
     &                                             nTheta)
                        End Do
                        Call mma_deallocate(Tmp)
*
*                       Modify from coefficients for normalized
*                       Gaussians to unnormalized Gaussians.
*
                        Do iCho_p = 1, NumCho_p
                           iTheta = iD_p(iCho_p)
                           ik=LTP(1,iTheta)
                           il=LTP(2,iTheta)
                           Fact = Shells(kShll)%Cff_p(ik,ik,1)
     &                          * Shells(lShll)%Cff_p(il,il,1)
                           Call DScal_(nPhi,Fact,
     &                            Shells(iShll)%Cff_c(iCho_p,1,1),
     &                                 nTheta)
                        End Do
                        Call mma_deallocate(LTP)
*
                        Call mma_deallocate(iD_p)
                        Call mma_deallocate(Indkl_p)
                        Call mma_deallocate(Scr)
                        Call mma_deallocate(QTmp)
#ifdef _DEBUG_
*                       If (iPrint.ge.49)
                        Call RecPrt('SLIM coeffcients',' ',
     &                              Shells(iShll)%Cff_c(1,1,1),
     &                              nTheta,nPhi)
#endif
*                                                                      *
************************************************************************
*                                                                      *
                     Else ! Do_aCD_Basis
*                                                                      *
************************************************************************
*                                                                      *
*                       Put in the selected set of coeffients. Note
*                       again that the order should be that according
*                       to the CD in order to prepivot, since the CD
*                       itself is implemented without pivoting.
*
                        Do iCntrc = 1, nCntrc
                           kC = ConR(1,iCntrc)
                           lC = ConR(2,iCntrc)
#ifdef _DEBUG_
                           If (iPrint.ge.49)
     &                     Write (6,*) 'kC,lC=',kC,lC
#endif
*                                                                      *
************************************************************************
*                                                                      *
*                          Form the unnormalized coefficients!
*
                           jkl = 0
                           If (Diagonal) Then
                              Do iExp_k = 1, nExp(kShll)
*
                               Coeff_kk=Shells(kShll)%Cff_c(iExp_k,kC,1)
                               Coeff_lk=Shells(lShll)%Cff_c(iExp_k,lC,1)

                                 Do iExp_l = 1 , iExp_k

                               Coeff_ll=Shells(lShll)%Cff_c(iExp_l,lC,1)
                               Coeff_kl=Shells(kShll)%Cff_c(iExp_l,kC,1)
                                    Coeff_  =Coeff_ll*Coeff_kk
     &                                      +Coeff_kl*Coeff_lk
                                    If (iExp_k.eq.iExp_l) Then
                                       Coeff_  =Coeff_  *Half
                                    End If
                                    jkl = jkl + 1

                               Shells(iShll)%Cff_c(jkl,iCntrc,1)=Coeff_

                                 End Do
                              End Do
                           Else
                              Do iExp_k = 1, nExp(kShll)

                               Coeff_k =Shells(kShll)%Cff_c(iExp_k,kC,1)

                                 Do iExp_l = 1 , nExp(lShll)

                               Coeff_l =Shells(lShll)%Cff_c(iExp_l,lC,1)

                                    Coeff_kl=Coeff_l*Coeff_k

                                    jkl = jkl + 1
                              Shells(iShll)%Cff_c(jkl,iCntrc,1)=Coeff_kl
                                 End Do
                              End Do
                           End If
*
                        End Do ! iCntrc
#ifdef _DEBUG_
                        If (iPrint.ge.49)
     &                  Call RecPrt('aCD Coefficients','(6G20.12)',
     &                              Shells(iShll)%Cff_c(1,1,1),
     &                              nPrim,nCntrc)
#endif
*                                                                      *
************************************************************************
*                                                                      *
                     End If
*                                                                      *
************************************************************************
*                                                                      *
                     Shells(iShll)%Cff_c(:,:,2)=
     &                       Shells(iShll)%Cff_c(:,:,1)
*
                     Call mma_deallocate(Con)
                     Call mma_deallocate(ConR)
                     If (Do_acCD_Basis) Call mma_deallocate(Prm)
*
*                    Put in unit matrix of uncontracted set
*
                     Shells(iShll)%Cff_p(:,:,1)=Zero
                     ForAll (i=1:nPrim) Shells(iShll)%Cff_p(i,i,1)=One
*
                     Shells(iShll)%Cff_p(:,:,2)=
     &                       Shells(iShll)%Cff_p(:,:,1)
                     Call Nrmlz(Shells(iShll)%Exp,nPrim,
     &                          Shells(iShll)%Cff_p(1,1,1),
     &                          nPrim ,lAng)
#ifdef _DEBUG_
                     If (iPrint.ge.99) Then
                        Call RecPrt('uncon1',' ',
     &                               Shells(iShll)%Cff_p(:,:,1),
     &                               nPrim,nPrim)
                        Call RecPrt('uncon2',' ',
     &                               Shells(iShll)%Cff_p(:,:,2),
     &                               nPrim,nPrim)
                     End If
#endif
*
*                    OK let's do the correction now!
*
#ifdef _DEBUG_
                     If (iPrint.ge.99) Then
                     Call RecPrt('Coefficients 10',' ',
     &                           Shells(iShll)%Cff_c(:,:,1),
     &                           nPrim,nCntrc)
                     iOff = nPrim*nCntrc
                     Call RecPrt('Coefficients 20',' ',
     &                           Shells(iShll)%Cff_c(:,:,2),
     &                           nPrim,nCntrc)
                     End If
#endif
                     iOff = nPrim*nCntrc
                     Call Fix_Coeff(nPrim,nCntrc,
     &                              Shells(iShll)%Cff_c(:,:,2),
     &                              Shells(iShll)%Cff_p(:,:,1),'F')
#ifdef _DEBUG_
                     If (iPrint.ge.99) Then
                     Call RecPrt('Coefficients 1',' ',
     &                            Shells(iShll)%Cff_c(:,:,1),
     &                            nPrim,nCntrc)
                     iOff = nPrim*nCntrc
                     Call RecPrt('Coefficients 2','(6G20.13)',
     &                            Shells(iShll)%Cff_c(:,:,2),
     &                            nPrim,nCntrc)
                     End If
#endif
*
*                    Now remove any primitives with all zero
*                    coefficents!
*
                     Call Fix_Exponents(nPrim,mPrim,nCntrc,
     &                                  Shells(iShll)%Exp,
     &                                  Shells(iShll)%Cff_c(:,:,1),
     &                                  Shells(iShll)%Cff_p(:,:,1))
                     nPrim=mPrim
                     Shells(iShll)%nExp=nPrim
                     nExp(iShll)=nPrim
#ifdef _DEBUG_
                     If (iPrint.ge.99) Then
                     Call RecPrt('Coefficients 1',' ',
     &                           Shells(iShll)%Cff_c(:,:,1),
     &                           nPrim,nCntrc)
                     iOff = nPrim*nCntrc
                     Call RecPrt('Coefficients 2',' ',
     &                           Shells(iShll)%Cff_c(:,:,2),
     &                           nPrim,nCntrc)
                     End If
#endif
                  End If
*
                  If (iShll.lt.MxShll) ipExp(iShll+1) = ipExp(iShll)
*
                  nBasis(iShll)=nBasis_Cntrct(iShll)
                  If (jAng.eq.0) Then
                     Transf(iShll)=Transf(kShll)
                     Prjct(iShll)=Prjct(kShll)
                  Else
                     Transf(iShll)=.True.
                     Prjct(iShll)=.False.
                  End If
                  AuxShell(iShll)=.True.
*
                  If (Do_acCD_Basis.and.Found) Then
                     Call mma_deallocate(iList2_p)
                     Call mma_deallocate(TInt_p)
                  End If
                  Shells(iShll)%pCff(:,:) = Shells(iShll)%Cff_c(:,:,1)
*
               End Do ! jAng
            End Do ! iAng
*
            nTot_Shells(nCnttp) = iShll-jShll
            nVal_Shells(nCnttp) = iShll-jShll
*                                                                      *
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
            If (Keep_Basis) Then
               If (Show.and.nPrint(2).ge.6) Then
                  Write (6,*)
                  Write (6,*)
                  Write(6,'(1X,A,I5,A,A)')
     &                  'Basis Set ',nCnttp,' Label: ', BSLbl(1:Indx-1)
                  Write(6,'(1X,A)') 'On-the-fly basis set generation'
               End if
*
*
*              Transfer the coordinate information
*
               nCnt = dbsc(iCnttp)%nCntr
               dbsc(nCnttp)%nCntr=nCnt
               mdciCnttp(nCnttp)=mdc
!              Call allocate(dbsc(nCnttp)%Coor(1:3,nCnt))
               Call mma_allocate(dbsc(nCnttp)%Coor,3,nCnt,
     &                           Label='dbsc:C')
               dbsc(nCnttp)%Coor(:,:) = dbsc(iCnttp)%Coor(:,:)
*
*              Compute the number of elements stored in the dynamic
*              memory so far.
               nInfo = ipExp(iShll+1) - Info
               Mx_Shll=iShll+1
               Max_Shells=Mx_Shll
               Mx_mdc=mdc
*
            Else
*
*              If all the shells are empty, skip the whole basis set!
*
               nCnttp=nCnttp-1
            End If
*
         End Do   ! iBS
*
*        Done for this valence basis set.
*
         Mx_Shll = iShll + 1
         Max_Shells=Mx_Shll
*                                                                      *
************************************************************************
*                                                                      *
*        Deallocate
*
         Call mma_deallocate(iD_c)
         Call mma_deallocate(iList2_c)
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
*        Let us now Gram-Schmidt orthonormalize the auxiliary basis for
*        better numerics and balance.
*
         Do jCnttp = nCnttp_start + 1, nCnttp
            Call Renorm2(jCnttp)
         End Do
*                                                                      *
************************************************************************
*                                                                      *
*        Optionally add auxiliary basis set to the end of the
*        temporary auxiliary basis set library.
*
      If (W2L) Then
         Lu_lib=17
         Lu_lib=IsFreeUnit(Lu_lib)
         call molcas_open(Lu_lib,'RICDLIB')
         ReWind(Lu_lib)
         Do ! For ever
            Read (Lu_lib,*,END=777)
         End Do
 777     Continue
         BACKSPACE(Lu_lib)
*
         Do jCnttp = nCnttp_start + 1, nCnttp
            If (jCnttp.eq.nCnttp_start+1) Then
               Write (Lu_lib,'(A)') '/'//Label
            Else
               Write (Lu_lib,'(A)') Label
            End If
            If (jCnttp.eq.nCnttp_start+1) Then
               Write (Lu_lib,'(F6.2,2I10)') Charge(jCnttp),
     &               nVal_Shells(jCnttp)-1,nCnttp-nCnttp_start
            Else
               Write (Lu_lib,'(F6.2, I10)') Charge(jCnttp),
     &               nVal_Shells(jCnttp)-1
            End If
            Write (Lu_lib,*) ' Dummy reference line.'
            Write (Lu_lib,*) ' Dummy reference line.'
            Do iAng = 0, nVal_Shells(jCnttp)-1
               iShll_ = ipVal(jCnttp) + iAng
               iSph=0
               If (Prjct(iShll_))  iSph=1
               If (Transf(iShll_)) iSph=iSph+2
               Write (Lu_lib,'(3I10)') nExp(iShll_), nBasis(iShll_),
     &                                 iSph
*
*              Skip if the shell is empty.
*
               If (nExp(iShll_).eq.0) Cycle
*
*              Write out the exponents
*
               Write (Lu_lib,'( 5(1X,D20.13))')
     &               (Shells(iShll_)%Exp(i),i=1,nExp(iShll_))
*
*              Write out the contraction coefficients
*
               Do i = 1, nExp(iShll_)
                  Write (Lu_lib,'( 5(1X,D20.13))')
     &                  (Shells(iShll_)%Cff_c(i,j,1),
     &                        j=1,nBasis(iShll_))
               End Do
*
            End Do
         End Do
         Close(Lu_lib)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
*
      End
