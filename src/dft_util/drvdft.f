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
      Subroutine DrvDFT(h1,TwoHam,D,RepNuc,nh1,First,Dff,
     &                  lRF,KSDFT,ExFac,Do_Grad,Grad,nGrad,iSpin,
     &                  D1I,D1A,nD1,DFTFOCK)
      Implicit Real*8 (a-h,o-z)
      External LSDA, Overlap, BLYP, BPBE, B3LYP, HFS, HFB,
     &         XAlpha, LSDA5, B3LYP5, B2PLYP, TLYP, NLYP,
     &         NucAtt, NEWF, NEWF1, OLYP, O3LYP, OPBE,
     &         PBE, PBE0, PBEsol, M06L, M06, M062X, HFO,
     &         M06HF, Checker, SSB, HFG, GLYP, GPBE,
     &         HFB86, B86LYP, B86PBE, BWIG, KT3,
     &         O2PLYP,  KT2,  RGE2, REVPBE,
     &         PTCA,S12G, S12H
#include "real.fh"
#include "WrkSpc.fh"
#include "nq_info.fh"
#include "debug.fh"
#include "pamint.fh"
#include "ksdft.fh"
      Real*8 h1(nh1), TwoHam(nh1), D(nh1,2), Grad(nGrad), Vxc_ref(2)
      Real*8 D1I(nD1),D1A(nD1)
      Logical First, Dff, lRF,  Do_Grad
      Logical Do_MO,Do_TwoEl, Found
      Character*(*) KSDFT
      Character*4 DFTFOCK
*
      KSDFA = KSDFT
      lKSDFT=LEN(KSDFT)
      Debug=.False.
*                                                                      *
************************************************************************
*                                                                      *
c     Call SetQue('Trace=on')
      Call QEnter('DrvDFT')
*                                                                      *
************************************************************************
*                                                                      *
      Call Put_iScalar('Multiplicity',iSpin)
      Call Get_iScalar('nSym',mIrrep)
      Call Get_iArray('nBas',mBas,mIrrep)
*
      Call Set_Basis_Mode('Valence')
      Call Setup_iSD()
*                                                                      *
************************************************************************
*                                                                      *
      If (Do_Grad) Call FZero(Grad,nGrad)
*                                                                      *
************************************************************************
*                                                                      *
      If (iSpin.eq.1) Then
         nD=1
      Else
         nD=2
      End If
*
*     What is this?
*
      If (DFTFOCK.eq.'DIFF') nD=2
      If (DFTFOCK.eq.'ROKS') nD=2
      Call GetMem('D-DS','Allo','Real',ip_D_DS,nh1*nD)
*
*---- Get the total density
*
      Call Get_D1ao(ipD1ao,nDens)
      If (nDens.ne.nh1) Then
         Call WarningMessage(2,'DrvDFT: nDens.ne.nh1')
         Write (6,*) 'nDens=',nDens
         Write (6,*) 'nh1  =',nh1
         Call Abend()
      End If
      call dcopy_(nh1,Work(ipD1ao),1,Work(ip_D_DS),1)
*      Call RecPrt('D1ao',' ',Work(ipD1ao),nh1,1)
*
      Call GetMem('DrvXV','Free','Real',ipD1ao,nDens)
*
*---- Get the spin density
*
      If (nD.ne.1) Then
         Call Get_D1Sao(ipD1Sao,nDens)
*        Call RecPrt('D1Sao',' ',Work(ipD1Sao),nh1,1)
         call dcopy_(nh1,Work(ipD1Sao),1,Work(ip_D_DS+nh1),1)
         Call GetMem('DrvXV','Free','Real',ipD1Sao,nDens)
      End If
*
*---- Compute alpha and beta densities
*
*     Call RecPrt('DTot',' ',Work(ip_D_DS),nh1,1)
*     Call RecPrt('DSpn',' ',Work(ip_D_DS+nh1),nh1,1)
      If (nD.eq.1) Then
         Do i = 1, nh1
            DTot=Work(ip_D_DS+i-1)
            d_Alpha=Half*DTot
            Work(ip_D_DS+i-1)=    d_Alpha
         End Do
      Else
         Do i = 1, nh1
            DTot=Work(ip_D_DS+i-1)
            DSpn=Work(ip_D_DS+i-1+nh1)
            d_Alpha=Half*(DTot+DSpn)
            d_Beta =Half*(DTot-DSpn)
            Work(ip_D_DS+i-1)=    d_Alpha
            Work(ip_D_DS+i-1+nh1)=d_Beta
         End Do
      End If
*     Call RecPrt('Da',' ',Work(ip_D_DS),nh1,1)
*     Call RecPrt('Db',' ',Work(ip_D_DS+nh1),nh1,1)
*
      If(KSDFT(1:3).ne.'SCF') Then
        Call Get_iArray('nIsh',nIsh,mIrrep)
        Call Get_iArray('nFro',nFro,mIrrep)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     DFT functionals, compute integrals over the potential
*
      Func            =Zero
      Funcaa          =Zero
      Funcbb          =Zero
      Funccc          =Zero
      Dens_I          =Zero
      Dens_a1         =Zero
      Dens_b1         =Zero
      Dens_a2         =Zero
      Dens_b2         =Zero
      Dens_t1         =Zero
      Dens_t2         =Zero
      Grad_I          =Zero
      Tau_I           =Zero
      Do_MO           =.False.
      Do_TwoEl        =.False.
*
      If (nD.eq.2.and.DFTFOCK.eq.'DIFF') Then
         numAO=0
         If(KSDFT(1:3).ne.'SCF') Then
           Do iIrrep=0,mIrrep-1
             nAsh(iIrrep)=0
           End Do
           Call qpg_iArray('nAsh',Found,nData)
           If(Found .and. nData.eq.mIrrep) Then
             Call Get_iArray('nAsh',nAsh,mIrrep)
           End If
           Do iIrrep=0,mIrrep-1
             numAO=numAO+nAsh(iIrrep)
           End Do
        End If
        If (numAO.ne.0)
     &  Do_TwoEl        =.True.
        Do_MO           =.True.

      End If
*
*     nFckDim: number of different types of Fock matrices. Normally for
*     conventional functionals we have one Fock matrix for closed shell
*     calculations and two (F_alpha and F_beta) for open shell systems.
*     For CASDFT we have always two (F_inactive and F_active)
*
************************************************************************
*                                                                      *
*      LSDA LDA SVWN  GLM stuff                                        *
*                                                                      *
       If (KSDFT.eq.'LSDA ' .or.
     &     KSDFT.eq.'LDA '  .or.
     &     KSDFT.eq.'TLSDA'  .or. !GLM
     &     KSDFT.eq.'FTLSDA'  .or. !AMS
     &     KSDFT.eq.'SVWN ') Then
         If(KSDFT.eq.'TLSDA'
     &     .or.KSDFT.eq.'FTLSDA') Do_MO=.true. !GLM
         If(KSDFT.eq.'TLSDA'
     &     .or.KSDFT.eq.'FTLSDA') Do_TwoEl=.true. !GLM
         ExFac=Get_ExFac(KSDFT)
         Functional_type=LDA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(LSDA   ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
c          write(6,*) 'Func in drvdft :', Func
*                                                                      *
************************************************************************
*                                                                      *
*      LSDA5 LDA5 SVWN5                                                *
*                                                                      *
       Else If (KSDFT.eq.'LSDA5' .or.
     &          KSDFT.eq.'LDA5'  .or.
     &          KSDFT.eq.'TLSDA5 '.or.
     &          KSDFT.eq.'SVWN5') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=LDA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(LSDA5  ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
c         write(6,*) 'Func in drvdft :', Func
*                                                                      *
************************************************************************
*                                                                      *
*     HFB                                                              *
*                                                                      *
       Else If (KSDFT.eq.'HFB') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(HFB    ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     HFO                                                              *
*                                                                      *
       Else If (KSDFT.eq.'HFO') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(HFO    ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     HFG                                                              *
*                                                                      *
       Else If (KSDFT.eq.'HFG') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(HFG    ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     HFB86                                                            *
*                                                                      *
       Else If (KSDFT.eq.'HFB86') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(HFB86   ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*      HFS                                                             *
*                                                                      *
       Else If (KSDFT.eq.'HFS') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=LDA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(HFS    ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*      XALPHA                                                          *
*                                                                      *
       Else If (KSDFT.eq.'XALPHA') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=LDA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(XAlpha ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     Overlap                                                          *
*                                                                      *
      Else If (KSDFT.eq.'Overlap') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=LDA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(Overlap,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     NucAtt                                                           *
*                                                                      *
      Else If (KSDFT.eq.'NucAtt') Then
         ExFac=One
         Functional_type=LDA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(NucAtt,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     BWIG                                                             *
*                                                                      *
      Else If (KSDFT.eq.'BWIG') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(BWIG   ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     BLYP                                                             *
*                                                                      *
      Else If (KSDFT.eq.'BLYP'
     &       .or.  KSDFT.eq.'TBLYP' !GLM
     &       .or.  KSDFT.eq.'FTBLYP' !AMS
     &        ) Then
       If(KSDFT.eq.'TBLYP'
     &       .or. KSDFT.eq.'FTBLYP') Do_MO=.true.
       If(KSDFT.eq.'TBLYP'
     &       .or. KSDFT.eq.'FTBLYP') Do_TwoEl=.true.
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(BLYP   ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     OLYP                                                             *
*                                                                      *
      Else If (KSDFT.eq.'OLYP') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(OLYP   ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     KT3                                                              *
*                                                                      *
      Else If (KSDFT.eq.'KT3') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(KT3   ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     KT2                                                              *
*                                                                      *
      Else If (KSDFT.eq.'KT2') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(KT2   ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     GLYP                                                             *
*                                                                      *
      Else If (KSDFT.eq.'GLYP') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(GLYP   ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     B86LYP                                                           *
*                                                                      *
      Else If (KSDFT.eq.'B86LYP') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(B86LYP   ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     BPBE                                                             *
*                                                                      *
      Else If (KSDFT.eq.'BPBE') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(BPBE   ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     OPBE                                                             *
*                                                                      *
      Else If (KSDFT.eq.'OPBE') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(OPBE   ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     GPBE                                                             *
*                                                                      *
      Else If (KSDFT.eq.'GPBE') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(GPBE   ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     B86PBE                                                           *
*                                                                      *
      Else If (KSDFT.eq.'B86PBE') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(B86PBE  ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     TLYP                                                             *
*                                                                      *
      Else If (KSDFT.eq.'TLYP') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(TLYP   ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     NLYP                                                             *
*                                                                      *
      Else If (KSDFT.eq.'NLYP') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(NLYP   ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     B3LYP                                                            *
*                                                                      *
      Else If (KSDFT.eq.'B3LYP ') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(B3LYP  ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     O3LYP                                                            *
*                                                                      *
      Else If (KSDFT.eq.'O3LYP ') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(O3LYP  ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     B2PLYP                                                           *
*                                                                      *
      Else If (KSDFT.eq.'B2PLYP') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(B2PLYP  ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     O2PLYP                                                           *
*                                                                      *
      Else If (KSDFT.eq.'O2PLYP') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(O2PLYP  ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     B3LYP5                                                           *
*                                                                      *
      Else If (KSDFT.eq.'B3LYP5') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(B3LYP5 ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     PBE                                                              *
*                                                                      *
      Else If (KSDFT.eq.'PBE'
     &     .or.KSDFT.eq.'TPBE' !GLM
     &     .or.KSDFT.eq.'FTPBE'!AMS
     &         ) then
         If(KSDFT.eq.'TPBE'
     &     .or.KSDFT.eq.'FTPBE') Do_MO=.true.
         If(KSDFT.eq.'TPBE'
     &     .or.KSDFT.eq.'FTPBE') Do_TwoEl=.true.
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(PBE   ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     revPBE                                                           *
*                                                                      *
      Else If (KSDFT.eq.'REVPBE'
     &     .or.KSDFT.eq.'TREVPBE' !GLM
     &     .or.KSDFT.eq.'FTREVPBE'!AMS
     &         ) then
         If(KSDFT.eq.'TREVPBE'
     &     .or.KSDFT.eq.'FTREVBPE') Do_MO=.true.
         If(KSDFT.eq.'TREVPBE'
     &     .or.KSDFT.eq.'FTREVPBE') Do_TwoEl=.true.
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(REVPBE,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     SSB                                                              *
*                                                                      *
      Else If (KSDFT.eq.'SSB'
     &     .or.KSDFT.eq.'TSSB') Then !GLM
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(SSB   ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     S12H                                                             *
*                                                                      *
      Else If (KSDFT.eq.'S12H') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(S12H  ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     S12G                                                             *
*                                                                      *
      Else If (KSDFT.eq.'S12G') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(S12G  ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     PBEsol                                                           *
*                                                                      *
      Else If (KSDFT.eq.'PBESOL') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(PBEsol   ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     RGE2                                                             *
*                                                                      *
      Else If (KSDFT.eq.'RGE2') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(RGE2   ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     PTCA                                                             *
*                                                                      *
      Else If (KSDFT.eq.'PTCA') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(PTCA   ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     PBE0                                                             *
*                                                                      *
      Else If (KSDFT.eq.'PBE0') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(PBE0  ,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     M06-L                                                            *
*                                                                      *
      Else If (KSDFT.eq.'M06L') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=meta_GGA_type1
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(M06L,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     M06                                                              *
*                                                                      *
      Else If (KSDFT.eq.'M06 ') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=meta_GGA_type1
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(M06,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     M06-2X                                                           *
*                                                                      *
      Else If (KSDFT.eq.'M062X') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=meta_GGA_type1
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(M062X,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     M06-HF                                                           *
*                                                                      *
      Else If (KSDFT.eq.'M06HF') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=meta_GGA_type1
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(M06HF,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*     Checker                                                          *
*                                                                      *
      Else If (KSDFT.eq.'CHECKER') Then
         ExFac=Zero
         Functional_type=meta_GGA_type2
         nFckDim = nD
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         Call DrvNQ(Checker,Work(ipF_DFT),nFckDim,Func,
     &              Work(ip_D_DS),nh1,nD,
     &              Do_Grad,
     &              Grad,nGrad,
     &              Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
*     CASDFT functionals:                                              *
*                                                                      *
      Else If (KSDFT(1:4).eq.'NEWF') Then
*                                                                      *
*        These functionals are still under construction.               *
*        The code, written by S.G. & C., will be now modified by       *
*                      Giovanni Ghigo (CGG)                            *
*                                                                      *
         If (DFTFOCK.ne.'DIFF') Then
            Call WarningMessage(2,
     &                 ' This is CASDFT type functional !!!;'
     &               //' You cannot use it in this calculation.')
            Call Quit_OnUserError()
         End If
         Do_Grad  = .False.
         Do_MO    = .True.
         Do_TwoEl = .True.
*
         ExFac=Get_ExFac(KSDFT)
         Functional_type=CASDFT_type
         nFckDim = 2
         Call Allocate_Work(ipF_DFT,nh1*nFckDim)
         Call FZero(Work(ipF_DFT),nh1*nFckDim)
         If ( KSDFT(5:5).eq.'0' )
     &      Call DrvNQ(NEWF ,Work(ipF_DFT),nFckDim,Func,
     &                 Work(ip_D_DS),nh1,nD,Do_Grad,
     &                 Grad,nGrad,
     &                 Do_MO,Do_TwoEl,DFTFOCK)
         If ( KSDFT(5:5).eq.'1' )
     &      Call DrvNQ(NEWF1 ,Work(ipF_DFT),nFckDim,Func,
     &                 Work(ip_D_DS),nh1,nD,Do_Grad,
     &                 Grad,nGrad,
     &                 Do_MO,Do_TwoEl,DFTFOCK)
*                                                                      *
************************************************************************
*                                                                      *
      Else
         Call WarningMessage(2,
     &               ' DrvDFT: Undefined functional type!')
         Write (6,*) '         Functional=',KSDFT(1:lKSDFT)
         Call Quit_OnUserError()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Energy_integrated=Func
*                                                                      *
************************************************************************
*                                                                      *
      If (KSDFT.eq.'Overlap'.or.KSDFT.eq.'NucAtt') Then
         call dcopy_(nh1,Work(ipF_DFT),1,h1,1)
         If (KSDFT.eq.'NucAtt') Energy_integrated=Func
      Else
*
*        Put out the integrated DFT energy and the DFT Fock matrices
*        on the RUNFILE
*
*        Call Put_DFT_Energy(Energy_integrated)
         Call Poke_dScalar('KSDFT energy',Energy_integrated)
         Call Put_dScalar('CASDFT energy',Energy_integrated)
         Call Put_dExcdRa(Work(ipF_DFT),nFckDim*nh1)
*         Write(6,'(a,f22.16)') " Energy in drvdft ",Energy_integrated
#ifdef _DEBUG_
         Write(6,'(a,f22.16)') " Energy ",Energy_integrated
         If (nFckDim.eq.1) Then
            Do i=1,nh1
               Write(6,'(i4,f22.16)') i,Work(ipF_DFT+i-1)
            End Do
         Else
            Do i=1,nh1
              Write(6,'(i4,3f22.16)') i,Work(ipF_DFT+i-1),
     &                                  Work(ipF_DFT+i-1+nh1),
     &        (Work(ipF_DFT+i-1)+Work(ipF_DFT+i-1+nh1))/2.0d0
            End Do
         End If
#endif

*
*        In the SCF program (traclc.f) the program computes the trace
*        of the one-electron hamiltonian over a set of densities. The
*        DFT contribution is not linear with respect to variations of
*        the density. However, with the following term we can include
*        the linear component in that code.
*
         Fact = Two
         If (nD.ne.1) Fact=One
         Vxc_ref(1)=Fact*DDot_(nh1,Work(ipF_DFT),1,Work(ip_D_DS),1)
         If (nD.ne.1) Then
           Vxc_ref(2)=DDot_(nh1,Work(ipF_DFT+nh1),1,Work(ip_D_DS+nh1),1)
         Else
            Vxc_ref(2)=Zero
         End If
         Call Put_Temp('Vxc_ref ',Vxc_ref,2)
      End If
*
      Call Free_Work(ipF_DFT)
      Call GetMem('D-DS','Free','Real',ip_D_DS,2*nh1)
      Call Free_iSD()
      Call QExit('DrvDFT')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(TwoHam)
         Call Unused_real_array(D)
         Call Unused_real(RepNuc)
         Call Unused_logical(First)
         Call Unused_logical(Dff)
         Call Unused_logical(lRF)
         Call Unused_real_array(D1I)
         Call Unused_real_array(D1A)
      End If
      End
