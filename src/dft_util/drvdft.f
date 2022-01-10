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
* Copyright (C) 2022, Roland Lindh                                     *
************************************************************************
      Subroutine DrvDFT(h1,TwoHam,D,RepNuc,nh1,First,Dff,
     &                  lRF,KSDFT,ExFac,Do_Grad,Grad,nGrad,iSpin,
     &                  D1I,D1A,nD1,DFTFOCK)
      use KSDFT_Info, only: KSDFA, funcaa, funcbb, funccc
      Implicit Real*8 (a-h,o-z)
      External LSDA, Overlap, BLYP, BPBE, B3LYP, HFS, HFB,
     &         XAlpha, LSDA5, B3LYP5, B2PLYP, TLYP,
     &         NucAtt, OLYP, O3LYP, OPBE,
     &         PBE, PBE0, PBEsol, M06L, M06, M062X, HFO,
     &         M06HF, SSBSW, SSBD, HFG, GLYP, GPBE,
     &         HFB86, B86LYP, B86PBE, BWIG, KT3,
     &         O2PLYP,  KT2,  RGE2, REVPBE,
     &         PTCA,S12G, S12H
#include "real.fh"
#include "stdalloc.fh"
#include "nq_info.fh"
#include "debug.fh"
#include "pamint.fh"
#include "ksdft.fh"
      Real*8 h1(nh1), TwoHam(nh1), D(nh1,2), Grad(nGrad), Vxc_ref(2)
      Real*8 D1I(nD1),D1A(nD1)
      Logical First, Dff, lRF,  Do_Grad
      Logical Do_MO,Do_TwoEl
      Character*(*) KSDFT
      Character*4 DFTFOCK
      Real*8, Allocatable:: D_DS(:,:), F_DFT(:,:)

      abstract interface
          Subroutine DFT_FUNCTIONAL(mGrid,nD)
          Integer mGrid, nD
          end subroutine
      end interface

      procedure(DFT_FUNCTIONAL), pointer :: sub => null()

*
      KSDFA = KSDFT
      lKSDFT=LEN(KSDFT)
      Debug=.False.
*                                                                      *
************************************************************************
*                                                                      *
c     Call SetQue('Trace=on')
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
      Call Get_dScalar('DFT exch coeff',CoefX)
      Call Get_dScalar('DFT corr coeff',CoefR)
*
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
      If (DFTFOCK.eq.'ROKS') nD=2
      Call mma_allocate(D_DS,nh1,nD,Label='D_DS')
*
*---- Get the total density
*
      Call Get_D1ao(D_DS,nh1)
*     Call RecPrt('D1ao',' ',D_DS(:,1),nh1,1)
*
*
*---- Get the spin density
*
      If (nD.ne.1) Then
         Call Get_D1Sao(D_DS(:,2),nh1)
*        Call RecPrt('D1Sao',' ',D_DS(:,2),nh1,1)
      End If
*
*---- Compute alpha and beta densities
*
*     Call RecPrt('DTot',' ',D_DS(:,1),nh1,1)
*     Call RecPrt('DSpn',' ',D_DS(:,2),nh1,1)
      If (nD.eq.1) Then
        D_DS(:,1)=Half*D_DS(:,1)
      Else
         Do i = 1, nh1
            DTot=D_DS(i,1)
            DSpn=D_DS(i,2)
            d_Alpha=Half*(DTot+DSpn)
            d_Beta =Half*(DTot-DSpn)
            D_DS(i,1)=d_Alpha
            D_DS(i,2)=d_Beta
         End Do
      End If
*     Call RecPrt('Da',' ',D_DS(:,1),nh1,1)
*     Call RecPrt('Db',' ',D_DS(:,2),nh1,1)
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
*     nFckDim: number of different types of Fock matrices. Normally for
*     conventional functionals we have one Fock matrix for closed shell
*     calculations and two (F_alpha and F_beta) for open shell systems.
*     For CASDFT we have always two (F_inactive and F_active)
*
*                                                                      *
************************************************************************
*                                                                      *
       Select Case(KSDFT)
*                                                                      *
************************************************************************
*                                                                      *
*      LSDA LDA SVWN                                                   *
*                                                                      *
      Case('LSDA ','LDA ','TLSDA','FTLSDA','SVWN ')
         Functional_type=LDA_type
         Sub => LSDA
         If(KSDFT.eq.'TLSDA'
     &     .or.KSDFT.eq.'FTLSDA') Do_MO=.true.
         If(KSDFT.eq.'TLSDA'
     &     .or.KSDFT.eq.'FTLSDA') Do_TwoEl=.true.

*                                                                      *
************************************************************************
*                                                                      *
*      LSDA5 LDA5 SVWN5                                                *
*                                                                      *
       Case('LSDA5','LDA5','TLSDA5 ','SVWN5')
         Functional_type=LDA_type
         Sub => LSDA5
*                                                                      *
************************************************************************
*                                                                      *
*     HFB                                                              *
*                                                                      *
       Case('HFB')
         Functional_type=GGA_type
         Sub => HFB
*                                                                      *
************************************************************************
*                                                                      *
*     HFO                                                              *
*                                                                      *
       Case('HFO')
         Functional_type=GGA_type
         Sub => HFO
*                                                                      *
************************************************************************
*                                                                      *
*     HFG                                                              *
*                                                                      *
       Case('HFG')
         Functional_type=GGA_type
         Sub => HFG
*                                                                      *
************************************************************************
*                                                                      *
*     HFB86                                                            *
*                                                                      *
       Case('HFB86')
         Functional_type=GGA_type
         Sub => HFB86
*                                                                      *
************************************************************************
*                                                                      *
*      HFS                                                             *
*                                                                      *
       Case('HFS')
         Functional_type=LDA_type
         Sub => HFS
*                                                                      *
************************************************************************
*                                                                      *
*      XALPHA                                                          *
*                                                                      *
       Case('XALPHA')
         Functional_type=LDA_type
         Sub => XAlpha
*                                                                      *
************************************************************************
*                                                                      *
*     Overlap                                                          *
*                                                                      *
      Case('Overlap')
         Functional_type=LDA_type
         Sub => Overlap
*                                                                      *
************************************************************************
*                                                                      *
*     NucAtt                                                           *
*                                                                      *
      Case('NucAtt')
         Functional_type=LDA_type
         Sub => NucAtt
*                                                                      *
************************************************************************
*                                                                      *
*     BWIG                                                             *
*                                                                      *
      Case('BWIG')
         Functional_type=GGA_type
         Sub => BWIG
*                                                                      *
************************************************************************
*                                                                      *
*     BLYP                                                             *
*                                                                      *
      Case('BLYP','TBLYP','FTBLYP')
         Functional_type=GGA_type
         Sub => BLYP
       If(KSDFT.eq.'TBLYP'
     &       .or. KSDFT.eq.'FTBLYP') Do_MO=.true.
       If(KSDFT.eq.'TBLYP'
     &       .or. KSDFT.eq.'FTBLYP') Do_TwoEl=.true.
*                                                                      *
************************************************************************
*                                                                      *
*     OLYP                                                             *
*                                                                      *
      Case ('OLYP')
         Functional_type=GGA_type
         Sub => OLYP
*                                                                      *
************************************************************************
*                                                                      *
*     KT3                                                              *
*                                                                      *
      Case('KT3')
         Functional_type=GGA_type
         Sub => KT3
*                                                                      *
************************************************************************
*                                                                      *
*     KT2                                                              *
*                                                                      *
      Case('KT2')
         Functional_type=GGA_type
         Sub => KT2
*                                                                      *
************************************************************************
*                                                                      *
*     GLYP                                                             *
*                                                                      *
      Case('GLYP')
         Functional_type=GGA_type
         Sub => GLYP
*                                                                      *
************************************************************************
*                                                                      *
*     B86LYP                                                           *
*                                                                      *
      Case('B86LYP')
         Functional_type=GGA_type
         Sub => B86LYP
*                                                                      *
************************************************************************
*                                                                      *
*     BPBE                                                             *
*                                                                      *
      Case('BPBE')
         Functional_type=GGA_type
         Sub => BPBE
*                                                                      *
************************************************************************
*                                                                      *
*     OPBE                                                             *
*                                                                      *
      Case('OPBE','TOPBE','FTOPBE')
         Functional_type=GGA_type
         Sub => OPBE
         If(KSDFT.eq.'TOPBE'
     &     .or.KSDFT.eq.'FTOPBE') Do_MO=.true.
         If(KSDFT.eq.'TOPBE'
     &     .or.KSDFT.eq.'FTOPBE') Do_TwoEl=.true.
*                                                                      *
************************************************************************
*                                                                      *
*     GPBE                                                             *
*                                                                      *
      Case('GPBE')
         Functional_type=GGA_type
         Sub => GPBE
*                                                                      *
************************************************************************
*                                                                      *
*     B86PBE                                                           *
*                                                                      *
      Case('B86PBE')
         Functional_type=GGA_type
         Sub => B86PBE
*                                                                      *
************************************************************************
*                                                                      *
*     TLYP                                                             *
*                                                                      *
      Case('TLYP')
         Functional_type=GGA_type
         Sub => TLYP
*                                                                      *
************************************************************************
*                                                                      *
*     B3LYP                                                            *
*                                                                      *
      Case('B3LYP ')
         Functional_type=GGA_type
         Sub => B3LYP
*                                                                      *
************************************************************************
*                                                                      *
*     O3LYP                                                            *
*                                                                      *
      Case('O3LYP ')
         Functional_type=GGA_type
         Sub => O3LYP
*                                                                      *
************************************************************************
*                                                                      *
*     B2PLYP                                                           *
*                                                                      *
      Case('B2PLYP')
         Functional_type=GGA_type
         Sub => B2PLYP
*                                                                      *
************************************************************************
*                                                                      *
*     O2PLYP                                                           *
*                                                                      *
      Case('O2PLYP')
         Functional_type=GGA_type
         Sub => O2PLYP
*                                                                      *
************************************************************************
*                                                                      *
*     B3LYP5                                                           *
*                                                                      *
      Case('B3LYP5')
         Functional_type=GGA_type
         Sub => B3LYP5
*                                                                      *
************************************************************************
*                                                                      *
*     PBE                                                              *
*                                                                      *
      Case('PBE','TPBE','FTPBE')
         Functional_type=GGA_type
         Sub => PBE
         If(KSDFT.eq.'TPBE'
     &     .or.KSDFT.eq.'FTPBE') Do_MO=.true.
         If(KSDFT.eq.'TPBE'
     &     .or.KSDFT.eq.'FTPBE') Do_TwoEl=.true.
*                                                                      *
************************************************************************
*                                                                      *
*     revPBE                                                           *
*                                                                      *
      Case('REVPBE','TREVPBE','FTREVPBE')
         Functional_type=GGA_type
         Sub => REVPBE
         If(KSDFT.eq.'TREVPBE'
     &     .or.KSDFT.eq.'FTREVBPE') Do_MO=.true.
         If(KSDFT.eq.'TREVPBE'
     &     .or.KSDFT.eq.'FTREVPBE') Do_TwoEl=.true.
*                                                                      *
************************************************************************
*                                                                      *
*     SSBSW                                                              *
*                                                                      *
      Case('SSBSW','TSSBSW')
         Functional_type=GGA_type
         Sub => SSBSW
*                                                                      *
************************************************************************
*                                                                      *
*     SSBD                                                             *
*                                                                      *
      Case('SSBD','TSSBD')
         Functional_type=GGA_type
         Sub => SSBD
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
*     S12H                                                             *
*                                                                      *
      Case('S12H')
         Functional_type=GGA_type
         Sub => S12H
*                                                                      *
************************************************************************
*                                                                      *
*     S12G                                                             *
*                                                                      *
      Case('S12G','TS12G')
         Functional_type=GGA_type
         Sub => S12G
*                                                                      *
************************************************************************
*                                                                      *
*     PBEsol                                                           *
*                                                                      *
      Case('PBESOL')
         Functional_type=GGA_type
         Sub => PBESol
*                                                                      *
************************************************************************
*                                                                      *
*     RGE2                                                             *
*                                                                      *
      Case('RGE2')
         Functional_type=GGA_type
         Sub => RGE2
*                                                                      *
************************************************************************
*                                                                      *
*     PTCA                                                             *
*                                                                      *
      Case('PTCA')
         Functional_type=GGA_type
         Sub => PTCA
*                                                                      *
************************************************************************
*                                                                      *
*     PBE0                                                             *
*                                                                      *
      Case('PBE0')
         Functional_type=GGA_type
         Sub => PBE0
*                                                                      *
************************************************************************
*                                                                      *
*     M06-L                                                            *
*                                                                      *
      Case('M06L')
         Functional_type=meta_GGA_type1
         Sub => M06L
*                                                                      *
************************************************************************
*                                                                      *
*     M06                                                              *
*                                                                      *
      Case('M06 ')
         Functional_type=meta_GGA_type1
         Sub => M06
*                                                                      *
************************************************************************
*                                                                      *
*     M06-2X                                                           *
*                                                                      *
      Case('M062X')
         Functional_type=meta_GGA_type1
         Sub => M062X
*                                                                      *
************************************************************************
*                                                                      *
*     M06-HF                                                           *
*                                                                      *
      Case('M06HF')
         Functional_type=meta_GGA_type1
         Sub => M06HF
*                                                                      *
************************************************************************
*                                                                      *
      Case default
         Call WarningMessage(2,
     &               ' DrvDFT: Undefined functional type!')
         Write (6,*) '         Functional=',KSDFT(1:lKSDFT)
         Call Quit_OnUserError()
       End Select
*                                                                      *
************************************************************************
*                                                                      *
      ExFac=Get_ExFac(KSDFT)
      nFckDim = nD
      Call mma_allocate(F_DFT,nh1,nFckDim,Label='F_DFT')
      F_DFT(:,:)=Zero
      Call DrvNQ(Sub,F_DFT,nFckDim,Func,
     &           D_DS,nh1,nD,
     &           Do_Grad,
     &           Grad,nGrad,
     &           Do_MO,Do_TwoEl,DFTFOCK)

      Sub => Null()
*                                                                      *
************************************************************************
*                                                                      *
      Energy_integrated=Func
*                                                                      *
************************************************************************
*                                                                      *
      If (KSDFT.eq.'Overlap'.or.KSDFT.eq.'NucAtt') Then
         call dcopy_(nh1,F_DFT,1,h1,1)
         If (KSDFT.eq.'NucAtt') Energy_integrated=Func
      Else
*
*        Put out the integrated DFT energy and the DFT Fock matrices
*        on the RUNFILE
*
*        Call Put_DFT_Energy(Energy_integrated)
         Call Poke_dScalar('KSDFT energy',Energy_integrated)
         Call Put_dScalar('CASDFT energy',Energy_integrated)
         Call Put_dExcdRa(F_DFT,nFckDim*nh1)
*         Write(6,'(a,f22.16)') " Energy in drvdft ",Energy_integrated
#ifdef _DEBUGPRINT_
         Write(6,'(a,f22.16)') " Energy ",Energy_integrated
         If (nFckDim.eq.1) Then
            Do i=1,nh1
               Write(6,'(i4,f22.16)') i,F_DFT(i,1)
            End Do
         Else
            Do i=1,nh1
              Write(6,'(i4,3f22.16)') i,F_DFT(i,1),
     &                                  F_DFT(i,2),
     &        F_DFT(i,1)+F_DFT(i,2)/2.0d0
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
         Vxc_ref(1)=Fact*DDot_(nh1,F_DFT(:,1),1,D_DS,1)
         If (nD.ne.1) Then
           Vxc_ref(2)=DDot_(nh1,F_DFT(:,2),1,D_DS(:,2),1)
         Else
            Vxc_ref(2)=Zero
         End If
         Call Put_Temp('Vxc_ref ',Vxc_ref,2)
      End If
*
      Call mma_deallocate(F_DFT)
      Call mma_deallocate(D_DS)
      Call Free_iSD()
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
