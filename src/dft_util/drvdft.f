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
         Sub => LSDA
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
         Sub => LSDA5
c         write(6,*) 'Func in drvdft :', Func
*                                                                      *
************************************************************************
*                                                                      *
*     HFB                                                              *
*                                                                      *
       Else If (KSDFT.eq.'HFB') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => HFB
*                                                                      *
************************************************************************
*                                                                      *
*     HFO                                                              *
*                                                                      *
       Else If (KSDFT.eq.'HFO') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => HFO
*                                                                      *
************************************************************************
*                                                                      *
*     HFG                                                              *
*                                                                      *
       Else If (KSDFT.eq.'HFG') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => HFG
*                                                                      *
************************************************************************
*                                                                      *
*     HFB86                                                            *
*                                                                      *
       Else If (KSDFT.eq.'HFB86') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => HFB86
*                                                                      *
************************************************************************
*                                                                      *
*      HFS                                                             *
*                                                                      *
       Else If (KSDFT.eq.'HFS') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=LDA_type
         Sub => HFS
*                                                                      *
************************************************************************
*                                                                      *
*      XALPHA                                                          *
*                                                                      *
       Else If (KSDFT.eq.'XALPHA') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=LDA_type
         Sub => XAlpha
*                                                                      *
************************************************************************
*                                                                      *
*     Overlap                                                          *
*                                                                      *
      Else If (KSDFT.eq.'Overlap') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=LDA_type
         Sub => Overlap
*                                                                      *
************************************************************************
*                                                                      *
*     NucAtt                                                           *
*                                                                      *
      Else If (KSDFT.eq.'NucAtt') Then
         ExFac=One
         Functional_type=LDA_type
         Sub => NucAtt
*                                                                      *
************************************************************************
*                                                                      *
*     BWIG                                                             *
*                                                                      *
      Else If (KSDFT.eq.'BWIG') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => BWIG
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
         Sub => BLYP
*                                                                      *
************************************************************************
*                                                                      *
*     OLYP                                                             *
*                                                                      *
      Else If (KSDFT.eq.'OLYP') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => OLYP
*                                                                      *
************************************************************************
*                                                                      *
*     KT3                                                              *
*                                                                      *
      Else If (KSDFT.eq.'KT3') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => KT3
*                                                                      *
************************************************************************
*                                                                      *
*     KT2                                                              *
*                                                                      *
      Else If (KSDFT.eq.'KT2') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => KT2
*                                                                      *
************************************************************************
*                                                                      *
*     GLYP                                                             *
*                                                                      *
      Else If (KSDFT.eq.'GLYP') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => GLYP
*                                                                      *
************************************************************************
*                                                                      *
*     B86LYP                                                           *
*                                                                      *
      Else If (KSDFT.eq.'B86LYP') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => B86LYP
*                                                                      *
************************************************************************
*                                                                      *
*     BPBE                                                             *
*                                                                      *
      Else If (KSDFT.eq.'BPBE') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => BPBE
*                                                                      *
************************************************************************
*                                                                      *
*     OPBE                                                             *
*                                                                      *
      Else If (KSDFT.eq.'OPBE'
     &     .or.KSDFT.eq.'TOPBE' !GLM
     &     .or.KSDFT.eq.'FTOPBE'!AMS
     &         ) then
         If(KSDFT.eq.'TOPBE'
     &     .or.KSDFT.eq.'FTOPBE') Do_MO=.true.
         If(KSDFT.eq.'TOPBE'
     &     .or.KSDFT.eq.'FTOPBE') Do_TwoEl=.true.
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => OPBE
*                                                                      *
************************************************************************
*                                                                      *
*     GPBE                                                             *
*                                                                      *
      Else If (KSDFT.eq.'GPBE') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => GPBE
*                                                                      *
************************************************************************
*                                                                      *
*     B86PBE                                                           *
*                                                                      *
      Else If (KSDFT.eq.'B86PBE') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => B86PBE
*                                                                      *
************************************************************************
*                                                                      *
*     TLYP                                                             *
*                                                                      *
      Else If (KSDFT.eq.'TLYP') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => TLYP
*                                                                      *
************************************************************************
*                                                                      *
*     B3LYP                                                            *
*                                                                      *
      Else If (KSDFT.eq.'B3LYP ') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => B3LYP
*                                                                      *
************************************************************************
*                                                                      *
*     O3LYP                                                            *
*                                                                      *
      Else If (KSDFT.eq.'O3LYP ') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => O3LYP
*                                                                      *
************************************************************************
*                                                                      *
*     B2PLYP                                                           *
*                                                                      *
      Else If (KSDFT.eq.'B2PLYP') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => B2PLYP
*                                                                      *
************************************************************************
*                                                                      *
*     O2PLYP                                                           *
*                                                                      *
      Else If (KSDFT.eq.'O2PLYP') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => O2PLYP
*                                                                      *
************************************************************************
*                                                                      *
*     B3LYP5                                                           *
*                                                                      *
      Else If (KSDFT.eq.'B3LYP5') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => B3LYP5
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
         Sub => PBE
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
         Sub => REVPBE
*                                                                      *
************************************************************************
*                                                                      *
*     SSBSW                                                              *
*                                                                      *
      Else If (KSDFT.eq.'SSBSW'
     &     .or.KSDFT.eq.'TSSBSW') Then !GLM
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => SSBSW
*                                                                      *
************************************************************************
*                                                                      *
*     SSBD                                                             *
*                                                                      *
      Else If (KSDFT.eq.'SSBD'
     &     .or.KSDFT.eq.'TSSBD') Then !GLM
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => SSBD
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
*     S12H                                                             *
*                                                                      *
      Else If (KSDFT.eq.'S12H') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => S12H
*                                                                      *
************************************************************************
*                                                                      *
*     S12G                                                             *
*                                                                      *
      Else If (KSDFT.eq.'S12G'
     &     .or.KSDFT.eq.'TS12G') Then !GLM
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => S12G
*                                                                      *
************************************************************************
*                                                                      *
*     PBEsol                                                           *
*                                                                      *
      Else If (KSDFT.eq.'PBESOL') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => PBESol
*                                                                      *
************************************************************************
*                                                                      *
*     RGE2                                                             *
*                                                                      *
      Else If (KSDFT.eq.'RGE2') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => RGE2
*                                                                      *
************************************************************************
*                                                                      *
*     PTCA                                                             *
*                                                                      *
      Else If (KSDFT.eq.'PTCA') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => PTCA
*                                                                      *
************************************************************************
*                                                                      *
*     PBE0                                                             *
*                                                                      *
      Else If (KSDFT.eq.'PBE0') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=GGA_type
         Sub => PBE0
*                                                                      *
************************************************************************
*                                                                      *
*     M06-L                                                            *
*                                                                      *
      Else If (KSDFT.eq.'M06L') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=meta_GGA_type1
         Sub => M06L
*                                                                      *
************************************************************************
*                                                                      *
*     M06                                                              *
*                                                                      *
      Else If (KSDFT.eq.'M06 ') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=meta_GGA_type1
         Sub => M06
*                                                                      *
************************************************************************
*                                                                      *
*     M06-2X                                                           *
*                                                                      *
      Else If (KSDFT.eq.'M062X') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=meta_GGA_type1
         Sub => M062X
*                                                                      *
************************************************************************
*                                                                      *
*     M06-HF                                                           *
*                                                                      *
      Else If (KSDFT.eq.'M06HF') Then
         ExFac=Get_ExFac(KSDFT)
         Functional_type=meta_GGA_type1
         Sub => M06HF
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
