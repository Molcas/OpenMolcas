!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2022, Roland Lindh                                     *
!***********************************************************************
Subroutine Driver(KSDFT,Do_Grad,Func,Grad,nGrad,Do_MO,Do_TwoEl,D_DS,F_DFT,nh1,nD,DFTFOCK)
      use libxc_parameters
      use OFembed, only: KEOnly, dFMD, Do_Core
      use libxc,   only: Only_exc
      Implicit None
#include "nq_info.fh"
      Character*(*) KSDFT
      Logical Do_Grad
      Integer :: nGrad, nh1, nD
      Real*8 :: Func, Grad(nGrad)
      Logical Do_MO, Do_TwoEl
      Real*8 :: D_DS(nh1,nD), F_DFT(nh1,nD)
      Character*4 DFTFOCK


      abstract interface
          Subroutine DFT_FUNCTIONAL(mGrid,nD)
          Integer mGrid, nD
          end subroutine
      end interface

!***********************************************************************
!     Define external functions not defined in LibXC. These are either
!     accessed through the procedure pointer sub or External_sub.

      External:: Overlap, NucAtt, ndsd_ts
!***********************************************************************
      procedure(DFT_FUNCTIONAL), pointer :: sub => null()
!     Sometime we need an external routine which covers something which
!     Libxc doesn't support.
      procedure(DFT_FUNCTIONAL), pointer :: External_sub => null()
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
!      Default is to use the libxc interface
!      Coefficient for the individual contibutions are defaulted to 1.0D0

       Sub => libxc_functionals     ! Default
       Coeffs(:)=1.0D0              ! Default
!                                                                      *
!***********************************************************************
!                                                                      *
       Select Case(KSDFT)
!                                                                      *
!***********************************************************************
!                                                                      *
!      LSDA LDA SVWN                                                   *
!                                                                      *
      Case('LSDA ','LDA ','TLSDA','FTLSDA','SVWN ')
         Functional_type=LDA_type

!----    Slater exchange
!
!----    Vosko-Wilk-Nusair correlation functional III

         If (Do_Core) Then
            nFuncs=1
            func_id(1:nFuncs)=[XC_LDA_C_VWN_RPA]
            Coeffs(1)=dFMD
         Else
            nFuncs=2
            func_id(1:nFuncs)=[XC_LDA_X,XC_LDA_C_VWN_RPA]
         End If

         If(KSDFT.eq.'TLSDA'.or.KSDFT.eq.'FTLSDA') Do_MO=.true.
         If(KSDFT.eq.'TLSDA'.or.KSDFT.eq.'FTLSDA') Do_TwoEl=.true.
!                                                                      *
!***********************************************************************
!                                                                      *
!      LSDA5 LDA5 SVWN5                                                *
!                                                                      *
       Case('LSDA5','LDA5','TLSDA5 ','SVWN5')
         Functional_type=LDA_type

!----    Slater exchange
!
!----    Vosko-Wilk-Nusair correlation functional V

         If (Do_Core) Then
            nFuncs=1
            func_id(1:nFuncs)=[XC_LDA_C_VWN]
            Coeffs(1)=dFMD
         Else
            nFuncs=2
            func_id(1:nFuncs)=[XC_LDA_X,XC_LDA_C_VWN]
         End If
!                                                                      *
!***********************************************************************
!                                                                      *
!     HFB                                                              *
!                                                                      *
       Case('HFB')
         Functional_type=GGA_type

!----    Slater exchange + Becke 88 exchange

         nFuncs=1
         func_id(1:nFuncs)=[XC_GGA_X_B88]
!                                                                      *
!***********************************************************************
!                                                                      *
!     HFO                                                              *
!                                                                      *
       Case('HFO')
         Functional_type=GGA_type

!----    Slater exchange + OPTx exchange

         nFuncs=1
         func_id(1:nFuncs)=[XC_GGA_X_OPTX]
!                                                                      *
!***********************************************************************
!                                                                      *
!     HFG                                                              *
!                                                                      *
       Case('HFG')
         Functional_type=GGA_type

!----    Slater exchange + G96 exchange

         nFuncs=1
         func_id(1:nFuncs)=[XC_GGA_X_G96]
!                                                                      *
!***********************************************************************
!                                                                      *
!     HFB86                                                            *
!                                                                      *
       Case('HFB86')
         Functional_type=GGA_type

!----    Slater exchange + B86 exchange

         nFuncs=1
         func_id(1:nFuncs)=[XC_GGA_X_B86]
!                                                                      *
!***********************************************************************
!                                                                      *
!      HFS                                                             *
!                                                                      *
       Case('HFS')
         Functional_type=LDA_type

!----    Slater exchange

         nFuncs=1
         func_id(1:nFuncs)=[XC_LDA_X]
!                                                                      *
!***********************************************************************
!                                                                      *
!      XALPHA                                                          *
!                                                                      *
       Case('XALPHA')
         Functional_type=LDA_type

!----    Slater exchange

         nFuncs=1
         func_id(1:nFuncs)=[XC_LDA_X]
         Coeffs(1)=0.70D0
!                                                                      *
!***********************************************************************
!                                                                      *
!     Overlap                                                          *
!                                                                      *
      Case('Overlap')
         Functional_type=LDA_type
         Sub => Overlap
!                                                                      *
!***********************************************************************
!                                                                      *
!     NucAtt                                                           *
!                                                                      *
      Case('NucAtt')
         Functional_type=LDA_type
         Sub => NucAtt
!                                                                      *
!***********************************************************************
!                                                                      *
!     BWIG                                                             *
!                                                                      *
      Case('BWIG')
         Functional_type=GGA_type

!----    Slater exchange + B88 exchange
!
!----    Wigner correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_B88,XC_LDA_C_OW_LYP]
!                                                                      *
!***********************************************************************
!                                                                      *
!     BLYP                                                             *
!                                                                      *
      Case('BLYP','TBLYP','FTBLYP')
         Functional_type=GGA_type

!----    Slater exchange + B88 exchange
!
!----    Lee-Yang-Parr correlation

         If (Do_Core) Then
            nFuncs=1
            func_id(1:nFuncs)=[XC_GGA_C_LYP]
            Coeffs(1)=dFMD
         Else
            nFuncs=2
            func_id(1:nFuncs)=[XC_GGA_X_B88,XC_GGA_C_LYP]
         End If

         If(KSDFT.eq.'TBLYP'.or. KSDFT.eq.'FTBLYP') Do_MO=.true.
         If(KSDFT.eq.'TBLYP'.or. KSDFT.eq.'FTBLYP') Do_TwoEl=.true.
!                                                                      *
!***********************************************************************
!                                                                      *
!     OLYP                                                             *
!                                                                      *
      Case ('OLYP')
         Functional_type=GGA_type

!----    Slater exchange + OPT exchange
!
!----    Lee-Yang-Parr correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_OPTX,XC_GGA_C_LYP]
!                                                                      *
!***********************************************************************
!                                                                      *
!     KT3                                                              *
!                                                                      *
      Case('KT3')
         Functional_type=GGA_type

!----    Slater exchange + + OPT exchange + Keal-Tozer exchange
!
!----    Lee-Yang-Parr correlation

         nFuncs=4
         func_id(1:nFuncs)=[XC_LDA_X,XC_GGA_X_OPTX,XC_GGA_X_KT1,XC_GGA_C_LYP]
         Coeffs(1)= (1.092d0-1.051510d0*(0.925452d0/1.431690d0)-(0.004d0/0.006d0))
         Coeffs(2)= (0.925452d0/1.431690d0)
         Coeffs(3)= (0.0040d0/0.006d0)
         Coeffs(4)= 0.864409d0

!                                                                      *
!***********************************************************************
!                                                                      *
!     KT2                                                              *
!                                                                      *
      Case('KT2')
         Functional_type=GGA_type

!----    Slater exchange + Keal-Tozer exchange
!
!----    Vosko-Wilk-Nusair correlation functional III

         nFuncs=3
         func_id(1:nFuncs)=[XC_LDA_X,XC_GGA_X_KT1,XC_LDA_C_VWN_RPA]
         Coeffs(1)= 0.07173d0
!        Coeffs(2)= 1.0D0
         Coeffs(3)= 0.576727d0
!                                                                      *
!***********************************************************************
!                                                                      *
!     GLYP                                                             *
!                                                                      *
      Case('GLYP')
         Functional_type=GGA_type

!----    Slater exchange + Gill 96 exchange
!
!----    Lee-Yang-Parr correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_G96,XC_GGA_C_LYP]
!                                                                      *
!***********************************************************************
!                                                                      *
!     B86LYP                                                           *
!                                                                      *
      Case('B86LYP')
         Functional_type=GGA_type

!----    Slater exchange + Becke 86 exchange
!
!----    Lee-Yang-Parr correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_B86,XC_GGA_C_LYP]
!                                                                      *
!***********************************************************************
!                                                                      *
!     BPBE                                                             *
!                                                                      *
      Case('BPBE')
         Functional_type=GGA_type

!----    Slater exchange + Becke 88 exchange
!
!----    Perdew-Burk-Ernzerhof correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_B88,XC_GGA_C_PBE]
!                                                                      *
!***********************************************************************
!                                                                      *
!     OPBE                                                             *
!                                                                      *
      Case('OPBE','TOPBE','FTOPBE')
         Functional_type=GGA_type

!----    Slater exchange + OPT exchange
!
!----    Perdew-Burk-Ernzerhof correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_OPTX,XC_GGA_C_PBE]

         If(KSDFT.eq.'TOPBE'.or.KSDFT.eq.'FTOPBE') Do_MO=.true.
         If(KSDFT.eq.'TOPBE'.or.KSDFT.eq.'FTOPBE') Do_TwoEl=.true.
!                                                                      *
!***********************************************************************
!                                                                      *
!     GPBE                                                             *
!                                                                      *
      Case('GPBE')
         Functional_type=GGA_type

!----    Slater exchange + Gill 96 exchange
!
!----    Perdew-Burk-Ernzerhof correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_G96,XC_GGA_C_PBE]
!                                                                      *
!***********************************************************************
!                                                                      *
!     B86PBE                                                           *
!                                                                      *
      Case('B86PBE')
         Functional_type=GGA_type

!----    Slater exchange + Becke 86 exchange
!
!----    Perdew-Burk-Ernzerhof correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_B86,XC_GGA_C_PBE]
!                                                                      *
!***********************************************************************
!                                                                      *
!     TLYP                                                             *
!                                                                      *
      Case('TLYP')
         Functional_type=GGA_type

!----    Lee-Yang-Parr correlation

         nFuncs=1
         func_id(1:nFuncs)=[XC_GGA_C_LYP]
!                                                                      *
!***********************************************************************
!                                                                      *
!     B3LYP                                                            *
!                                                                      *
      Case('B3LYP ')
         Functional_type=GGA_type

!----    Slater exchange + Becke 88 exchange
!
!----    Perdew-Burk-Ernzerhof correlation
!
!----    Vosko-Wilk-Nusair correlation functional III
!
!----    Lee-Yang-Parr correlation

         nFuncs=4
         func_id(1:nFuncs)=[XC_LDA_X,XC_GGA_X_B88,XC_LDA_C_VWN_RPA,XC_GGA_C_LYP]
         Coeffs(1)=0.08D0
         Coeffs(2)=0.72D0
         Coeffs(3)=1.0D0-0.81D0
         Coeffs(4)=0.81D0
!                                                                      *
!***********************************************************************
!                                                                      *
!     O3LYP                                                            *
!                                                                      *
      Case('O3LYP ')
         Functional_type=GGA_type

!----    Slater exchange + OPT exchange
!
!----    Perdew-Burk-Ernzerhof correlation
!
!----    Vosko-Wilk-Nusair correlation functional III
!
!----    Lee-Yang-Parr correlation

         nFuncs=4
         func_id(1:nFuncs)=[XC_LDA_X,XC_GGA_X_OPTX,XC_LDA_C_VWN_RPA,XC_GGA_C_LYP]
         Coeffs(1)=(0.9262D0-1.051510d0*(0.8133D0/1.431690d0))
         Coeffs(2)=(0.8133D0/1.431690d0)
         Coeffs(3)=1.0D0-0.81D0
         Coeffs(4)=0.81D0
!                                                                      *
!***********************************************************************
!                                                                      *
!     B2PLYP                                                           *
!                                                                      *
      Case('B2PLYP')
         Functional_type=GGA_type

!----    Slater exchange + Becke 88 exchange
!
!----    Lee-Yang-Parr correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_B88,XC_GGA_C_LYP]
         Coeffs(1)=0.470D0
         Coeffs(2)=0.73D0
!                                                                      *
!***********************************************************************
!                                                                      *
!     O2PLYP                                                           *
!                                                                      *
      Case('O2PLYP')
         Functional_type=GGA_type

!----    Slater exchange + OPT exchange
!
!----    OPT exchange
!
!----    Lee-Yang-Parr correlation

         nFuncs=3
         func_id(1:nFuncs)=[XC_LDA_X,XC_GGA_X_OPTX,XC_GGA_C_LYP]
         Coeffs(1)=(0.525755D0-1.051510d0*(0.715845D0/1.431690d0))
         Coeffs(2)=(0.715845D0/1.431690d0)
         Coeffs(3)=0.75D0
!                                                                      *
!***********************************************************************
!                                                                      *
!     B3LYP5                                                           *
!                                                                      *
      Case('B3LYP5')
         Functional_type=GGA_type

!----    Slater exchange + Becke 88 exchange
!
!----    Perdew-Burk-Ernzerhof correlation
!
!----    Vosko-Wilk-Nusair correlation functional V
!
!----    Lee-Yang-Parr correlation

         nFuncs=4
         func_id(1:nFuncs)=[XC_LDA_X,XC_GGA_X_B88,XC_LDA_C_VWN,XC_GGA_C_LYP]
         Coeffs(1)=0.08D0
         Coeffs(2)=0.72D0
         Coeffs(3)=1.0D0-0.81D0
         Coeffs(4)=0.81D0
!                                                                      *
!***********************************************************************
!                                                                      *
!     PBE                                                              *
!                                                                      *
      Case('PBE','TPBE','FTPBE')
         Functional_type=GGA_type

!----    Perdew-Burk-Ernzerhof exchange
!
!----    Perdew-Burk-Ernzerhof correlation

         If (Do_Core) Then
            nFuncs=1
            func_id(1:nFuncs)=[XC_GGA_C_PBE]
            Coeffs(1)=dFMD
         Else
            nFuncs=2
            func_id(1:nFuncs)=[XC_GGA_X_PBE,XC_GGA_C_PBE]
         End If

         If(KSDFT.eq.'TPBE'.or.KSDFT.eq.'FTPBE') Do_MO=.true.
         If(KSDFT.eq.'TPBE'.or.KSDFT.eq.'FTPBE') Do_TwoEl=.true.
!                                                                      *
!***********************************************************************
!                                                                      *
!     revPBE                                                           *
!                                                                      *
      Case('REVPBE','TREVPBE','FTREVPBE')
         Functional_type=GGA_type

!----    Revised Perdew-Burk-Ernzerhof exchange
!
!----    Perdew-Burk-Ernzerhof correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_PBE_R,XC_GGA_C_PBE]

         If(KSDFT.eq.'TREVPBE'.or.KSDFT.eq.'FTREVBPE') Do_MO=.true.
         If(KSDFT.eq.'TREVPBE'.or.KSDFT.eq.'FTREVPBE') Do_TwoEl=.true.
!                                                                      *
!***********************************************************************
!                                                                      *
!     SSBSW                                                              *
!                                                                      *
      Case('SSBSW','TSSBSW')
         Functional_type=GGA_type

!----    Swarts-Solo-Bickelhaupt sw exchange
!
!----    Perdew-Burk-Ernzerhof correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_SSB_SW,XC_GGA_C_PBE]
!                                                                      *
!***********************************************************************
!                                                                      *
!     SSBD                                                             *
!                                                                      *
      Case('SSBD','TSSBD')
         Functional_type=GGA_type

!----    Swarts-Solo-Bickelhaupt-D exchange
!
!----    Perdew-Burk-Ernzerhof correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_SSB_D,XC_GGA_C_PBE]
!                                                                      *
!                                                                      *
!***********************************************************************
!                                                                      *
!     S12H                                                             *
!                                                                      *
      Case('S12H')
         Functional_type=GGA_type

!----    Swarts hybrid exchange
!
!----    Perdew-Burk-Ernzerhof correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_HYB_GGA_X_S12H,XC_GGA_C_PBE]
!                                                                      *
!***********************************************************************
!                                                                      *
!     S12G                                                             *
!                                                                      *
      Case('S12G','TS12G')
         Functional_type=GGA_type

!----    Swarts GGA exchange
!
!----    Perdew-Burk-Ernzerhof correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_S12G,XC_GGA_C_PBE]
!                                                                      *
!***********************************************************************
!                                                                      *
!     PBEsol                                                           *
!                                                                      *
      Case('PBESOL')
         Functional_type=GGA_type

!----    Perdew-Burk-Ernzerhof SOL exchange
!
!----    Perdew-Burk-Ernzerhof SOL correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_PBE_SOL,XC_GGA_C_PBE_SOL]
!                                                                      *
!***********************************************************************
!                                                                      *
!     RGE2                                                             *
!                                                                      *
      Case('RGE2')
         Functional_type=GGA_type

!----    Ruzsinzsky-Csonka-Scusceria exchange
!
!----    Perdew-Burk-Ernzerhof SOL correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_RGE2,XC_GGA_C_PBE_SOL]
!                                                                      *
!***********************************************************************
!                                                                      *
!     PTCA                                                             *
!                                                                      *
      Case('PTCA')
         Functional_type=GGA_type

!----    Perdew-Burk-Ernzerhof exchange
!
!----    Tognotti-Cortona-Adamo correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_PBE,XC_GGA_C_TCA]
!                                                                      *
!***********************************************************************
!                                                                      *
!     PBE0                                                             *
!                                                                      *
      Case('PBE0')
         Functional_type=GGA_type

!----    Perdew-Burk-Ernzerhof exchange
!
!----    Perdew-Burk-Ernzerhof correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_GGA_X_PBE,XC_GGA_C_PBE]
         Coeffs(1)=0.75D0
!        Coeffs(2)=1,00D0
!                                                                      *
!***********************************************************************
!                                                                      *
!     M06-L                                                            *
!                                                                      *
      Case('M06L')
         Functional_type=meta_GGA_type1

!----    Minnesota 2006 L exchange
!
!----    Minnesota 2006 L correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_MGGA_X_M06_L,XC_MGGA_C_M06_L]
!                                                                      *
!***********************************************************************
!                                                                      *
!     M06                                                              *
!                                                                      *
      Case('M06 ')
         Functional_type=meta_GGA_type1

!----    Minnesota 2006 exchange
!
!----    Minnesota 2006 correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_HYB_MGGA_X_M06,XC_MGGA_C_M06]
!                                                                      *
!***********************************************************************
!                                                                      *
!     M06-2X                                                           *
!                                                                      *
      Case('M062X')
         Functional_type=meta_GGA_type1

!----    Minnesota 2006 2X exchange
!
!----    Minnesota 2006 2X correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_HYB_MGGA_X_M06_2X,XC_MGGA_C_M06_2X]
!                                                                      *
!***********************************************************************
!                                                                      *
!     M06-HF                                                           *
!                                                                      *
      Case('M06HF')
         Functional_type=meta_GGA_type1

!----    Minnesota 2006 HF exchange
!
!----    Minnesota 2006 HF correlation

         nFuncs=2
         func_id(1:nFuncs)=[XC_HYB_MGGA_X_M06_HF,XC_MGGA_C_M06_HF]
!                                                                      *
!***********************************************************************
!                                                                      *
!      LDTF/LSDA (Thomas-Fermi for KE)                                 *
!                                                                      *
       Case('LDTF/LSDA ','LDTF/LDA  ')
         Functional_type=LDA_type

         If (Do_Core) Then
            nFuncs=1
            func_id(1:nFuncs)=[XC_LDA_C_VWN_RPA]
            Coeffs(1)=dFMD
         Else
            If (KEOnly) Then
               nFuncs=1
               func_id(1:nFuncs)=[XC_LDA_K_TF]
            Else
               nFuncs=3
               func_id(1:nFuncs)=[XC_LDA_K_TF,XC_LDA_X,XC_LDA_C_VWN_RPA]
            End If
         End If
!                                                                      *
!***********************************************************************
!                                                                      *
!      LDTF/LSDA5 (Thomas-Fermi for KE)                                *
!                                                                      *
       Case('LDTF/LSDA5','LDTF/LDA5 ')
         Functional_type=LDA_type

         If (Do_Core) Then
            nFuncs=1
            func_id(1:nFuncs)=[XC_LDA_C_VWN]
            Coeffs(1)=dFMD
         Else
              If (KEOnly) Then
                 nFuncs=1
                 func_id(1:nFuncs)=[XC_LDA_K_TF]
              Else
                 nFuncs=3
                 func_id(1:nFuncs)=[XC_LDA_K_TF,XC_LDA_X,XC_LDA_C_VWN]
             End If
         End If
!                                                                      *
!***********************************************************************
!                                                                      *
!      LDTF/PBE   (Thomas-Fermi for KE)                                *
!                                                                      *
       Case('LDTF/PBE  ')
         Functional_type=GGA_type

         If (Do_Core) Then
            nFuncs=1
            func_id(1:nFuncs)=[XC_GGA_C_PBE]
            Coeffs(1)=dFMD
         Else
            If (KEOnly) Then
               nFuncs=1
               func_id(1:nFuncs)=[XC_LDA_K_TF]
            Else
               nFuncs=3
               func_id(1:nFuncs)=[XC_LDA_K_TF,XC_GGA_X_PBE,XC_GGA_C_PBE]
            End If
         End If
!                                                                      *
!***********************************************************************
!                                                                      *
!      NDSD/PBE   (NDSD for KE)                                        *
!                                                                      *
       Case('NDSD/PBE  ')
         Functional_type=meta_GGA_type2

         If (Do_Core) Then
            nFuncs=1
            func_id(1:nFuncs)=[XC_GGA_C_PBE]
            Coeffs(1)=dFMD
         Else
            If (KEOnly) Then
               Sub => ndsd_ts
            Else
               Only_exc=.True.
               nFuncs=2
               func_id(1:nFuncs)=[XC_GGA_X_PBE,XC_GGA_C_PBE]
               External_Sub => ndsd_ts
            End If
         End If
!                                                                      *
!***********************************************************************
!                                                                      *
!      LDTF/BLYP  (Thomas-Fermi for KE)                                *
!                                                                      *
       Case('LDTF/BLYP ')
         Functional_type=GGA_type

         If (Do_Core) Then
            nFuncs=1
            func_id(1:nFuncs)=[XC_GGA_C_LYP]
            Coeffs(1)=dFMD
         Else
            If (KEOnly) Then
               nFuncs=1
               func_id(1:nFuncs)=[XC_LDA_K_TF]
            Else
               nFuncs=3
               func_id(1:nFuncs)=[XC_LDA_K_TF,XC_GGA_X_B88,XC_GGA_C_LYP]
            End If
         End If
!                                                                      *
!***********************************************************************
!                                                                      *
!      NDSD/BLYP  (NDSD for KE)                                        *
!                                                                      *
       Case('NDSD/BLYP ')
         Functional_type=meta_GGA_type2

         If (Do_Core) Then
            nFuncs=1
            func_id(1:nFuncs)=[XC_GGA_C_LYP]
            Coeffs(1)=dFMD
         Else
            If (KEOnly) Then
               Sub => ndsd_ts
            Else
               Only_exc=.True.
               nFuncs=2
               func_id(1:nFuncs)=[XC_GGA_X_B88,XC_GGA_C_LYP]
               External_Sub => ndsd_ts
            End If
         End If
!                                                                      *
!***********************************************************************
!                                                                      *
!      Kinetic only  (Thomas-Fermi)                                    *
!                                                                      *
       Case('TF_only')
         Functional_type=LDA_type

         nFuncs=1
         func_id(1:nFuncs)=[XC_LDA_K_TF]
!                                                                      *
!***********************************************************************
!                                                                      *
!      HUNTER  (von Weizsacker KE, no calc of potential)               *
!                                                                      *
       Case('HUNTER')
         Functional_type=GGA_type

         nFuncs=1
         func_id(1:nFuncs)=[XC_GGA_K_TFVW]
         Only_exc=.True.
!                                                                      *
!***********************************************************************
!                                                                      *
      Case default
         Call WarningMessage(2,' Driver: Undefined functional type!')
         Write (6,*) '         Functional=',KSDFT(1:LEN(KSDFT))
         Call Quit_OnUserError()
       End Select
!                                                                      *
!***********************************************************************
!                                                                      *
!     Now let's do some integration!
!     If the libxc interface is used do the proper initialization and closure.

      If (Associated(Sub,libxc_functionals)) Call Initiate_libxc_functionals(nD)

      Call DrvNQ(Sub,F_DFT,nD,Func,D_DS,nh1,nD,                         &
     &           Do_Grad,Grad,nGrad,Do_MO,Do_TwoEl,DFTFOCK)

      If (Associated(Sub,libxc_functionals)) Call Remove_libxc_functionals()

      If (Associated(External_Sub)) Call DrvNQ(External_Sub,F_DFT,nD,Func,D_DS,nh1,nD,         &
     &           Do_Grad,Grad,nGrad,Do_MO,Do_TwoEl,DFTFOCK)

      Sub          => Null()
      External_Sub => Null()
      Only_exc=.False.
!                                                                      *
!***********************************************************************
!                                                                      *
      End Subroutine Driver
