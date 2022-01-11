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
      Implicit None
#include "nq_info.fh"
      Character*(*) KSDFT
      Logical Do_Grad
      Integer :: nGrad, nh1, nD
      Real*8 :: Func, Grad(nGrad)
      Logical Do_MO, Do_TwoEl
      Real*8 :: D_DS(nh1,nD), F_DFT(nh1,nD)
      Character*4 DFTFOCK

      External LSDA, Overlap, BLYP, BPBE, B3LYP, HFS, HFB,XAlpha, LSDA5, B3LYP5, B2PLYP, TLYP
      External NucAtt, OLYP, O3LYP, OPBE,PBE, PBE0, PBEsol, M06L, M06, M062X, HFO
      External M06HF, SSBSW, SSBD, HFG, GLYP, GPBE, HFB86, B86LYP, B86PBE, BWIG, KT3
      External O2PLYP, KT2, RGE2, REVPBE, PTCA, S12G, S12H

      abstract interface
          Subroutine DFT_FUNCTIONAL(mGrid,nD)
          Integer mGrid, nD
          end subroutine
      end interface

      procedure(DFT_FUNCTIONAL), pointer :: sub => null()

!                                                                      *
!***********************************************************************
!                                                                      *

!      Default is to use the libxc interface
!      Coefficient for the individual contibutions are defaulted to 1.0D0

       Sub => libxc_functionals     ! Default
       Coeffs(:)=1.0D0              ! Default
       Select Case(KSDFT)
!                                                                      *
!***********************************************************************
!                                                                      *
!      LSDA LDA SVWN                                                   *
!                                                                      *
      Case('LSDA ','LDA ','TLSDA','FTLSDA','SVWN ')
         Functional_type=LDA_type
         If(KSDFT.eq.'TLSDA'.or.KSDFT.eq.'FTLSDA') Do_MO=.true.
         If(KSDFT.eq.'TLSDA'.or.KSDFT.eq.'FTLSDA') Do_TwoEl=.true.

!----    Slater exchange
!
!----    Vosko-Wilk-Nusair correlation functional III

         nFuncs=2
         func_id(1:nFuncs)=[int(1,4),int(8,4)]
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

         nFuncs=2
         func_id(1:nFuncs)=[int(1,4),int(7,4)]
!                                                                      *
!***********************************************************************
!                                                                      *
!     HFB                                                              *
!                                                                      *
       Case('HFB')
         Functional_type=GGA_type

!----    Slate exchange + Becke 88 exchange

         nFuncs=1
         func_id(1:nFuncs)=[int(106,4)]
!                                                                      *
!***********************************************************************
!                                                                      *
!     HFO                                                              *
!                                                                      *
       Case('HFO')
         Functional_type=GGA_type

!----    Slate exchange + OPTx exchange

         nFuncs=1
         func_id(1:nFuncs)=[int(110,4)]
!                                                                      *
!***********************************************************************
!                                                                      *
!     HFG                                                              *
!                                                                      *
       Case('HFG')
         Functional_type=GGA_type

!----    Slate exchange + G96 exchange

         nFuncs=1
         func_id(1:nFuncs)=[int(107,4)]
!                                                                      *
!***********************************************************************
!                                                                      *
!     HFB86                                                            *
!                                                                      *
       Case('HFB86')
         Functional_type=GGA_type

!----    Slate exchange + B86 exchange

         nFuncs=1
         func_id(1:nFuncs)=[int(103,4)]
!                                                                      *
!***********************************************************************
!                                                                      *
!      HFS                                                             *
!                                                                      *
       Case('HFS')
         Functional_type=LDA_type

!----    Slate exchange

         nFuncs=1
         func_id(1:nFuncs)=[int(1,4)]
!                                                                      *
!***********************************************************************
!                                                                      *
!      XALPHA                                                          *
!                                                                      *
       Case('XALPHA')
         Functional_type=LDA_type

!----    Slate exchange

         nFuncs=1
         func_id(1:nFuncs)=[int(1,4)]
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
         func_id(1:nFuncs)=[int(106,4),int(573,4)]
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

         nFuncs=2
         func_id(1:nFuncs)=[int(106,4),int(131,4)]

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
         func_id(1:nFuncs)=[int(110,4),int(131,4)]
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
         func_id(1:nFuncs)=[int(1,4),int(110,4),int(145,4),int(131,4)]
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
         func_id(1:nFuncs)=[int(1,4),int(145,4),int(8,4)]
         Coeffs(1)= 0.07173d0
!        Coeffs(2)= One
         Coeffs(3)= 0.576727d0
!                                                                      *
!***********************************************************************
!                                                                      *
!     GLYP                                                             *
!                                                                      *
      Case('GLYP')
         Functional_type=GGA_type
         Sub => GLYP
!                                                                      *
!***********************************************************************
!                                                                      *
!     B86LYP                                                           *
!                                                                      *
      Case('B86LYP')
         Functional_type=GGA_type
         Sub => B86LYP
!                                                                      *
!***********************************************************************
!                                                                      *
!     BPBE                                                             *
!                                                                      *
      Case('BPBE')
         Functional_type=GGA_type
         Sub => BPBE
!                                                                      *
!***********************************************************************
!                                                                      *
!     OPBE                                                             *
!                                                                      *
      Case('OPBE','TOPBE','FTOPBE')
         Functional_type=GGA_type
         Sub => OPBE
         If(KSDFT.eq.'TOPBE'.or.KSDFT.eq.'FTOPBE') Do_MO=.true.
         If(KSDFT.eq.'TOPBE'.or.KSDFT.eq.'FTOPBE') Do_TwoEl=.true.
!                                                                      *
!***********************************************************************
!                                                                      *
!     GPBE                                                             *
!                                                                      *
      Case('GPBE')
         Functional_type=GGA_type
         Sub => GPBE
!                                                                      *
!***********************************************************************
!                                                                      *
!     B86PBE                                                           *
!                                                                      *
      Case('B86PBE')
         Functional_type=GGA_type
         Sub => B86PBE
!                                                                      *
!***********************************************************************
!                                                                      *
!     TLYP                                                             *
!                                                                      *
      Case('TLYP')
         Functional_type=GGA_type
         Sub => TLYP
!                                                                      *
!***********************************************************************
!                                                                      *
!     B3LYP                                                            *
!                                                                      *
      Case('B3LYP ')
         Functional_type=GGA_type
         Sub => B3LYP
!                                                                      *
!***********************************************************************
!                                                                      *
!     O3LYP                                                            *
!                                                                      *
      Case('O3LYP ')
         Functional_type=GGA_type
         Sub => O3LYP
!                                                                      *
!***********************************************************************
!                                                                      *
!     B2PLYP                                                           *
!                                                                      *
      Case('B2PLYP')
         Functional_type=GGA_type
         Sub => B2PLYP
!                                                                      *
!***********************************************************************
!                                                                      *
!     O2PLYP                                                           *
!                                                                      *
      Case('O2PLYP')
         Functional_type=GGA_type
         Sub => O2PLYP
!                                                                      *
!***********************************************************************
!                                                                      *
!     B3LYP5                                                           *
!                                                                      *
      Case('B3LYP5')
         Functional_type=GGA_type
         Sub => B3LYP5
!                                                                      *
!***********************************************************************
!                                                                      *
!     PBE                                                              *
!                                                                      *
      Case('PBE','TPBE','FTPBE')
         Functional_type=GGA_type
         Sub => PBE
         If(KSDFT.eq.'TPBE'.or.KSDFT.eq.'FTPBE') Do_MO=.true.
         If(KSDFT.eq.'TPBE'.or.KSDFT.eq.'FTPBE') Do_TwoEl=.true.
!                                                                      *
!***********************************************************************
!                                                                      *
!     revPBE                                                           *
!                                                                      *
      Case('REVPBE','TREVPBE','FTREVPBE')
         Functional_type=GGA_type
         Sub => REVPBE
         If(KSDFT.eq.'TREVPBE'.or.KSDFT.eq.'FTREVBPE') Do_MO=.true.
         If(KSDFT.eq.'TREVPBE'.or.KSDFT.eq.'FTREVPBE') Do_TwoEl=.true.
!                                                                      *
!***********************************************************************
!                                                                      *
!     SSBSW                                                              *
!                                                                      *
      Case('SSBSW','TSSBSW')
         Functional_type=GGA_type
         Sub => SSBSW
!                                                                      *
!***********************************************************************
!                                                                      *
!     SSBD                                                             *
!                                                                      *
      Case('SSBD','TSSBD')
         Functional_type=GGA_type
         Sub => SSBD
!                                                                      *
!                                                                      *
!***********************************************************************
!                                                                      *
!     S12H                                                             *
!                                                                      *
      Case('S12H')
         Functional_type=GGA_type
         Sub => S12H
!                                                                      *
!***********************************************************************
!                                                                      *
!     S12G                                                             *
!                                                                      *
      Case('S12G','TS12G')
         Functional_type=GGA_type
         Sub => S12G
!                                                                      *
!***********************************************************************
!                                                                      *
!     PBEsol                                                           *
!                                                                      *
      Case('PBESOL')
         Functional_type=GGA_type
         Sub => PBESol
!                                                                      *
!***********************************************************************
!                                                                      *
!     RGE2                                                             *
!                                                                      *
      Case('RGE2')
         Functional_type=GGA_type
         Sub => RGE2
!                                                                      *
!***********************************************************************
!                                                                      *
!     PTCA                                                             *
!                                                                      *
      Case('PTCA')
         Functional_type=GGA_type
         Sub => PTCA
!                                                                      *
!***********************************************************************
!                                                                      *
!     PBE0                                                             *
!                                                                      *
      Case('PBE0')
         Functional_type=GGA_type
         Sub => PBE0
!                                                                      *
!***********************************************************************
!                                                                      *
!     M06-L                                                            *
!                                                                      *
      Case('M06L')
         Functional_type=meta_GGA_type1
         Sub => M06L
!                                                                      *
!***********************************************************************
!                                                                      *
!     M06                                                              *
!                                                                      *
      Case('M06 ')
         Functional_type=meta_GGA_type1
         Sub => M06
!                                                                      *
!***********************************************************************
!                                                                      *
!     M06-2X                                                           *
!                                                                      *
      Case('M062X')
         Functional_type=meta_GGA_type1
         Sub => M062X
!                                                                      *
!***********************************************************************
!                                                                      *
!     M06-HF                                                           *
!                                                                      *
      Case('M06HF')
         Functional_type=meta_GGA_type1
         Sub => M06HF
!                                                                      *
!***********************************************************************
!                                                                      *
      Case default
         Call WarningMessage(2,' DrvDFT: Undefined functional type!')
         Write (6,*) '         Functional=',KSDFT(1:LEN(KSDFT))
         Call Quit_OnUserError()
       End Select
!                                                                      *
!***********************************************************************
!                                                                      *
      If (Associated(Sub,libxc_functionals)) Call Initiate_libxc_functionals(nD)

      Call DrvNQ(Sub,F_DFT,nD,Func,                                     &
     &           D_DS,nh1,nD,                                           &
     &           Do_Grad,                                               &
     &           Grad,nGrad,                                            &
     &           Do_MO,Do_TwoEl,DFTFOCK)

      If (Associated(Sub,libxc_functionals)) Call Remove_libxc_functionals()

      Sub => Null()
!                                                                      *
!***********************************************************************
!                                                                      *
      End Subroutine Driver
