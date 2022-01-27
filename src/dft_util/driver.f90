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
Subroutine Driver(KSDFA,Do_Grad,Func,Grad,nGrad,Do_MO,Do_TwoEl,D_DS,F_DFT,nh1,nD,DFTFOCK)

use libxc_parameters
use OFembed, only: KEOnly, dFMD, Do_Core
use libxc,   only: Only_exc
use nq_Grid, only: l_casdft
use nq_pdft, only: lft
Implicit None
#include "real.fh"
#include "nq_info.fh"
Character*(*) KSDFA
Logical Do_Grad
Integer :: nGrad, nh1, nD
Real*8 :: ExFac, Func, Grad(nGrad)
Logical Do_MO, Do_TwoEl
Real*8 :: D_DS(nh1,nD), F_DFT(nh1,nD)
Character*4 DFTFOCK
logical :: LDTF=.False., NDSD=.False.
character(LEN=12) :: FLabel=''

abstract interface
   Subroutine DFT_FUNCTIONAL(mGrid,nD)
      Integer mGrid, nD
   end subroutine
end interface

!***********************************************************************
!     Define external functions not defined in LibXC. These are either
!     accessed through the procedure pointer sub or External_sub.

procedure(DFT_FUNCTIONAL) :: Overlap, NucAtt, ndsd_ts
!***********************************************************************
procedure(DFT_FUNCTIONAL), pointer :: sub => null()
!     Sometime we need an external routine which covers something which
!     Libxc doesn't support.
procedure(DFT_FUNCTIONAL), pointer :: External_sub => null()
!                                                                      *
!***********************************************************************
! Global variable for MCPDFT functionals                               *
FLabel=KSDFA ! The user could be passing an explicit string! Hence, the local copy.

!
!     Set some flags and clean up the label to be just the label of the
!     underlaying DFT functional.
!
l_casdft = FLabel(1:2).eq.'T:' .or. FLabel(1:3).eq.'FT:'

lft      = FLabel(1:3).eq.'FT:'

If (l_casdft) Then
   If (lft) Then
      FLabel=FLabel(4:)
   Else
      FLabel=FLabel(3:)
   End If
   Do_MO=.true.
   Do_TwoEl=.true.
   If (.NOT.Do_PDFTPOT .and. .Not.DO_Grad) Only_exc=.True.
End If

If (FLabel(1:5)=='LDTF/')Then
   LDTF=.true.
   FLabel=FLabel(6:)
End If
If (FLabel(1:5)=='NDSD/')Then
   NDSD=.true.
   FLabel=FLabel(6:)
End If
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
    Select Case(FLabel)
!                                                                      *
!***********************************************************************
!                                                                      *
!      LSDA LDA SVWN                                                   *
!                                                                      *
      Case('LSDA ','LDA ','SVWN ')
         Functional_type=LDA_type

!----    Slater exchange
!
!----    Vosko-Wilk-Nusair correlation functional III

         If (Do_Core) Then
            nFuncs=1
            func_id(1:nFuncs)=[XC_LDA_C_VWN_RPA]
            Coeffs(1)=dFMD
         Else
            If (LDTF) Then
               If (KEOnly) Then
                  nFuncs=1
                  func_id(1:nFuncs)=[XC_LDA_K_TF]
               Else
                  nFuncs=3
                  func_id(1:nFuncs)=[XC_LDA_K_TF,XC_LDA_X,XC_LDA_C_VWN_RPA]
               End If
            Else
               nFuncs=2
               func_id(1:nFuncs)=[XC_LDA_X,XC_LDA_C_VWN_RPA]
            End If
         End If
!                                                                      *
!***********************************************************************
!                                                                      *
!      LSDA5 LDA5 SVWN5                                                *
!                                                                      *
      Case('LSDA5','LDA5','SVWN5')
         Functional_type=LDA_type

!----    Slater exchange
!
!----    Vosko-Wilk-Nusair correlation functional V

         If (Do_Core) Then
            nFuncs=1
            func_id(1:nFuncs)=[XC_LDA_C_VWN]
            Coeffs(1)=dFMD
         Else
            If (LDTF) Then
               If (KEOnly) Then
                  nFuncs=1
                  func_id(1:nFuncs)=[XC_LDA_K_TF]
               Else
                  nFuncs=3
                  func_id(1:nFuncs)=[XC_LDA_K_TF,XC_LDA_X,XC_LDA_C_VWN]
              End If
            Else
               nFuncs=2
               func_id(1:nFuncs)=[XC_LDA_X,XC_LDA_C_VWN]
            End If
         End If
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
!     BLYP                                                             *
!                                                                      *
      Case('BLYP')
         Functional_type=GGA_type

!----    Slater exchange + B88 exchange
!
!----    Lee-Yang-Parr correlation

         If (Do_Core) Then
            nFuncs=1
            func_id(1:nFuncs)=[XC_GGA_C_LYP]
            Coeffs(1)=dFMD
         Else
            If (LDTF) Then
               If (KEOnly) Then
                  nFuncs=1
                  func_id(1:nFuncs)=[XC_LDA_K_TF]
               Else
                  nFuncs=3
                  func_id(1:nFuncs)=[XC_LDA_K_TF,XC_GGA_X_B88,XC_GGA_C_LYP]
               End If
            Else If (NDSD) Then
               If (KEOnly) Then
                  Sub => ndsd_ts
               Else
                  Only_exc=.True.
                  nFuncs=2
                  func_id(1:nFuncs)=[XC_GGA_X_B88,XC_GGA_C_LYP]
                  External_Sub => ndsd_ts
               End If
            Else
               nFuncs=2
               func_id(1:nFuncs)=[XC_GGA_X_B88,XC_GGA_C_LYP]
            End If
         End If
!                                                                      *
!***********************************************************************
!                                                                      *
!     PBE                                                              *
!                                                                      *
      Case('PBE ')
         Functional_type=GGA_type

!----    Perdew-Burk-Ernzerhof exchange
!
!----    Perdew-Burk-Ernzerhof correlation

         If (Do_Core) Then
            nFuncs=1
            func_id(1:nFuncs)=[XC_GGA_C_PBE]
            Coeffs(1)=dFMD
         Else
            If (LDTF) Then
               If (KEOnly) Then
                  nFuncs=1
                  func_id(1:nFuncs)=[XC_LDA_K_TF]
               Else
                  nFuncs=3
                  func_id(1:nFuncs)=[XC_LDA_K_TF,XC_GGA_X_PBE,XC_GGA_C_PBE]
               End If
            Else If (NDSD) Then
               If (KEOnly) Then
                  Sub => ndsd_ts
               Else
                  Only_exc=.True.
                  nFuncs=2
                  func_id(1:nFuncs)=[XC_GGA_X_PBE,XC_GGA_C_PBE]
                  External_Sub => ndsd_ts
               End If
            Else
               nFuncs=2
               func_id(1:nFuncs)=[XC_GGA_X_PBE,XC_GGA_C_PBE]
            End If
         End If
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
         Call Find_Functional(FLabel,ExFac)

       End Select
!                                                                      *
!***********************************************************************
!                                                                      *
       If (Functional_type/=LDA_type.and.Functional_type/=GGA_type.and.l_CasDFT) Then
          Write (6,*) ' MC-PDFT combined with invalid functional class'
          Call Abend()
       End If
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
      LDTF=.False.
      NDSD=.False.
!                                                                      *
!***********************************************************************
!                                                                      *
      End Subroutine Driver
