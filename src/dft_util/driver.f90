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
use xc_f03_lib_m
use Functionals, only: Get_Funcs
use OFembed, only: KEOnly, dFMD, Do_Core
use libxc,   only: Only_exc
use nq_Grid, only: l_casdft
use nq_pdft, only: lft
use Definitions, only: LibxcInt
use nq_Info
Implicit None
#include "real.fh"
#include "ksdft.fh"
Character*(*) KSDFA
Logical Do_Grad
Integer :: i, j, nGrad, nh1, nD
Real*8 :: Func, Grad(nGrad)
Logical Do_MO, Do_TwoEl
Real*8 :: D_DS(nh1,nD), F_DFT(nh1,nD)
Character*4 DFTFOCK
logical :: LDTF, NDSD
character(LEN=80) :: FLabel
type(xc_f03_func_t) :: func_
type(xc_f03_func_info_t) :: info_

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
   FLabel=FLabel(Index(FLabel,'T:')+2:)
   Do_MO=.true.
   Do_TwoEl=.true.
   If (.NOT.Do_PDFTPOT .and. .Not.DO_Grad) Only_exc=.True.
End If

If (FLabel(1:5)=='LDTF/')Then
   LDTF=.true.
   FLabel=FLabel(6:)
Else
   LDTF=.false.
End If
If (FLabel(1:5)=='NDSD/')Then
   NDSD=.true.
   FLabel=FLabel(6:)
Else
   NDSD=.false.
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
! The names TF_only and HUNTER are hardcoded in some parts of the code,*
! so we define them explicitly instead of relying on the external file *
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
         Call Get_Funcs(FLabel)

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
      If (Do_Core) Then
         ! Keep only correlation
         Do i=1,nFuncs
            Call xc_f03_func_init(func_,func_id(i),0_LibxcInt)
            info_ = xc_f03_func_get_info(func_)
            If (xc_f03_func_info_get_kind(info_) == XC_CORRELATION) Then
               Coeffs(i) = Coeffs(i)*dFMD
            Else
               Coeffs(i) = Zero
            End If
            Call xc_f03_func_end(func_)
         End Do
      Else If (LDTF) Then
         ! Add TF kinetic with same coeff as exchange
         ! and optionally kill everything else
         Do i=1,nFuncs
            Call xc_f03_func_init(func_,func_id(i),0_LibxcInt)
            info_ = xc_f03_func_get_info(func_)
            If (xc_f03_func_info_get_kind(info_) == XC_EXCHANGE) Then
               If (nFuncs == nFuncs_max) Then
                  Write (6,*) ' Too many functionals for LDTF'
                  Call Abend()
               End If
               func_id(nFuncs+1) = XC_LDA_K_TF
               Coeffs(nFuncs+1) = Coeffs(i)
               nFuncs = nFuncs+1
            End If
            If (KEOnly) Coeffs(i) = Zero
            Call xc_f03_func_end(func_)
         End Do
      Else If (NDSD) Then
         ! Add ndsd_ts, and optionally kill everything else
         If (KEOnly) Then
            Coeffs(:) = Zero
            Sub => ndsd_ts
         Else
            Only_exc=.True.
            External_Sub => ndsd_ts
         End If
      End If
      ! Reduce list
      j = 0
      Do i=1,nFuncs
        If (Coeffs(i) == Zero) Cycle
        j = j+1
        If (j == i) Cycle
        Coeffs(j) = Coeffs(i)
        func_id(j) = func_id(i)
      End Do
      nFuncs = j
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
