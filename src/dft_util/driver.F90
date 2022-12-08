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

subroutine Driver(KSDFA,Do_Grad,Func,Grad,nGrad,Do_MO,Do_TwoEl,D_DS,F_DFT,nh1,nD,DFTFOCK)

use libxc_parameters, only: Coeffs, func_id, initiate_libxc_functionals, libxc_functionals, nFuncs, nFuncs_max, &
                            remove_libxc_functionals
use xc_f03_lib_m, only: XC_CORRELATION, XC_EXCHANGE, xc_f03_func_end, xc_f03_func_get_info, xc_f03_func_info_get_kind, &
                        xc_f03_func_init, xc_f03_func_t, xc_f03_func_info_t, XC_GGA_K_TFVW, XC_LDA_K_TF, XC_UNPOLARIZED
use Functionals, only: Get_Funcs
use KSDFT_Info, only: Do_PDFTPOT
use OFembed, only: dFMD, Do_Core, KEOnly
use libxc, only: Only_exc
use nq_Grid, only: l_casdft
use nq_Info, only: Functional_type, GGA_Type, LDA_Type
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(in) :: KSDFA
logical(kind=iwp), intent(in) :: Do_Grad
integer(kind=iwp), intent(in) :: nGrad, nh1, nD
real(kind=wp), intent(inout) :: Func, Grad(nGrad), F_DFT(nh1,nD)
logical(kind=iwp), intent(inout) :: Do_MO, Do_TwoEl
real(kind=wp), intent(in) :: D_DS(nh1,nD)
character(len=4), intent(in) :: DFTFOCK
integer(kind=iwp) :: i, j
logical(kind=iwp) :: IsFT, LDTF, NDSD
character(len=80) :: FLabel
type(xc_f03_func_t) :: func_
type(xc_f03_func_info_t) :: info_
!***********************************************************************
! Define external functions not defined in LibXC. These are either
! accessed through the procedure pointer Sub or External_sub.
abstract interface
  subroutine DFT_FUNCTIONAL(mGrid,nD)
    import :: iwp
    integer(kind=iwp), intent(in) :: mGrid, nD
  end subroutine DFT_FUNCTIONAL
end interface
procedure(DFT_FUNCTIONAL) :: Overlap, NucAtt, ndsd_ts
procedure(DFT_FUNCTIONAL), pointer :: Sub, External_sub
!***********************************************************************

!                                                                      *
!***********************************************************************
! Global variable for MCPDFT functionals                               *
FLabel = KSDFA ! The user could be passing an explicit string! Hence, the local copy.

! Set some flags and clean up the label to be just the label of the
! underlaying DFT functional.

l_casdft = (FLabel(1:2) == 'T:') .or. (FLabel(1:3) == 'FT:')

IsFT = FLabel(1:3) == 'FT:'

if (l_casdft) then
  FLabel = FLabel(index(FLabel,'T:')+2:)
  Do_MO = .true.
  Do_TwoEl = .true.
  if ((.not. Do_PDFTPOT) .and. (.not. DO_Grad)) Only_exc = .true.
end if

if (FLabel(1:5) == 'LDTF/') then
  LDTF = .true.
  FLabel = FLabel(6:)
else
  LDTF = .false.
end if
if (FLabel(1:5) == 'NDSD/') then
  NDSD = .true.
  FLabel = FLabel(6:)
else
  NDSD = .false.
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Default is to use the libxc interface
! Coefficient for the individual contibutions are defaulted to 1.0

Sub => libxc_functionals     ! Default
External_Sub => null()       ! Default
Coeffs(:) = One              ! Default
!                                                                      *
!***********************************************************************
!                                                                      *
select case (FLabel)
!                                                                      *
!***********************************************************************
!                                                                      *
! Overlap

  case ('Overlap')
    Functional_type = LDA_type
    Sub => Overlap
!                                                                      *
!***********************************************************************
!                                                                      *
! NucAtt

  case ('NucAtt')
    Functional_type = LDA_type
    Sub => NucAtt
!                                                                      *
!***********************************************************************
!                                                                      *
! The names TF_only and HUNTER are hardcoded in some parts of the code,*
! so we define them explicitly instead of relying on the external file *
!                                                                      *
!***********************************************************************
!                                                                      *
! Kinetic only (Thomas-Fermi)

  case ('TF_only')
    Functional_type = LDA_type

    nFuncs = 1
    func_id(1:nFuncs) = [XC_LDA_K_TF]
!                                                                      *
!***********************************************************************
!                                                                      *
!  HUNTER (von Weizsacker KE, no calc of potential)

  case ('HUNTER')
    Functional_type = GGA_type

    nFuncs = 1
    func_id(1:nFuncs) = [XC_GGA_K_TFVW]
    Only_exc = .true.
!                                                                      *
!***********************************************************************
!                                                                      *
  case default
    call Get_Funcs(FLabel)

end select
!                                                                      *
!***********************************************************************
!                                                                      *
if ((Functional_type /= LDA_type) .and. (Functional_type /= GGA_type) .and. l_CasDFT) then
  write(u6,*) ' MC-PDFT combined with invalid functional class'
  call Abend()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (Do_Core) then
  ! Keep only correlation
  do i=1,nFuncs
    call xc_f03_func_init(func_,func_id(i),XC_UNPOLARIZED)
    info_ = xc_f03_func_get_info(func_)
    if (xc_f03_func_info_get_kind(info_) == XC_CORRELATION) then
      Coeffs(i) = Coeffs(i)*dFMD
    else
      Coeffs(i) = Zero
    end if
    call xc_f03_func_end(func_)
  end do
else if (LDTF) then
  ! Add TF kinetic with same coeff as exchange
  ! and optionally kill everything else
  do i=1,nFuncs
    call xc_f03_func_init(func_,func_id(i),XC_UNPOLARIZED)
    info_ = xc_f03_func_get_info(func_)
    if (xc_f03_func_info_get_kind(info_) == XC_EXCHANGE) then
      if (nFuncs == nFuncs_max) then
        write(u6,*) ' Too many functionals for LDTF'
        call Abend()
      end if
      func_id(nFuncs+1) = XC_LDA_K_TF
      Coeffs(nFuncs+1) = Coeffs(i)
      nFuncs = nFuncs+1
    end if
    if (KEOnly) Coeffs(i) = Zero
    call xc_f03_func_end(func_)
  end do
else if (NDSD) then
  ! Add ndsd_ts, and optionally kill everything else
  if (KEOnly) then
    Coeffs(:) = Zero
    Sub => ndsd_ts
  else
    Only_exc = .true.
    External_Sub => ndsd_ts
  end if
end if
! Reduce list
j = 0
do i=1,nFuncs
  if (Coeffs(i) == Zero) cycle
  j = j+1
  if (j == i) cycle
  Coeffs(j) = Coeffs(i)
  func_id(j) = func_id(i)
end do
nFuncs = j
!                                                                      *
!***********************************************************************
!                                                                      *
! Now let's do some integration!
! If the libxc interface is used do the proper initialization and closure.

if (associated(Sub,libxc_functionals)) call Initiate_libxc_functionals(nD)

call DrvNQ(Sub,F_DFT,nD,Func,D_DS,nh1,nD,Do_Grad,Grad,nGrad,Do_MO,Do_TwoEl,DFTFOCK,IsFT)

if (associated(Sub,libxc_functionals)) call Remove_libxc_functionals()

if (associated(External_Sub)) call DrvNQ(External_Sub,F_DFT,nD,Func,D_DS,nh1,nD,Do_Grad,Grad,nGrad,Do_MO,Do_TwoEl,DFTFOCK,IsFT)

Only_exc = .false.
LDTF = .false.
NDSD = .false.
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Driver
