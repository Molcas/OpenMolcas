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
! Copyright (C) 2022, Ignacio Fdez. Galvan                             *
!***********************************************************************

module Functionals

use libxc_parameters, only: nFuncs_max
use nq_Info, only: Functional_Type, GGA_Type, LDA_Type, meta_GGA_Type1, meta_GGA_Type2, Other_Type
use Constants, only: Zero
use Definitions, only: wp, iwp, LibxcInt

implicit none
private

! Def_* variables are a cache of last read values, so the database file
! does not need to be re-read over and over

integer(kind=iwp) :: Def_Functional_Type = Other_Type, Def_nFuncs = 0
integer(kind=LibxcInt) :: Def_func_id(nFuncs_max) = -1_LibxcInt
real(kind=wp) :: Def_Coeffs(nFuncs_max) = Zero, Def_ExFac = Zero
character(len=80) :: Def_Label = ''
character(len=*), parameter :: Custom_File = 'CUSTFUNC', Custom_Func = '-999_CUSTOM_FUNCTIONAL'

public :: Custom_File, Custom_Func, Get_Func_ExFac, Get_Funcs, Init_Funcs, Print_Info

contains

subroutine Init_Funcs(Label)

  use fortran_strings, only: to_upper

  character(len=*), intent(in) :: Label
  character(len=len(Label)) :: UpLabel

  UpLabel = to_upper(Label)
  if (UpLabel /= Def_Label) call Find_Functional(UpLabel)

end subroutine Init_Funcs

function Get_Func_ExFac(Label)

  real(kind=wp) :: Get_Func_Exfac
  character(len=*), intent(in) :: Label

  call Init_Funcs(Label)
  Get_Func_ExFac = Def_ExFac

end function Get_Func_ExFac

subroutine Get_Funcs(Label)

  use libxc_parameters, only: Coeffs, func_id, nFuncs

  character(len=*), intent(in) :: Label

  call Init_Funcs(Label)
  nFuncs = Def_nFuncs
  Coeffs(1:Def_nFuncs) = Def_Coeffs(1:Def_nFuncs)
  func_id(1:Def_nFuncs) = Def_func_id(1:Def_nFuncs)
  Functional_Type = Def_Functional_Type

end subroutine Get_Funcs

subroutine Print_Info()

  use xc_f03_lib_m, only: xc_f03_func_end, xc_f03_func_get_info, xc_f03_func_info_get_name, xc_f03_func_info_get_references, &
                          xc_f03_func_info_t, xc_f03_func_init, xc_f03_func_reference_get_doi, xc_f03_func_reference_get_ref, &
                          xc_f03_func_reference_t, xc_f03_func_t, XC_UNPOLARIZED
  use Definitions, only: u6

  integer(kind=iwp) :: i, old_j
  integer(kind=LibxcInt) :: j
  type(xc_f03_func_t) :: func
  type(xc_f03_func_info_t) :: info
  type(xc_f03_func_reference_t) :: ref

  if (Def_nFuncs < 1) return

  write(u6,*)
  do i=1,Def_nFuncs
    call xc_f03_func_init(func,Def_func_id(i),xc_unpolarized)
    info = xc_f03_func_get_info(func)
    write(u6,100) trim(xc_f03_func_info_get_name(info))
    ! old_j is a workaround for a bug in Libxc 5.2.0
    old_j = -1
    j = 0
    do while ((j >= 0) .and. (j /= old_j))
      old_j = j
      ref = xc_f03_func_info_get_references(info,j)
      write(u6,101) trim(xc_f03_func_reference_get_ref(ref)),trim(xc_f03_func_reference_get_doi(ref))
    end do
    call xc_f03_func_end(func)
  end do

100 format(6x,'* ',a)
101 format(8x,'- ',a,' doi:',a)

end subroutine Print_Info

! Read functional database file and assign DFT parameters according to the chosen functional
subroutine Find_Functional(Label)

  use xc_f03_lib_m, only: xc_f03_func_end, xc_f03_func_get_info, xc_f03_func_info_get_family, xc_f03_func_info_get_flags, &
                          xc_f03_func_info_t, xc_f03_func_init, xc_f03_func_t, xc_f03_hyb_exx_coef, XC_FAMILY_GGA, &
                          XC_FAMILY_HYB_GGA, XC_FAMILY_HYB_LDA, XC_FAMILY_HYB_MGGA, XC_FAMILY_LDA, XC_FAMILY_MGGA, &
                          XC_FLAGS_NEEDS_LAPLACIAN, XC_UNPOLARIZED
  use fortran_strings, only: to_upper
  use Constants, only: One
  use Definitions, only: u6

  character(len=*), intent(in) :: Label
  integer(kind=iwp) :: i, istatus, Lu, nComp
  integer(kind=LibxcInt) :: flags(nFuncs_max)
  real(kind=wp) :: Coeff
  character(len=256) :: Line
  character(len=80) :: Labels(nFuncs_max), Word1, Word2
  integer(kind=iwp), external :: IsFreeUnit
  type(xc_f03_func_t) :: func(nFuncs_max)
  type(xc_f03_func_info_t) :: info(nFuncs_max)

  Def_ExFac = Zero
  Def_nFuncs = 0

  ! First test if this is already a Libxc functional
  Def_func_id(1) = get_func(Label,test=.true.)
  if (Def_func_id(1) >= 0) then

    Def_nFuncs = 1
    Def_Coeffs(1) = One
    Labels(1) = Label

  else

    ! If not, we have to read the database file, or the custom functional file
    Lu = IsFreeUnit(11)
    if (Label == Custom_Func) then
      call molcas_open(Lu,Custom_File)
    else
      call molcas_open(Lu,'FUNCDATA')
    end if

    ! Find the line that starts with the keyword name
    do
      read(Lu,'(A)',iostat=istatus) Line
      if (istatus /= 0) then
        call WarningMessage(2,' Find_Functional: Undefined functional type!')
        write(u6,*) '         Functional=',trim(Label)
        call Quit_OnUserError()
      end if
      Line = adjustl(Line)
      if ((Line == '') .or. (Line(1:1) == '#')) cycle
      read(Line,*) Word1,Word2
      if (to_upper(Word1) == Label) exit
    end do

    ! Once found, read the second word
    read(Word2,*,iostat=istatus) nComp
    if (istatus == 0) then
      ! If it's a number, read the component functionals and factors
      if (Def_nFuncs > nFuncs_max) then
        call WarningMessage(2,' Find_Functional: Too many components!')
        write(u6,*) '         nFuncs=',Def_nFuncs
        call Quit_OnUserError()
      end if
      i = 0
      do while (i < nComp)
        read(Lu,'(A)',iostat=istatus) Line
        if (istatus /= 0) then
          call WarningMessage(2,' Find_Functional: Error in functional definition!')
          write(u6,*) '         Functional=',trim(Label)
          call Quit_OnUserError()
        end if
        Line = adjustl(Line)
        if ((Line == '') .or. (Line(1:1) == '#')) cycle
        i = i+1
        read(Line,*,iostat=istatus) Coeff,Word2
        if (istatus /= 0) then
          call WarningMessage(2,' Find_Functional: Error in functional definition!')
          write(u6,*) '         Functional=',trim(Label)
          call Quit_OnUserError()
        end if
        ! HF_X means exact exchange
        if (to_upper(Word2) == 'HF_X') then
          Def_ExFac = Def_ExFac+Coeff
        else
          Def_nFuncs = Def_nFuncs+1
          Def_Coeffs(Def_nFuncs) = Coeff
          Def_func_id(Def_nFuncs) = get_func(Word2)
          Labels(Def_nFuncs) = Word2
        end if
      end do
    else
      ! Otherwise, this is just an alias for a Libxc functional
      Def_nFuncs = 1
      Def_Coeffs(1) = One
      Def_func_id(1) = get_func(Word2)
      Labels(1) = Word2
    end if

    close(Lu)

  end if

  ! Now process the functional(s)

  do i=1,Def_nFuncs
    call xc_f03_func_init(func(i),Def_func_id(i),XC_UNPOLARIZED)
    info(i) = xc_f03_func_get_info(func(i))
    flags(i) = xc_f03_func_info_get_flags(info(i))
  end do

  Def_Functional_Type = Other_Type
  do i=1,Def_nFuncs
    ! Check whether the functional uses some unsupported feature
    call check_supported(Labels(i),flags(i))
    ! Add exact exchange from components
    Def_ExFac = Def_ExFac+Def_Coeffs(i)*xc_f03_hyb_exx_coef(func(i))
    ! Assign functional type ("maximum" of its components)
    select case (xc_f03_func_info_get_family(info(i)))
      case (XC_FAMILY_LDA,XC_FAMILY_HYB_LDA)
        Def_Functional_type = max(Def_Functional_type,LDA_type)
      case (XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
        Def_Functional_type = max(Def_Functional_type,GGA_type)
      case (XC_FAMILY_MGGA,XC_FAMILY_HYB_MGGA)
        if (iand(flags(i),XC_FLAGS_NEEDS_LAPLACIAN) > 0) then
          Def_Functional_type = max(Def_Functional_type,meta_GGA_type2)
        else
          Def_Functional_type = max(Def_Functional_type,meta_GGA_type1)
        end if
    end select
  end do

  do i=1,Def_nFuncs
    call xc_f03_func_end(func(i))
  end do

  Def_Label = Label

end subroutine Find_Functional

function get_func(xcLabel,test)

  use xc_f03_lib_m, only: xc_f03_functional_get_number
  use Definitions, only: u6

  integer(kind=LibxcInt) :: get_func
  character(len=*), intent(in) :: xcLabel
  logical(kind=iwp), intent(in), optional :: test
  logical(kind=iwp) :: do_test

  if (present(test)) then
    do_test = test
  else
    do_test = .false.
  end if

  get_func = xc_f03_functional_get_number(xcLabel)
  if ((get_func < 0) .and. (.not. do_test)) then
    call WarningMessage(2,' Find_Functional: Undefined functional in Libxc!')
    write(u6,*) '         Functional=',trim(xcLabel)
    call Quit_OnUserError()
  end if

end function get_func

subroutine check_supported(Label,flags)

  use xc_f03_lib_m, only: XC_FLAGS_HYB_CAM, XC_FLAGS_HYB_CAMY, XC_FLAGS_HYB_LC, XC_FLAGS_HYB_LCY, XC_FLAGS_VV10
  use Definitions, only: u6

  character(len=*), intent(in) :: Label
  integer(kind=LibxcInt), intent(in) :: flags
  integer(kind=iwp) :: lev, lt, maxlev

  maxlev = 0
  lt = len_trim(Label)

  if ((iand(flags,XC_FLAGS_HYB_CAM) > 0) .or. (iand(flags,XC_FLAGS_HYB_CAMY) > 0) .or. (iand(flags,XC_FLAGS_HYB_LC) > 0) .or. &
      (iand(flags,XC_FLAGS_HYB_LCY) > 0)) then
    lev = 2
    maxlev = max(lev,maxlev)
    call WarningMessage(lev,' Find_Functional: Range separation is not supported!')
  end if
  if ((iand(flags,XC_FLAGS_VV10) > 0)) then
    lev = 2
    maxlev = max(lev,maxlev)
    call WarningMessage(lev,' Find_Functional: Non-local correlation is not supported!')
  end if
  if ((Label(lt-1:lt) == '_D') .or. (Label(lt-2:lt) == '_D3')) then
    lev = 1 ! make it just a warning for now
    maxlev = max(lev,maxlev)
    call WarningMessage(lev,' Find_Functional: Dispersion corrections are not implemented!')
  end if
  if (maxlev > 0) write(u6,*) '         Functional=',trim(Label)
  if (maxlev > 1) call Quit_OnUserError()

end subroutine check_supported

end module Functionals
