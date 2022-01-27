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
!
! Read functional database file and assign DFT parameters according
! to the chosen functional

subroutine Find_Functional(Label,ExFac)

use libxc_parameters, only: Coeffs, func_id, nFuncs, nFuncs_max
use xc_f03_lib_m, only: xc_f03_func_end, xc_f03_family_from_id, xc_f03_func_get_info, xc_f03_func_info_get_flags, xc_f03_func_t, &
                        xc_f03_functional_get_number, xc_f03_func_init, xc_f03_hyb_exx_coef, XC_FAMILY_GGA, XC_FAMILY_HYB_GGA, &
                        XC_FAMILY_HYB_LDA, XC_FAMILY_HYB_MGGA, XC_FAMILY_LDA, XC_FAMILY_MGGA, XC_FLAGS_NEEDS_LAPLACIAN
use fortran_strings, only: to_upper
use Constants, only: Zero
use Definitions, only: wp, iwp, u6, LibxcInt

implicit none
character(len=*), intent(in) :: Label
real(kind=wp), intent(out) :: ExFac
#include "nq_info.fh"
integer(kind=iwp) :: i, istatus, Lu
character(len=len(Label)) :: UpLabel
character(len=256) :: Line
character(len=80) :: Word1, Word2, Word3
integer(kind=iwp), external :: IsFreeUnit

ExFac = Zero

Lu = IsFreeUnit(11)
call molcas_open(Lu,'FUNCDATA')

UpLabel = to_upper(Label)
do
  ! First find the line that starts with the keyword name
  read(Lu,'(A)',iostat=istatus) Line
  if (istatus /= 0) then
    call WarningMessage(2,' Find_Functional: Undefined functional type!')
    write(u6,*) '         Functional=',trim(Label)
    call Quit_OnUserError()
  end if
  Line = adjustl(Line)
  if ((Line == '') .or. (Line(1:1) == '#')) cycle
  read(Line,*,iostat=istatus) Word1,Word2,Word3
  if (to_upper(Word1) == UpLabel) then
    ! Once found, read the second word
    read(Word2,*,iostat=istatus) nFuncs
    if (istatus == 0) then
      ! If it's a number, read the component functionals and factors
      if (nFuncs > nFuncs_max) then
        call WarningMessage(2,' Find_Functional: Too many components!')
        write(u6,*) '         nFuncs=',nFuncs
        call Quit_OnUserError()
      end if
      do i=1,nFuncs
        read(Lu,*) Coeffs(i),Word2
        func_id(i) = get_func(Word2)
      end do
      read(Word3,*) ExFac
    else
      ! Otherwise, this is just an alias for a Libxc functional
      nFuncs = 1
      func_id(1) = get_func(Word2)
      ExFac = get_exx(func_id(1))
    end if
    ! Assign functional type ("maximum" of its components)
    Functional_Type = -1
    do i=1,nFuncs
      select case (xc_f03_family_from_id(func_id(i)))
        case (XC_FAMILY_LDA,XC_FAMILY_HYB_LDA)
          Functional_type = max(Functional_type,LDA_type)
        case (XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
          Functional_type = max(Functional_type,GGA_type)
        case (XC_FAMILY_MGGA,XC_FAMILY_HYB_MGGA)
          if (needs_laplacian(func_id(i))) then
            Functional_type = max(Functional_type,meta_GGA_type2)
          else
            Functional_type = max(Functional_type,meta_GGA_type1)
          end if
      end select
    end do
    exit
  end if
end do

close(Lu)

contains

function get_func(xcLabel)

  integer(kind=LibxcInt) :: get_func
  character(len=*), intent(in) :: xcLabel

  get_func = xc_f03_functional_get_number(xcLabel)
  if (get_func < 0) then
    call WarningMessage(2,' Find_Functional: Undefined functional in Libxc!')
    write(u6,*) '         Functional=',trim(xcLabel)
    call Quit_OnUserError()
  end if

end function get_func

function get_exx(id)

  real(kind=wp) :: get_exx
  integer(kind=LibxcInt), intent(in) :: id
  type(xc_f03_func_t) :: func

  call xc_f03_func_init(func,id,0_LibxcInt)
  get_exx = xc_f03_hyb_exx_coef(func)
  call xc_f03_func_end(func)

end function get_exx

function needs_laplacian(id)

  logical(kind=iwp) :: needs_laplacian
  integer(kind=LibxcInt), intent(in) :: id
  integer(kind=LibxcInt) :: flags
  type(xc_f03_func_t) :: func

  call xc_f03_func_init(func,id,0_LibxcInt)
  flags = xc_f03_func_info_get_flags(xc_f03_func_get_info(func))
  needs_laplacian = iand(flags,XC_FLAGS_NEEDS_LAPLACIAN) > 0
  call xc_f03_func_end(func)

end function needs_laplacian

end subroutine Find_Functional
